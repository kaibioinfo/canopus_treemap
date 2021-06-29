from .ontology import SiriusWorkspace,Compound,Formula,Category
from .visualization import CanopusRenderer
from .network import MolecularNetwork
from .quantify import quantileAndCompoundNormalization,binnify,subtractBlank,permutationTest
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import RandomForestClassifier 
from IPython.display import HTML, display, Javascript
from scipy.stats import trim_mean
import seaborn as sbn
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import math
import glob
from pathlib import Path
import re
import copy

class Canopus(object):

  def __init__(self, sirius, gnps=None):
    self.__cleanNorm__()
    self.probabilityThreshold=0.33
    self.sirius = SiriusWorkspace(sirius)
    self.Quant = self.sirius.make_quant()
    self.conditions = dict()
    if gnps:
      networkfile = list(Path(gnps,"gnps_molecular_network_graphml").glob("*.graphml" ))
      if not networkfile:
        print("Cannot find gnps_molecular_network_graphml file in the given directory")
      else:
        self.gnps = MolecularNetwork.parse(networkfile[0])
        self.gnps.feedSirius(self.sirius)
      identificationFile = list(Path(gnps, "DB_result").glob("*.tsv"))
      if identificationFile:
        self.gnps_hits = pd.read_csv(identificationFile[0],sep="\t")
      else:
        self.gnps_hits = None
      clusterinfo = list(Path(gnps, "clusterinfo_summary").glob("*.tsv"))
      if clusterinfo:
        self.gnps.feedClusterInfo(clusterinfo[0])
    self.all = Condition("all", ".*", "black")
    self.all.samples = [c for c in self.Quant.columns]
    self.all.compounds = list(self.sirius.compounds)
    self.blank = None
    self.__keywords__ = None

  def __getitem__(self, key):
    """
    lookup the given key in the conditions or categories or compounds. Use this in combination with Jupyter autocompletion
    """
    if key in self.conditions:
      return self.condition(key)
    if key in self.sirius.ontology.categoriesByName:
      return self.sirius.ontology.categoriesByName[key]
    if key in self.sirius.compounds:
      return self.sirius.compounds[key]
    raise ValueError("Unknown keyword '" + key + "'")

  def _ipython_key_completions_(self):
    return self.__all_keywords__()

  def exportCytoscape(self, file, selection=None):
    """
    export the data into a cytoscape graphml file.
    file      -- name of the output file
    selection -- a list/set of class names which should be used for annotating nodes in the network.
                 if omitted, CANOPUS will try to find useful classes itself.  
    """
    self.gnps.write(file, selection)

  def defineCondition(self,name, regexp=None, color=None):
    """
    define a condition (group, organ, treatment, ...)
    name    -- the name of the condition. It will be used as label in plots.
    regexp  -- a regular expression string to find samples belonging to this condition. 
              If left None, all samples that contain the name of the condition in their 
              filename will be used instead.
    color   -- the color used in plots for this condition
    """
    self.conditions[name] = Condition(name,regexp,color)
    self.__conditionsHaveChanged__()
    self.__assignSamples__(self.condition(name))
    return self.conditions[name]

  def defineBlank(self,regexp=".*(?i:blank).*"):
    """
    specify that there exists blank run samples. The given regular expression string is
    used to identify the blank run samples. Otherwise, all filenames containing
    the substring "blank" are used as blank run
    """
    self.blank = Condition("blank", regexp)
    self.__cleanNorm__()
    self.__assignSamples__(self.blank)

  def condition(self, name, color=None):
    """
    return the condition with the given name
    """
    c=None
    if type(name) is Condition:
      c=name
    elif name is None:
      c = self.all
    else: 
      c = self.conditions[name]
    if color:
      c = copy.copy(c)
      c.color = color
    return c

  def invert(self, condition, color, name=None):
    """invert a condition"""
    condition = self.condition(condition)
    x = Condition(name if name else "not " + condition.name, 
      "(?!:"+condition.regexpString+")",
      color)
    self.__assignSamples__(x)
    return x

  def join(self, *conditions,name=None,color=None):
    """
    Select all samples that are part of at least one of the given conditions
    """
    if name is None:
      name = "|".join([self.condition(c).name for c in conditions])
    regexp = "(?:" + "|".join([self.condition(c).regexpString for c in conditions]) + ")"
    if color is None:
      for c in conditions:
        d=self.condition(c).coloring()
        if d:
          color = d
          break
    joined = Condition(name,regexp,color)
    self.__assignSamples__(joined)
    return joined

  def select(self,*conditions,name=None,color=None):
    """
    Select all samples that are part of ALL the given conditions
    """
    if name is None:
      name = "&".join([self.condition(c).name for c in conditions])
    # we cannot define a proper regexp =/
    aset = set(self.Quant.columns)
    for c in conditions:
      aset = aset & set(c.samples)
    samples = list(aset)
    pseudoReg = "^(?:" + "|".join(samples) + ")$"
    if color is None:
      for c in conditions:
        d=self.condition(c).coloring()
        if d:
          color = d
          break
    condition = Condition(name,pseudoReg, color)
    self.__assignSamples__(condition)
    return condition

  def differential(self, conditionLeft, conditionRight=None, method="robust_forest",thresholding=True, binning=False):
    """
    Start a differential expression analysis: Find the compounds that differentiate between the two
    conditions conditionLeft and conditionRight. Order all compounds by how strongly they differentiate
    between both conditions. Order all compound categories according to how strongly they differentiate between
    both conditions: A category differentiates, if there are more highly differentiating compounds than it would
    be expected from this category (using a whitney-u permutation test on the ranked compound list).
    conditionLeft -- first condition
    conditionRight -- either the second condition, or, when left None, all other samples that do not belong to conditionLeft
    method -- the method to use, one of:
    'forest' -- order compounds by a random forest, that differentiates between the two conditions
    'robust_forest' -- the same, but use logarithmized and rounded intensities, such that the forest cannot
                       overfit on small intensity deviations
    'fold_change' -- order compounds by the average fold change between the conditions
    'highest_expression' -- order compounds by fold change, but only with respect to compounds which are higher
                            expressed in the left condition than in the right condition
    """
    return DifferentialApi(self, self.condition(conditionLeft), self.condition(conditionRight) if conditionRight else self.invert(conditionLeft, "red"), method,thresholding,binning=binning)

  def treemap(self, condition=None, method="binary",probabilities=True):
    """
    draw a sunburst tree map. Either of all compounds in all samples, OR for samples following the specific condition
    condition   -- if given, the condition for which we want to display the compound classes
    method      -- one of the following keywords:
      "binary"    -- select all compounds which are expressed with at least 10% relative intensity in at least one sample
      "quantify"  -- use the Quant table to display the relative amount of compound classes in intensity 
    probabilities -- if True (default), use expected class counts, otherwise decide by a 50% threshold
    """
    self.__ensureNormalization__()
    condition = self.condition(condition)
    r = CanopusRenderer(self.sirius)
    r.use_probabilities = probabilities
    stats = self.sirius.statistics
    if method == "binary":
      stats = self.sirius.select(condition.compounds)
    elif method == "quantify":
      r.useQuantification(self.QuantNormalized)
    else:
      raise ValueError("Unknown keyword '" + method + "'") 
    r.addTreemap(stats)
    r.render()

  def differentialTreemap(self, conditionLeft, conditionRight, method="binary",probabilities=True):
    """
    draw two sunburst tree maps, one with conditionLeft and one with conditionRight
    method      -- one of the following keywords:
      "binary"    -- select all compounds which are expressed with at least 10% relative intensity in at least one sample
      "quantify"  -- use the Quant table to display the relative amount of compound classes in intensity 
    probabilities -- if True (default), use expected class counts, otherwise decide by a 50% threshold
    """
    self.__ensureNormalization__()
    conditionLeft = self.condition(conditionLeft)
    conditionRight = self.condition(conditionRight)
    r = CanopusRenderer(self.sirius)
    r.use_probabilities = probabilities
    statsLeft = self.sirius.statistics
    statsRight = self.sirius.statistics
    if method == "binary":
      statsLeft = self.sirius.select(conditionLeft.compounds)
      statsRight = self.sirius.select(conditionRight.compounds)
    elif method == "quantify":
      raise ValueError("Not implemented yet!")
    else:
      raise ValueError("Unknown keyword '" + method + "'") 
    r.addTreemap(statsLeft)
    r.addTreemap(statsRight)
    r.render()

  def treemapFromTwoDatasets(self, otherDataset, method="binary",probabilities=True):
    raise ValueError("Not implemented yet!")

  def describe(self, compound, conditions=None):
    """
    Plot several statistics/information about the given compound.
    compound -- id of the compound
    conditions -- optional, a list of conditions which should be plotted in a histogram 
    """
    compound = self.__fid__(compound)
    self.featureHeader(compound)
    self.gnpsHit(compound)
    self.identification(compound)
    self.classification(compound)
    self.histogram(compound,conditions)

  def molecularNetwork(self):
    self.gnps.render()

  def gnpsHit(self, compound=None, category=None):
    """
    return the GNPS identification table (with all database hits of the dataset)
    compound -- if given, only return identifications for this compound
    category -- if given, only return identifications for compounds of this category
    """
    t = self.gnps_hits
    if compound:
      compound = self.__fid__(compound)
    if category:
      t = t[t["#Scan#"].isin(set(self.allFromCategory(category)))]
    if compound:
      t = t[t["#Scan#"]==int(compound)]
    display(t)
    return t

  def featureHeader(self,compound):
    """
    display the feature id and formula of the given compound, as well as its ZODIAC score
    """
    compound = self.__compound__(compound)
    display(HTML("<h3>%s (%s)</h3><br />Zodiac Score: %d<br />Feature-ID: %s" % (
      Formula(compound.formula).to_html(), compound.adduct, int(compound.zodiacScore*100) if (compound.zodiacScore is not None and not math.isnan(compound.zodiacScore) )else 0,
      compound.name
    )))

  def classification(self,compound):
    """
    display the compound categories which are predicted for this compound
    """
    compound = self.__fid__(compound)
    R = CanopusRenderer(self.sirius)
    R.canopusTreeTable(compound)

  def npcClassification(self, compound):
    """
    display the natural product categories (NPC) which are predicted for this compound
    """
    compound = self.__fid__(compound)
    R = CanopusRenderer(self.sirius)
    R.canopusTreeTable(compound, npc=True)


  def npcSummary(self, condition=None):
    """
    outputs a Pandas dataframe containing all NPC classes
    """
    compounds = [self.sirius.compounds[x] for x in self.condition(condition).compounds] if condition else self.sirius.compounds.values()
    rows = []
    for compound in compounds:
      vec = self.sirius.statistics.npcPrimaryAssignmentVector(compound, self.probabilityThreshold)
      probs = self.sirius.statistics.npcCategoriesAndProbabilitiesFor(compound, self.probabilityThreshold)
      ary = [compound.name,str(compound.directory)]
      for name in vec:
        if name is None:
          ary.append("N/A")
        else:
          ary.append(name.name)
        ary.append(probs[name] if name in probs else 0.0)
      key=self.sirius.statistics.assignments[compound]
      ary.append(key.name)
      ary.append(compound.canopusfp[self.sirius.statistics.workspace.revmap[key]])
      rows.append(ary)
    return pd.DataFrame(data=rows, columns=["name", "directoryName", "pathway", 
                                            "pathwayProbability", "superclass", "superclassProbability", 
                                            "class", "classProbability", "ClassyFirePrediction", "ClassyFirePredictionProbability"]).set_index("name", drop=True)

  def identification(self,compound):
    """
    display the CSI identifications for this compound
    """
    compound = self.__compound__(compound)
    f = Path("%s/structure_candidates.tsv" % compound.directory)
    if f.exists():
      table = pd.read_csv(f,sep="\t")
      display(table.sort_values(by="CSI:FingerIDScore",ascending=False).head(10))

  def heatmap(self, category, conditions=None, logarithmic=False):
    """
    display a heatmap with the intensities of all compounds belonging to the given category
    category -- name of the category we are interested in
    conditions -- list of conditions we want to compare. By default, compare all conditions
    logarithmic -- if True, use the logarithmic scale
    """
    if conditions is None:
      conditions = list(self.conditions.values())
    else:
      conditions = [self.condition(c) for c in conditions]
    cmps=self.allFromCategory(category)
    Q=self.QuantLogarithmic if logarithmic else self.QuantNormalized
    s=Q.loc[cmps,:]
    m = s.min().min()
    m2 = s.max().max()
    for (i,grp) in enumerate(conditions):
      plt.subplot(1,len(conditions),i+1)
      cb=plt.imshow(Q.loc[cmps,grp.samples],vmin=m,vmax=m2,aspect="auto")
      plt.title(grp.name)
    

  def histogram(self, compound, conditions=None):
    """
    display a bar plot of the compounds intensities among the given conditions
    compound -- feature id of the compound
    conditions -- if given, plot only the intensities for the given intensities
    """
    compound = self.__fid__(compound)
    if conditions is None:
      conditions = list(self.conditions.values())
    else:
      conditions = [self.condition(c) for c in conditions]
    plt.subplot(1, 2, 1)
    self.__drawQuant__(self.QuantNormalized,compound,conditions)
    plt.ylabel("normalized intensity")
    plt.legend()
    plt.subplot(1, 2, 2)
    self.__drawQuant__(self.QuantLogarithmic,compound,conditions)
    plt.ylabel("logarithmized intensity")
    

  def __drawQuant__(self, quant, compound, conditions):
    compound = self.__fid__(compound)
    offset = 0
    ticks=[]
    pos=[]
    for condition in conditions:
      plt.bar(np.arange(offset, offset+len(condition.samples)), quant.loc[compound,condition.samples],color=condition.coloring(),width=1,label=condition.name)
      if len(condition.samples) < 10:
        for i,s in enumerate(condition.samples):
          pos.append(offset+i)
          ticks.append(s)
      offset += len(condition.samples)
    #plt.xticks(pos,ticks,rotation='vertical')

  def __compound__(self,c):
    if type(c) is Compound:
      return c
    else: 
      return self.sirius.compounds[self.__fid__(c)]

  def allFromCategory(self, category, threshold=None):
    if threshold is None:
      threshold = self.probabilityThreshold
    i = self.__findCategoryIndex__(category)
    return [c for c in self.QuantNormalized.index if self.sirius.compounds[c].canopusfp[i] >= threshold]

  def __findCategoryIndex__(self,category):
    category = self.__cat__(category)
    if category in self.sirius.revmap:
      return self.sirius.revmap[category]
    raise ValueError("Unknown category '" + category + "'")

  def __assignSamples__(self,condition):
    self.__ensureNormalization__()  
    condition.samples = [index for index in self.Quant.columns if condition.regexp.match(index)]
    condition.compounds = []
    self.__keywords__ = None
    n = self.QuantNormalized.max(axis=1)
    for compound in self.compounds:
      if np.max(self.QuantNormalized.loc[compound,condition.samples]/n[compound])>=0.05:
        condition.compounds.append(compound)

  def __normalize__(self):
    if self.blank is None:
      self.defineBlank()
    q = subtractBlank(self.Quant, self.blank.regexp)
    self.compounds = q.index
    self.QuantNormalized = quantileAndCompoundNormalization(q)
    self.QuantLogarithmic = binnify(q)

  def __conditionsHaveChanged__(self):
    colors = sbn.color_palette(n_colors=len(self.conditions))
    for (i,c) in enumerate(self.conditions):
      self.conditions[c].fallbackColor = colors[i]

  def __ensureNormalization__(self):
    if self.QuantNormalized is None:
      self.__normalize__()

  def __all_keywords__(self):
    if self.__keywords__ is None:
      self.__keywords__ = list(self.Quant.index)
      self.__keywords__.extend(list(self.sirius.ontology.categoriesByName))
      self.__keywords__.extend(list(self.conditions))
    return self.__keywords__

  def __cleanNorm__(self):
      self.QuantNormalized = None
      self.QuantLogarithmic = None

  def __cat__(self,key):
    if type(key) is Category:
      return key
    if type(key) is str:
      if key in self.sirius.ontology.categories:
        return self.sirius.ontology.categories[key]
      elif key in self.sirius.ontology.categoriesByName:
        return self.sirius.ontology.categoriesByName[key]
    raise ValueError("Unknown category '" + key + "'")

  def __fid__(self,fid):
    if type(fid) is str:
      return fid
    if type(fid) is int:
      return str(fid)
    if type(fid) is Compound:
      return compound.name
    return str(fid)

class Condition(object):

  def __init__(self, name, regexp=None, color=None):
    self.name = name
    self.regexpString = regexp if regexp else ".*"+name+".*"
    self.regexp = re.compile(self.regexpString)
    self.color = color
    self.fallbackColor = None

  def coloring(self):
    return self.color if self.color else self.fallbackColor

class DifferentialApi(object):

  def __init__(self,canopus, conditionLeft, conditionRight, method, thresholding,binning=False):
    self.conditionLeft = conditionLeft
    self.conditionRight = conditionRight
    self.canopus = canopus
    self.binning=binning
    self.usedMethod = method
    self.orderCompounds(method)
    if thresholding:
        if type(thresholding) is bool:
            self.threshold()
        else:
            self.threshold(thresholding)

  def threshold(self, userdefined=None):
    if userdefined:
      self.ordering.loc[self.ordering[self.ordering.weight<=userdefined].index, "weight"]=0.0
    else:
      if self.usedMethod == "fold_change" or self.usedMethod == "highest_expression":
        self.ordering.loc[self.ordering[self.ordering.weight<1].index, "weight"]=0.0
      else:
        cutoff = np.minimum(10*self.ordering.weight.median(), np.percentile(self.ordering.weight.values, 90))
        self.ordering.loc[self.ordering[self.ordering.weight<=cutoff].index, "weight"]=0.0

  def _ipython_display_(self):
    display(HTML("<h5>Top differentiating compounds</h5>"))
    self.topCompounds()
    display(HTML("<h5>Top differentiating categories</h5>"))
    self.topCategories()
    display(HTML("<h5>Category heatmaps</h6>"))
    t = permutationTest(self.canopus.sirius,self.ordering,self.canopus.probabilityThreshold)
    for row in t.head(5).index:
      display(HTML("<h6>" + row + "</h6>"))
      self.canopus.heatmap(row, [self.conditionLeft, self.conditionRight])
      plt.figure()

  def topCompounds(self,n=None,category=None):
    """
    display the top differentiating compounds
    n -- if given, only display the top n compounds. By default displays 20 compounds
    """
    t = self.ordering
    if category:
      cmps = set(self.canopus.allFromCategory(category))
      t = t.query("compound in @cmps")
    s=t.head(n if n else 20)
    display(s)
    return t

  def topCategories(self,n=20):
    """
    display the top differentiating categories
    n -- if given, only display the top n categories. By default displays 20 categories
    """
    t = permutationTest(self.canopus.sirius,self.ordering,self.canopus.probabilityThreshold)
    display(t.head(n if n else 20))
    return t.head(n) if n else t

  def orderCompounds(self,method):
    """
    reorder the compounds by the given method, which can be one of the following:
    'forest' -- order compounds by a random forest, that differentiates between the two conditions
    'robust_forest' -- the same, but use logarithmized and rounded intensities, such that the forest cannot
                       overfit on small intensity deviations
    'fold_change' -- order compounds by the average fold change between the conditions
    'highest_expression' -- order compounds by fold change, but only with respect to compounds which are higher
                            expressed in the left condition than in the right condition
    """
    self.method = method
    if method == "robust_forest":
      self.orderByRobustForest()
    elif method == "forest":
      self.orderByForest()
    elif method == "fold_change":
      self.orderByFoldChange()
    elif method == "highest_expression":
      self.orderByFoldChange(False)
    else:
      raise ValueError("method '" + self.method + "' is unknown. It has to be one of the following: 'robust_forest', 'forest', 'fold_change', 'highest_expression")

  def orderByForest(self):
    self.__forest_ordering__(self.canopus.QuantNormalized)

  def orderByRobustForest(self):
    self.__forest_ordering__(self.canopus.QuantLogarithmic)

  def orderByFoldChange(self,bidirection=True):

    Quant = Quant = self.canopus.QuantNormalized
    pseudoCount = np.percentile(Quant.where(Quant>0).stack().values,0.1)
    A = self.conditionLeft.samples
    B = self.conditionRight.samples
    fold_changes = (trim_mean(Quant.loc[:,A],0.05,axis=1)+pseudoCount) / (trim_mean(Quant.loc[:,B],0.05,axis=1)+pseudoCount)
    #fold_changes = (np.mean(Quant.loc[:,A],axis=1)+pseudoCount) / (np.mean(Quant.loc[:,B],axis=1)+pseudoCount)

    W = np.log10(fold_changes)

    if self.binning:
      binsize = 0.5 if self.binning==True else self.binning
      scale = 1.0/binsize
      W = np.round(W*scale)/scale
    if bidirection:
      W = np.abs(W)
    table = pd.DataFrame(dict(compound=Quant.index, weight=W, fold_change=fold_changes, category=self.__assign_specific_class__(Quant.index)))
    table.sort_values(by="weight",ascending=False, inplace=True)
    table.set_index("compound",drop=True,inplace=True)
    #table[table.weight < 1] = 0.0 # we do not trust the lower values anyways
    self.ordering = table

  def __forest_ordering__(self,Quant):
    X = np.concatenate([Quant.loc[:,self.conditionLeft.samples], Quant.loc[:,self.conditionRight.samples]],axis=1).transpose()
    Y = np.concatenate([np.ones(len(self.conditionLeft.samples)),np.zeros(len(self.conditionRight.samples))])
    f=ExtraTreesClassifier(n_estimators=1000)
    f.fit(X,Y)
    bestFeatures = pd.DataFrame(dict(compound=Quant.index, weight=f.feature_importances_, category=self.__assign_specific_class__(Quant.index)))
    bestFeatures.sort_values(by="weight",ascending=False,inplace=True)
    bestFeatures.set_index("compound",drop=True,inplace=True)
    cutoff = np.minimum(10*bestFeatures.weight.median(), np.percentile(bestFeatures.weight.values, 90))
    bestFeatures[bestFeatures.weight < cutoff] = 0.0 # we do not trust the lower values anyways
    self.ordering =  bestFeatures

  def __assign_specific_class__(self, compounds):
    return [self.canopus.sirius.statistics.assignments[self.canopus.__compound__(c)] for c in compounds]