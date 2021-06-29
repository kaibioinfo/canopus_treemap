from IPython.display import HTML, display, Javascript
import uuid
from .ontology import Ontology, SiriusWorkspace,CanopusStatistics,Compound,Formula
import json
import pkg_resources
from pathlib import Path
import pandas as pd
import re
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math


class CanopusRenderer:
  def __init__(self, workspace, uid=None):
    self.workspace = workspace
    self.uid = uid if uid is not None else str(uuid.uuid4().hex)
    self.treemaps = []
    self.use_probabilities = False
    self.groups = []

  def __compound__(self,c):
    if type(c) is not Compound:
      return self.workspace.compounds[c]
    else:
      return c

  def useQuantification(self,quant):
    self.quantification = quant

  def defineGroup(self,name, regexp, groupcolor):
    """add a new discriminative group. All samples matching the pattern are part of this group. They are plotted together and in the same color."""
    self.groups.append((name,re.compile(regexp),groupcolor))

  def shortdesc(self,compound,form=None, threshold=0.33, filter=None):
    compound = self.__compound__(compound)
    # print formula and adduct
    display(HTML("<h3>%s (%s)</h3><br />Zodiac Score: %d<br />Feature-ID: %s" % (
      Formula(compound.formula).to_html(), compound.adduct, int(compound.zodiacScore*100) if (compound.zodiacScore is not None and not math.isnan(compound.zodiacScore) )else 0,
      compound.name
    )))
    # plot top structure results
    f = Path("%s/structure_candidates.tsv" % compound.directory)
    if f.exists():
      table = pd.read_csv(f,sep="\t")
      display(table.sort_values(by="score",ascending=False).head(10))
    # plot classification
    self.canopusTreeTable(compound)
    # plot quantification
    self.quantplot(compound, filter)

  def quantplot(self,compound,filter=None,Quant=None):
    if Quant is None:
      Quant = self.quantification
    if filter:
      filter = re.compile(filter)
      m = [x for x in Quant.columns if filter.match(x)]
      Quant = Quant.loc[:,m] 
    compound = self.__compound__(compound)
    fid = compound.name
    grp = self.groups if self.groups else [('all',re.compile('.*'),'steelblue')]
    n=0
    mids=[]
    names = []
    for (name,pattern,color) in grp:
      members = [column for column in Quant.columns if pattern.match(column)]
      indizes = np.arange(n,n+len(members))
      mid = n+len(indizes)//2
      n += len(members)
      plt.bar(
        indizes, Quant.loc[fid,members], width=1,color=color
      )
      mids.append(mid)
      names.append(name)
    plt.xticks(ticks=mids,labels=names)

  def canopusTreeTable(self,compound,threshold=0.25, npc=False):
    compound = self.__compound__(compound)
    fp = self.workspace.statistics.npcCategoriesAndProbabilitiesFor(compound,threshold) if npc else self.workspace.statistics.categoriesAndProbabilitiesFor(compound,threshold)
    klasses = set(fp.keys())
    astree = self.workspace.npc_ontology.root if npc else self.workspace.ontology.root
    # add root
    fp[astree] = 1.0
    display(HTML(self.__getTreeAsHTML__(fp,astree,0)))

  def __tooltip__(self,text, tipp):
    return "<span title='" + tipp + "'>" + text + "</span>"

  def __getTreeAsHTML__(self,fp,root,indent):
    s=""
    weight = "bold" if fp[root]>0.5 else "italic"
    if root.name:
      s += "<li>" + self.__tooltip__("<span style='font-weight:" + weight +"'>" + root.name + "</span>" + "<span style='margin-left:15pt;font-weight:" + weight + "'>" + 
      (str(int(fp[root]*100)) + " %" if fp[root] > 0 else "") + "</span>", root.description)
    t = "<ul>"
    c=0
    for child in root.children:
      if child in fp:
        t += self.__getTreeAsHTML__(fp,child,indent+1)
        c +=1
    t += "</ul>"
    if c>0:
      s += t
    if root.name:
      s += "</li>"
    return s
    
  def addTreemap(self,statistics = None):
      if statistics is None:
          statistics = self.workspace.statistics
      self.treemaps.append(statistics)
      
  
  def renderCSS(self):
    display(HTML(self.getCSS()))

  def getCSS(self):      
      return r"""
<style>
#main {
  width: 100%;
}
.sidebar {
  position: absolute;
  margin-right:20px;
  top: 10%;
  bottom:50%;
  right: 0%;
  width:200px;
  height:300px;
  padding: 12px;
  margin-left:1px; margin-right:1px;
  /*background-color: steelblue;*/
  display:inline-block;
  border-radius: 5px;
  border: 2px solid #000000;
  min-width: 100px;
  min-height: 100px;
}
.description {
  font-weight: 600;
  fill: #fff;
  font-size: 14px;
}

.sequence {
  width: 600px;
  height: 70px;
}


.node {
  font-weight: 900;
  color: #fff;
}
div.left {
  float:left;
}
div.right {
  float:right;
}

.chart path, .chart path {
  stroke: #fff;
}
div.chart {
position:relative;
}

h1 {
  text-align: center;
}

.chart > .explanation {
  position: absolute;
  top: 280px;
  left: 250px;
  width: 180px;
  text-align: center;
  z-index: -1;
}


.percentage {
  font-size: 2.5em;
}

.catemph {
  font-weight: bolder;
  font-size:125%;
}

.treelegend {
  min-height:55px;
}
div.node {
  padding: 3px;
  margin-left:1px; margin-right:1px;
  display:inline-block;
  border-radius: 5px;
  border: 2px solid #000000;
}
.sidebar div {
  font-weight: bold;
  fill: #fff;
}
</style>
"""

  def toFile(self, filename):
      with open(filename,"w") as fhandle:
        fhandle.write("<html>\n")
        fhandle.write(self.getCSS())
        fhandle.write("<script>\n")
        fhandle.write(self.getJavascript())
        fhandle.write("</script>\n")
        fhandle.write("<body>\n")
        fhandle.write(self.getHTML())
        fhandle.write("</body>\n")
        fhandle.write("</html>\n")

  def render(self):
      self.renderCSS()
      self.renderHTML()
      self.renderJavascript()

  def renderHTML(self):
    display(HTML(self.getHTML()))

  def getHTML(self):
    specialClasses = ["" for _ in self.treemaps]
    if len(self.treemaps)==2:
      specialClasses[0] = " left"
      specialClasses[1] = " right"

    charts = "\n".join(["<div class=\"chart" + specialClasses[i] + "\" id=\"chart" + self.uid + "_" + str(i) + "\">"
          """
          
              <div class="explanation" style="visibility: hidden;">
                <span class="percentage"></span><br/>
                <span class="category_name"></span>
              </div>
            </div>
          """ for (i,_) in enumerate(self.treemaps)])
    return "<div id=\"main" + "_" +self.uid + "\">" + """
    <div class="treelegend"></div>
    <div class="sequence"></div>
""" + charts + """<div class="sidebar">
    <div>
    <p class="description"></p>
    </div>
</div>"""
  
  def renderJavascript(self):
    display(Javascript(self.getJavascript()))

  def getJavascript(self):
      customCode = "var json = " + json.dumps(self.workspace.json_treemap(use_probabilities=self.use_probabilities)) + """;
var nodes = [];
allnodes(json, nodes);
loadColors(nodes);
      """;
      charts = "\n".join(["createVisualization(" + json.dumps(self.workspace.json_treemap(x,use_probabilities=self.use_probabilities)) + ", \"" + self.uid + "\""  + ",\"" + self.uid + "_" + str(i) + "\");" for (i, x) in enumerate(self.treemaps)])
      json_code = pkg_resources.resource_string("canopus.resources", "treemap.js").decode("utf-8")
      return json_code.replace('"<CUSTOM-CODE>";', customCode + "\n" + charts)
      