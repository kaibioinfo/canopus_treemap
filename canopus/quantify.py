import sklearn
from sklearn.ensemble import ExtraTreesClassifier
import re
import numpy as np
import pandas as pd
import scipy
from scipy.stats import trim_mean
from .visualization import CanopusRenderer

def subtractBlank(Quant, blank="^.+(?i:blank).+$"):
    bl = re.compile(blank)
    blanks = Quant.loc[:, [m for m in Quant.columns if bl.match(m)]]
    blankEx = 10*blanks.max(axis=1)
    return Quant[~(Quant.max(axis=1)<=blankEx)]

def binnify(Quant):
    S = normalizeByQuantiles(Quant)
    pseudoCount = np.percentile(S.where(S>0).stack().values,1)
    S += pseudoCount
    S = np.round(2*np.log10(S))/2.0
    minimum = S.min().min()
    S -= minimum
    return S

def normalizeByMean(Quant):
    Normed = Quant.copy()
    pseudoCount = Quant.where(Quant>0).min()[1:].min()
    asum = Quant.sum(axis=0)
    for feature in Quant.index:
        s = Quant.loc[feature,:].to_numpy()
        Normed.loc[feature,:] = (s+pseudoCount) / asum
    return Normed

def normalizeByQuantiles(Quant):
    df = Quant.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df

def quantileAndCompoundNormalization(Quant):
    q = normalizeByQuantiles(Quant)
    return q / q.max(axis=0)

def normalizeByLogGeom(Quant):
    Normed = Quant.copy()
    # divide by total ion count
    Normed /= Normed.sum(axis=0)
    pseudoCount = np.percentile(M.where(M>0).stack().values,1)
    for feature in Quant.index:
        s = Quant.loc[feature,:].to_numpy()
        geom = np.mean(np.log10(s[np.nonzero(s)]))
        v = np.zeros(len(s)) + pseudoCount
        s = s + v
        Normed.loc[feature,:] = np.log10(s) - geom
    return Normed

def orderByAbsoluteDifference(Quant, group1, group2):
    group1 = re.compile(group1)
    group2 = re.compile(group2)
    A = [m for m in Quant.columns if group1.match(m)]
    B = [m for m in Quant.columns if group2.match(m)]
    fold_changes = trim_mean(Quant.loc[:,A],0.1,axis=1) - trim_mean(Quant.loc[:,B],0.1,axis=1)
    table = pd.DataFrame(dict(compound=Quant.index, weight=np.abs(fold_changes), difference=fold_changes))
    table.sort_values(by="weight",ascending=False, inplace=True)
    return table

def orderByFoldChange(Quant,group1,group2):
    pseudoCount = np.percentile(Quant.where(Quant>0).stack().values,1)
    group1 = re.compile(group1)
    group2 = re.compile(group2)
    A = [m for m in Quant.columns if group1.match(m)]
    B = [m for m in Quant.columns if group2.match(m)]
    fold_changes = (trim_mean(Quant.loc[:,A],0.1,axis=1)+pseudoCount) / (trim_mean(Quant.loc[:,B],0.1,axis=1)+pseudoCount)
    table = pd.DataFrame(dict(compound=Quant.index, weight=np.abs(np.log10(fold_changes)), fold_change=fold_changes))
    table.sort_values(by="weight",ascending=False, inplace=True)
    return t

def orderByDiscrimination(Quant, group1, group2,ntrees=1000):
    group1 = re.compile(group1)
    group2 = re.compile(group2)
    gf = [m for m in Quant.columns if group1.match(m)]
    wild = [m for m in Quant.columns if group2.match(m)]
    X = np.concatenate([np.stack([Quant.loc[:,s] for s in gf]),np.stack([Quant.loc[:,s] for s in wild])])
    Y = np.concatenate([np.ones(len(gf)),np.zeros(len(wild))])
    f=ExtraTreesClassifier(n_estimators=ntrees)
    f.fit(X,Y)
    bestFeatures = pd.DataFrame(dict(compound=Quant.index, weight=f.feature_importances_))
    bestFeatures.sort_values(by="weight",ascending=False,inplace=True)
    return bestFeatures

import scipy
from scipy.stats import mannwhitneyu
def permutationTest(sirius, table, threshold=0.33):
    compounds = [
        sirius.compounds[compound] for compound in table.index
    ]
    tf={}
    for index in sirius.mapping:
        category = sirius.mapping[index]
        mp = np.array([compound.canopusfp[index]>=threshold for compound in compounds])
        A = table.weight.loc[mp]
        B = table.weight.loc[~mp]
        if (len(A)<10 or len(B)<10):
            continue
        pvalue = mannwhitneyu(A, B, alternative="greater")
        tf[category.name] = pvalue.pvalue
        
    return pd.DataFrame.from_dict(tf,orient="index",columns=["pvalue"]).sort_values(by="pvalue",ascending=True)

import re
functionType = type(lambda x: x)
def differentialAnalysis(sirius, Quant, nameLeft,regexpLeft, nameRight,regexpRight,threshold, restrictTo=".+",method="forest",n=10):
    restrictTo = re.compile(restrictTo)
    subset = [n for n in Quant.columns if restrictTo.match(n)]
    Q = Quant.loc[:,subset]
    R = CanopusRenderer(sirius)
    R.useQuantification(Q)
    R.defineGroup(nameLeft,regexpLeft,"seagreen")
    R.defineGroup(nameRight,regexpRight,"steelblue")
    if type(method) is functionType:
        compounds = method(Q,regexpLeft,regexpRight)
    elif method == "forest":
        compounds = orderByDiscrimination(Q,regexpLeft, regexpRight)
    elif method == "log":
        compounds = orderByFoldChange(Q,regexpLeft, regexpRight)
    elif method == "abs":
        compounds = orderByAbsoluteDifference(Q,regexpLeft, regexpRight)
    else:
        return None
    compounds.weight[compounds.weight<threshold] = 0
    compounds["mz"] = [sirius.compounds[n].mz for n in compounds["compound"]]
    compounds.set_index("compound",drop=True,inplace=True)
    display(compounds).head(n)
    display("%d discriminating compounds" % np.nonzero(compounds.weight)[0].shape[0])
    display(permutationTest(sirius,compounds))
    for row in compounds.head(n).values:
        display("Weight = %f" % row[1])
        R.shortdesc(row[0])