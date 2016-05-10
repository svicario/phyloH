#GoodTest
#Setting Test
from Bio import Phylo
from pandas import DataFrame, MultiIndex
from numpy import random,array
import scipy
from collections import OrderedDict
import sys
sys.path.append('../')
#esecutorePhyloHPandas = imp.load_source("esecutorePhyloHPandas","../esecutorePhyloHPandas.py")
from esecutorePhyloHPandas import *
from lib_script.DiversityMatrixTreeKL import *

def CreatingDataset(nsample=4):
    N=1000000
    com={"-f":"Echinodermata.tree","-g":"GroupTest","-s":"sampleTest"}
    #Creating DataSet
    t=Phylo.read(com["-f"], "newick")
    otu=[x.name for x in t.get_terminals()]
    #countsA=random.random(len(otu))
    #countsA=countsA/sum(countsA)
    #countsB=random.random(len(otu))
    #countsB=countsB/sum(countsB)
    FullCounts=scipy.stats.dirichlet.rvs([1]*len(otu),size=nsample)
    freqE=scipy.stats.dirichlet.rvs([1]*nsample,size=1)
    FullCounts=FullCounts.transpose()*freqE[0]
    #freqE=random.random(1)[0]
    #freqE=array([1-freqE,freqE])
    #FullCounts=array([countsA,countsB]).transpose()*freqE
    #TT=FullCounts.round(3)
    #TT=FullCounts.applymap(lambda x: round(x,ndigits=3))
    #N=round(1/(TT.min())+10,ndigits=-2)
    FullCounts=(FullCounts*N).astype("int")
    #countsA1=FullCounts[:,0]
    #countsA2=FullCounts[:,1]
    #countsB1=FullCounts[:,2]
    #countsB2=FullCounts[:,3]
    freqE=sum(FullCounts)
    grouplevels=["AA","BB"]*(nsample/2)
    grouplevels.sort()
    #grouplevels=["AA","AA","BB","BB"]
    samplelevels=[y[0]+str(x+1) for x,y in zip(range(nsample/2)*2,grouplevels)]
    
    handle=open(com["-s"],"w")
    for i,name in enumerate(samplelevels):
        handle.write("\n".join(["\t".join([name,str(x[0]),x[1]]) for x in zip(FullCounts[:,i],otu)]))
        handle.write("\n")
    handle.close()
    handle=open(com["-g"],"w")
    handle.write("\n".join([x+"\t"+y for x,y in zip(samplelevels,grouplevels)]))
    handle.close()
    
    FullCounts=DataFrame(FullCounts, index=otu, columns=samplelevels)
    #grouplevels=["AA","AA","BB","BB"]
    #samplelevels=["A1","A2","B1","B2"]
    FullCounts.columns=MultiIndex.from_tuples(list(zip(grouplevels,samplelevels)),names=["Group","Sample"])
    db=DBdata()
    db.readTreePandas(com['-f'])
    db.readSampleTable(com["-s"])
    db.readGroupTable(com["-g"])
    return FullCounts, t, otu, db
def traversing(c,count=[0]):
     """
     Putting Names on internal nodes, starting on root with 'L0'
     """
     if not c.is_terminal():
         c.name="L"+str(count[0])
     count[0]+=1
     for C in c.clades:
         C=traversing(C,count)
     return c

def traversingForDescendant(c, desc,anc=None):
     """
     getting from tree, name of  node, if is terminal, branch length, and list of descedant 
     """
     desc.append([c.name,anc,c.is_terminal(),[x.name for x in c.get_terminals()]])
     for C in c.clades:
         C=traversingForDescendant(C,desc,c.name)
     return c



def GetFromFiles(com):
    db=DBdata()
    db.readTreePandas(com['-f'])
    db.readSampleTable(com["-s"])
    db.readGroupTable(com["-g"])
    A=DataFrame(["\t".join(x) for x in db.itemTable])
    AA=A.iloc[:,0].value_counts()
    AAindex=[x.split("\t") for x in AA.index]
    z=MultiIndex.from_tuples(AAindex, names=["Taxon","Sample","Group"])
    AA.index=z
    D=AA.unstack(level=0).transpose()
    ch={}
    ch.update(zip(map(str,range(len(db.samplesNames))),db.samplesNames))
    D.rename(columns=ch,inplace=True)
    #T=self.compressTable()
    Depths, desc, L=db.TreeSummary
    D.fillna(value=0,inplace=True)
    zz=set(D.index.tolist())
    for i in zz.symmetric_difference([x.name for x in db.tree.get_terminals()]):
        D.loc[i]=[0]*D.shape[1]
    return D, db.tree, D.index.tolist(), db

def UltraTreeTest(FullCounts, t, otu, Equal=False):
    #freqS=FullCounts.sum(axis=0)/FullCounts.values.sum()
    if not Equal:
        freqS=FullCounts.sum(axis=0)/FullCounts.values.sum()
        #freqE=FullCounts.sum(axis=1, level="Group").sum(axis=0)/FullCounts.values.sum()
        FullCounts=FullCounts/FullCounts.sum(axis=0)
    else:
        FullCounts=FullCounts/FullCounts.sum(axis=0)
        freqS=FullCounts.sum(axis=0)/FullCounts.values.sum()
        #print freqS
    
    traversing(t.clade,count=[0])
    desc=[]
    traversingForDescendant(t.clade,desc)
    Counts=[x[0:3]+list(FullCounts.loc[x[3]].values.sum(axis=0)) for x in desc]
    #Counts=[x[0:3]+[A.loc[x[3]].values.sum(),B.loc[x[3]].values.sum()] for x in desc]
    Counts=DataFrame.from_records(Counts,columns=["Name","anc","Is_Leaf"]+FullCounts.columns.tolist())
    Counts=Counts.set_index("Name")
    #print Counts.columns
    #Counts=DataFrame({"Name":Name, "Is_leaf":Is_leaf, "anc":anc, "freqA":freqA,"freqB":freqB}, index=Name)
    
    Depths=t.clade.depths()
    D=[[x[0].name,x[1]] for x in Depths.items()]
    Depths={}
    Depths.update(D)
    D=[[x[0],Depths.setdefault(Counts.loc[x[0]].anc,0),x[1]] for x in Depths.items()]
    D.sort(key=lambda x: x[-1])
    
    from scipy.stats import entropy
    OldDepthNode=0
    Halphas=[]
    Hgammas=[]
    WindowSizes=[]
    tots=[]
    freqE=freqS.sum(axis=0, level="Group")
    #print freqE 
    for y in D:
        name, DepthAnc, DepthNode=y
        s=[x[0] for x in D if (min(x[-1], DepthNode)-max(x[1], OldDepthNode))>0]
        if s:
            WindowSizes.append(DepthNode-OldDepthNode)
            OldDepthNode=DepthNode
            temp=Counts.loc[s,FullCounts.columns.tolist()]
            temp.columns=FullCounts.columns
            tots.append(temp.sum())
            tempa=(temp*freqS).sum(axis=1,level="Group")
            Ha=[entropy(tempa.loc[:,x]) for x in tempa]
            #Ha=[entropy(temp.sum(axis=1,level="Group").iloc[:,x].values) for x in range(2)]
            Halphas.append((Ha*freqE).sum())
            Hgammas.append(entropy((temp*freqS).values.sum(axis=1)))
    
    WindowSizes=array(WindowSizes)
    Halpha=sum(array(Halphas)*WindowSizes)/sum(WindowSizes)
    Hgamma=sum(array(Hgammas)*WindowSizes)/sum(WindowSizes)
    Hbeta=Hgamma-Halpha
    HE=entropy(freqE)
    return Halpha,Hgamma,Hbeta,HE,tots




def Report(FullCounts, t, otu,db):
    Halpha,Hgamma,Hbeta,HE,tots=UltraTreeTest(FullCounts, t, otu)
    H=db.GetEntropiesPandas(q="1", Pairwise=1, EqualEffort=0)
    result=DataFrame.from_items([["Hgamma",[Hgamma,H["Hgamma"]]],
        ["Halpha",[Halpha,H["HalphaByEnvironment"]]],
         ["Hbeta",[Hbeta,H["MI_treeAndEnvironment"]]],              
        ["DistTurnover",[Hbeta/HE,H["DistTurnover"].iloc[1,0]]],
        ["DistTurnoverbySample",[Hbeta/HE,H["DistTurnoverBySample"].iloc[1,0]]]],
        columns=["Test","RegularRoutine"],
        orient="index")
    result["Dif"]=result.Test-result.RegularRoutine
    print(result)
    
    H=db.GetEntropiesPandas(q="1", Pairwise=1, EqualEffort=1)
    #Halpha,Hgamma,Hbeta,HE=UltraTreeTest(countsA,countsB,[0.5,0.5])
    Halpha,Hgamma,Hbeta,HE,tots=UltraTreeTest(FullCounts, t, otu, Equal=True)
    result=DataFrame.from_items([["Hgamma",[Hgamma,H["Hgamma"]]],
        ["Halpha",[Halpha,H["HalphaByEnvironment"]]],
         ["Hbeta",[Hbeta,H["MI_treeAndEnvironment"]]],              
        ["DistTurnover",[Hbeta/HE,H["DistTurnover"].iloc[1,0]]],
        ["DistTurnoverbySample",[Hbeta/HE,H["DistTurnoverBySample"].iloc[1,0]]]],
        columns=["Test","RegularRoutine"],
        orient="index")
    result["Dif"]=result.Test-result.RegularRoutine
    print(result)
    return None


if __name__ == '__main__':
    com={"-k":0}
    count=1
    key=None
    for i in sys.argv:
        if i[0]=="-":
            key=i
            com[key]=None
        elif key!=None:
            com[key]=i
            key=None
    print com
    if com.has_key("-f"):
        print "Test from user supplied sample. DistTurnover is wrong if more than two environment are present"
        FullCounts, t, otu,db=GetFromFiles(com)
    else:
        print "Test from Automatically generated data: 4 sample, two environments. Data from 4 identical dirichlets"
        FullCounts, t, otu,db=CreatingDataset(nsample=2)
    Report(FullCounts, t, otu,db)