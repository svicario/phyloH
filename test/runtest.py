#GoodTest
#Setting Test
from Bio import Phylo
from pandas import DataFrame
from numpy import random,array
from collections import OrderedDict
import sys
sys.path.append('../')
#esecutorePhyloHPandas = imp.load_source("esecutorePhyloHPandas","../esecutorePhyloHPandas.py")
from esecutorePhyloHPandas import *
com={"-f":"Echinodermata.tree","-g":"GroupTest","-s":"sampleTest"}
t=Phylo.read(com["-f"], "newick")
otu=[x.name for x in t.get_terminals()]
countsA=random.random(len(otu))
countsA=countsA/sum(countsA)
countsB=random.random(len(otu))
countsB=countsB/sum(countsB)
freqE=random.random(1)[0]
freqE=array([1-freqE,freqE])
FullCounts=array([countsA,countsB]).transpose()*freqE
N=round(1/(FullCounts.min())+10,ndigits=-2)
FullCounts=(FullCounts*N).astype("int")
countsA=FullCounts[:,0]/float(sum(FullCounts[:,0]))
countsB=FullCounts[:,1]/float(sum(FullCounts[:,1]))
freqE=sum(FullCounts)/float(FullCounts.sum())

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

traversing(t.clade)
A=DataFrame(countsA, index=otu)
B=DataFrame(countsB, index=otu)
desc=[]
traversingForDescendant(t.clade,desc)
Counts=[x[0:3]+[A.loc[x[3]].values.sum(),B.loc[x[3]].values.sum()] for x in desc]
Name, anc,Is_leaf, freqA,freqB=zip(*Counts)
Counts=DataFrame({"Name":Name, "Is_leaf":Is_leaf, "anc":anc, "freqA":freqA,"freqB":freqB}, index=Name)

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
totAs=[]
totBs=[]
for y in D:
    name, DepthAnc, DepthNode=y
    s=[x[0] for x in D if (min(x[-1], DepthNode)-max(x[1], OldDepthNode))>0]
    WindowSizes.append(DepthNode-OldDepthNode)
    OldDepthNode=DepthNode
    freqA=Counts.loc[s].freqA.values
    freqB=Counts.loc[s].freqB.values
    totAs.append(sum(freqA))
    totBs.append(sum(freqB))
    HA=entropy(freqA)
    HB=entropy(freqB)
    Halphas.append(HA*freqE[0]+HB*freqE[1])
    Hgammas.append(entropy(freqA*freqE[0]+freqB*freqE[1]))

WindowSizes=array(WindowSizes)
Halpha=sum(array(Halphas)*WindowSizes)/sum(WindowSizes)
Hgamma=sum(array(Hgammas)*WindowSizes)/sum(WindowSizes)
Hbeta=Hgamma-Halpha
HE=entropy(freqE)
 

handle=open("sampleTest","w")
handle.write("\n".join(["\t".join(["A",str(x[0]),x[1]]) for x in zip(FullCounts[:,0],otu)]))
handle.write("\n")
handle.write("\n".join(["\t".join(["B",str(x[0]),x[1]]) for x in zip(FullCounts[:,1],otu)]))
handle.close()
handle=open("GroupTest","w")
handle.write("A\tAA\nB\tBB")
handle.close()

db=DBdata()
db.readTreePandas(com['-f'])
db.readSampleTable(com["-s"])
db.readGroupTable(com["-g"])
H=db.GetEntropiesPandas(q="1", Pairwise=0, EqualEffort=0)
result=DataFrame.from_items([["Hgamma",[Hgamma,H["Hgamma"]]],
    ["Halpha",[Halpha,H["HalphaByEnvironment"]]],
     ["Hbeta",[Hbeta,H["MI_treeAndEnvironment"]]],              
    ["DistTurnover",[Hbeta/HE,H["DistTurnover"].loc["BB","AA"]]]],
    columns=["Test","RegularRoutine"],
    orient="index")
result["Dif"]=result.Test-result.RegularRoutine
print(result)