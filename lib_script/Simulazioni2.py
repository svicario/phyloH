from ete2 import Tree, nexml
from lib_script.entropyPerm5 import *
from lib_script.DiversityTreeKLFabris import *
#from lib_script import Itol,ItolExport


filename='/Users/saverio/work/Saccone/PRIN2007/FinalFusariumDB_2.nex.con'
RES={}
for b in range(1,10):
    db=DBdata()
    names,DB,Ns, table, t =db.fromSimulator(filename,10000, ambienti=2,replicates=2, bonus=float(b)/3.0,Xs=[], prob=prob)
    len(sum([x.get_leaf_names() for x in Ns],[]))
    Effect=0
    for n in Ns:
        K=n.get_leaf_names()
        for k in K:
            Effect+= abs(sum(table[k][db.groups[0]])-sum(table[k][db.groups[1]]))
    
    Effect
    Theo=db.GetEntropies(q=1)
    
    reps=1000
    reps=10
    RarG=[]
    SN={}
    N={}
    #Position={}
    R=[2500,5000,10000]
    R=[5000]
    #Pvalues={}
    for r in R:
        for f in ["Chao", "Shannon"]:
            db.tree=t
            if f=="Shannon":
                flatTree(db)
            cleanTree(db.tree)
            RarG+=db.RarefactionOfRead([r])
            Per=db.PermutationTest(reps=reps)
            SN[r]=getSignNodeMult(db.tree)
            N[r]={}
            N[r][f]=getFPFN(db.tree,Ns,SN[r])
            N[r][f]["Pvalue"]=sum([x["MI_treeAndEnvironment"]>RarG[-1][6][1] for x in Per])/float(len(Per))
            pos=[x for x, y in enumerate(SN[r]) if y.get_leaf_names()==Ns[0].get_leaf_names()]
            try:
                N[r][f]["Position"]=pos[0]
            except IndexError:
                N[r][f]["Position"]="NA"
    
    RES[Effect]=N


##header=["Effect","SampleSize"]
##TAB=[]
##for i in RES:
##    for r in RES[i][0]:
##        temp=[i,r]
##        for j in RES[i][0][r]:
##            temp.append(RES[i][0][r][j])
##        
####        temp.append(RES[i][1][r])
####        try:
####            temp.append(RES[i][2][r][0])
####        except IndexError:
####            temp.append("NA")
##        TAB.append(temp)
##
##header+=RES[i][0][r].keys()#+["Pvalue","Position"]
##
    
header=["Effect","SampleSize", "StatType"]
TAB=[]
for effect in RES:
    for Size in RES[effect]:
        for statType in RES[effect][Size]:
            Line=[effect,r,statType]
            for i in RES[effect][Size][statType]:
                Line.append(RES[effect][Size][statType][i])
            
            TAB.append(Line)

header+=RES[i][r].keys()


handle=open("Simulation.txt","w")
handle.write(TAB)
handle.close()

