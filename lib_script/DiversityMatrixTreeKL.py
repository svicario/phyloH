#from ete2 import Tree
from Bio import Phylo
from numpy import log, exp, arange
import numpy
from pandas import DataFrame, MultiIndex, read_csv, Series, cut






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


def traversingForDescendant(c, desc):
    """
    getting from tree, name of  node, if is terminal, branch length, and list of descedant 
    """
    desc.append([c.name,c.is_terminal(), c.branch_length,[x.name for x in c.get_terminals()]])
    for C in c.clades:
        C=traversingForDescendant(C,desc)
    return c    

def GetDepths(t):
    """
    getting distance from the root, formatted as dataframe
    """
    Depths=t.clade.depths()
    Depths=Depths.items()
    Depths=[[x[0].name,x[1]] for x in Depths]
    Depths=DataFrame.from_records(Depths)
    return Depths

def TreeExtractInfoMultiple(t):
    DepthsList=[]
    descList=[]
    LList=[]
    for tt in t:
        Depths, desc, L=TreeExtractInfo(t)
        DepthsList.append(Depths)
        descList.append(desc)
        LList.append(L)
    
    return DepthsList, descList, LList
def TreeExtractInfo(t):
    """
    Extracting Info from the tree to make statistics over matrix
    """
    traversing(t.clade,count=[0])
    desc=[]
    traversingForDescendant(t.clade, desc)
    Depths=GetDepths(t)
    Name,Is_Leaf,BL,descArr=zip(*desc)
    L=DataFrame.from_records(zip(Name,Is_Leaf,BL), columns=["Name","Is_Leaf","BL"])
    L.index=L.Name
    desc=zip(Name,descArr)
    return Depths, desc, L
def TreeMatrixMultiple(DList,descList,LList, Env=None,DshapeLarge=True):
    DtreeList=[]
    for D,desc,L in zip(DList,descList,LList):
        DtreeList.append(TreeMatrix(D,desc,L))
    return DtreeList
def TreeMatrix(D,desc,L, Env=None,DshapeLarge=True):
    """
    Applying tree information (desc,L) on a given count matrix (D) and columns grouping (Env) to obtain a matrix of count over the tree
    """
    #I assume that D has correct columns name
    if not DshapeLarge:
        Z=(D[["Sample","Taxon"]]).values
        Z=MultiIndex.from_tuples(map(tuple,tuple(Z)), names=["Sample","Taxon"])
        D.index=Z
        Dlarge=D.Count.unstack(level=0)
        Dlarge.fillna(value=0,inplace=True)
        #I assume that Environment has correct index and columns names
        ExperimentalDesignColumns=MultiIndex.from_tuples(
        map(tuple,tuple(Env.ix[Dlarge.columns].values))
        , names=["Sample","Group"])
    else:
        # if D is already Large Environment information is already included
        Dlarge=D
        ExperimentalDesignColumns=Dlarge.columns
    #if taxon only present in tree but not in table, access mode .ix correctly report NA for that line, that later will be converted to zero.
    NodeTableLarge=[[x[0],Dlarge.ix[x[-1]].sum()] for x in desc]
    Dtree=DataFrame.from_items(NodeTableLarge).transpose()
    NodeAndLeafNamesIndex=MultiIndex.from_tuples(
        map(tuple,tuple(L.loc[:,["Name","Is_Leaf"]].ix[Dtree.index].values))
        , names=["Name","Is_Leaf"])
    Dtree.index=NodeAndLeafNamesIndex
    Dtree.columns=ExperimentalDesignColumns
    Dtree.columns=Dtree.columns.reorder_levels(["Group", "Sample"])
    return Dtree

def ChaoShannon(P,L,T, contribution=True):
    """
    Phylogenetic Entropy calculator 
    """
    W=(numpy.array([L.BL]).transpose())/(numpy.array([T]))
    #this to allow T and P to be unidimensional vector and not matrix with a single column
    if W.shape[1]==1:
        W=W[:,0]
    H=-((W*P*log(P)).fillna(value=0))
    if not contribution:
        H=H.sum()
    return H

def ChaoKL(P, Pbase,L,T, contribution=True):
    """
    Phylogenetic KL calculator 
    """
    W=(numpy.array([L.BL]).transpose())/(numpy.array([T]))
    #this to allow T and P to be unidimensional vector and not matrix with a single column
    if W.shape[1]==1:
        W=W[:,0]
    #H=DataFrame(W*P*log(P.values/DataFrame(3*[Pbase]).values.transpose())).sum()
    H=((W*P*log(P.values/Pbase.values)).fillna(value=0))
    if not contribution:
        H=H.sum()
    return H

def MatrixMultiple(DtreeList,Perm, Pairwise, EqualEffort):
    HList=[]
    for Dtree in DtreeList:
        HList.append(MatrixKL(Depths=None, desc=None,L=None,D=None, Perm=Perm, Pairwise=Pairwise, EqualEffort=EqualEffort, Dtree=Dtree))
        del HList[-1]["ITEi"]
        del HList[-1]["ITSgivenEi"]
    return HList
def MatrixKL(Depths, desc,L,D, Env=None, Perm=False, Pairwise=True, EqualEffort=False, Dtree=None):
    """
    Calculate phylogenetic entropy partitioning without correction
    """
    #Expand observation to nodes
    if (not Dtree):
        Dtree=TreeMatrix(D=D, desc=desc,L=L,Env=Env)
    
    #here something better should exist!! to get levels with a given name in ther order of the level
    groupLevels=list(Dtree.columns.levels[[x for x,y in enumerate(Dtree.columns.names) if y=="Group"][0]])
    sampleLevels=list(Dtree.columns.levels[[x for x,y in enumerate(Dtree.columns.names) if y=="Sample"][0]])
    
    #count by experimental design
    Tot=Dtree.loc["L0",:].values.sum()
    Tots=Dtree.loc["L0",:]
    Tote=Dtree.loc["L0",:].sum(axis=1, level="Group")
    
    #Frequency by group
    #frequency of the environment
    Pe=Tote/Tot
    #frequency of the samples
    Ps=Tots/Tot
    #Frequency of the sample within the environment
    Pise=Tots/Tote[Dtree.columns.get_level_values("Group")].values
    #print "EqualEffort", EqualEffort
    if EqualEffort:
        ng=len(groupLevels)
        ns=len(sampleLevels)
        Ps=ns*[1/float(ns)]
        Ps=DataFrame(Ps, index=Tots.columns, columns=Tots.index).transpose()
        Pe=Series(Dtree.columns.get_level_values("Group")).value_counts()
        Pe=DataFrame(Pe).transpose()
        Pe=Pe[Tote.columns]
        #Pe.columns=Tote.columns
        Pise=1/(1.0*Pe[Dtree.columns.get_level_values("Group")].values)
        Pise=DataFrame(Pise, index=Tots.index, columns=Tots.columns)
        Pe=Pe/ns
    
    #Frequency by group and nodes
    #note that only frequency calcualted from count is Pis
    #the other are derived such to use artificial Pe and Pise value for equal counts 
    Pis=Dtree/Dtree.loc["L0",:].values[0]
    Pie=(Pis*Pise.values[0]).sum(axis=1, level="Group")
    #Pie=Dtree.sum(axis=1, level="Environment")/Dtree.loc["L0",:].sum(axis=1, level="Environment")
    Pi=(Pie*Pe.values[0]).sum(axis=1)
    #Pi=Dtree.sum(axis=1)/Dtree.loc["L0",:].sum(axis=1)
    
    #Chao Normalization factor
    Ts=(L.BL.values*Pis.transpose()).transpose().sum()
    Te=(L.BL.values*Pie.transpose()).transpose().sum()
    T=(L.BL.values*Pi).sum()
    #print T, Te, Ts
    
    #Estimate statistics
    HTi=ChaoShannon(Pi,L,T)
    Hgamma=HTi.sum()
    HTgivenEi=(Pe.values[0]*ChaoShannon(Pie,L,Te)).sum(axis=1)
    HgammaForEachEnvironment=ChaoShannon(Pie,L,Te).sum(axis=0)
    Halpha=HTgivenEi.sum()
    HTgivenSi=(Ps.values[0]*ChaoShannon(Pis,L,Ts)).sum(axis=1)
    HalphaBySamples=HTgivenSi.sum()
    ITEi=HTi-HTgivenEi
    KL_TgivenETe=ChaoKL(Pie, DataFrame(Pi),L,Te).sum()
    ITE=ITEi.sum()
    ITSgivenEi=HTgivenEi-HTgivenSi
    ITSgivenE=ITSgivenEi.sum()
    res={"tot":Tot,"tag":Tote,
            "Hgamma":Hgamma, "HalphaBySamples":HalphaBySamples,
            "MI_treeAndEnvironment":ITE, "MI_treeAndSampleGivenEnvironment":ITSgivenE,
            "MI_KL":KL_TgivenETe,"HalphaByEnvironment":Halpha, "HgammaEachEnvironment":HgammaForEachEnvironment,
            "ITEi":ITEi, "ITSgivenEi":ITSgivenEi}
    if (not Perm):
        res["HE"]=-(Pe*log(Pe)).fillna(value=0).values.sum()
        res["HS"]=-(Ps*log(Ps)).fillna(value=0).values.sum()
        res["Pie"]=Pie
        res["Pis"]=Pis
        res["counts"]=Tots
        if Pe.shape!=Ps.shape:
            DistMatrix=DataFrame( columns=groupLevels,index=groupLevels)
            while groupLevels:
                G=groupLevels.pop()
                for g in groupLevels:
                    Pisd=Pis[[G,g]]
                    Ped=Pe[[G,g]]/Pe[[G,g]].values.sum()
                    Pied=(Pisd*Pise[[G,g]].values[0]).sum(axis=1, level="Group").fillna(0)[[G,g]]
                    Pid=(Pied*Ped[[G,g]].values[0]).sum(axis=1).fillna(0)
                    Ted=(L.BL.values*Pied.transpose()).transpose().sum()
                    Td=(L.BL.values*Pid).sum()
                    HTd=ChaoShannon(Pid,L,Td).values.sum()
                    HTgivenEd=(Ped*(ChaoShannon(Pied,L,Ted).sum())).values.sum()
                    hed=-(Ped*log(Ped)).fillna(value=0).values.sum()
                    DistMatrix.loc[G,g]=(HTd-HTgivenEd)/hed
            res["DistTurnover"]=DistMatrix
        if Pairwise:
            DistMatrix=DataFrame( columns=sampleLevels,index=sampleLevels)
            #make a copy
            template=DataFrame((Pis+0).fillna(0).values, index=Pis.index, columns=sampleLevels)
            #template.columns=template.columns.droplevel(level="Group")
            Etemplate=DataFrame((Ps+0).values, columns=sampleLevels)
            #print Etemplate
            while sampleLevels:
                S=sampleLevels.pop()
                for s in sampleLevels:
                    Pisd=template[[S,s]]
                    Ped=Etemplate[[S,s]]
                    #print "ciao"
                    #print Ped.sum().sum()
                    #print Ped
                    Ped=Ped/Ped.sum().sum()
                    #print Pisd
                    #Ped=Ps+0
                    #Ped.columns=Ped.columns.droplevel(level="Group")
                    #Ped=Ped[[S,s]]
                    #Ped=Pisd.sum(axis=0)
                    #Ped=Ped/sum(Ped)
                    #Pied=(Pisd*Pise[[G,g]]).sum(axis=1, level="Group").fillna(0)
                    #print Pisd
                    #print Pisd.shape, Ped.values
                    Pid=(Pisd*Ped.values[0]).sum(axis=1).fillna(0)
                    Tisd=(L.BL.values*Pisd.transpose()).transpose().sum()
                    Td=(L.BL.values*Pid).sum()
                    HTd=ChaoShannon(Pid,L,Td).values.sum()
                    HTgivenEd=(Ped*(ChaoShannon(Pisd,L,Tisd).sum())).values.sum()
                    hed=-(Ped*log(Ped)).fillna(value=0).values.sum()
                    #print HTd,HTgivenEd,hed
                    DistMatrix.loc[S,s]=(HTd-HTgivenEd)/hed
            res["DistTurnoverBySample"]=DistMatrix 
    return res

if __name__ == '__main__':
    t=Phylo.read("/Users/saverio/Downloads/test_acqua/rep_phylo_h2o.tre", "newick")
    #Convention of phylo to put 1 for branch length of root is not good.
    t.clade.branch_length=0
    D=read_csv("/Users/saverio/Downloads/test_acqua/test_phyloH_sample.txt", sep="\t", header=None, index_col=None,names=["Sample","Count","Taxon"])
    Env=read_csv("/Users/saverio/Downloads/test_acqua/map_h2o_M.txt", sep="\t", index_col=None, header=None,names=["Sample","Environment"])
    Env.index=Env.Sample
    #def ExtractInfo(treeete2):
    #    desc=[[x.name, x.dist, x.get_distance(t),x.get_leaf_names()] for x in treeete2.traverse() ]
    #Extraction is meant to be done once per tree and results converserved in the main data structure
    Depths, desc,L=TreeExtractInfo(t)
    #Statistics calculation is meant to be done after each permutation of the count matrix
    res=MatrixKL(Depths, desc,L,D, Env, DshapeLarge=False)
    print res
