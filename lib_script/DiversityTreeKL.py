from ete2 import Tree
from numpy import log, exp
import numpy

"""
funzioni per stimare la diversita'p di una comunita tenendo conto dell'albero filogenetico che li lega
si assume che l'oggetto albero e' un oggetto ete2 
La funzione getAbbundance assume che le abbondanze sono nell'ultimo token separato da underscore del nome 
La funzione puo' essere sostituita con altri sistemi, l'importante e' che il risultato sia uno scalare o un array numpy
"""

def getAbbundance(name,table=None):
    #print table.keys()
    #print table[name]
    try:
        return table[name]
    except TypeError:
        return int(name.split('_')[-1])

def traversingAbbundance(n, table=None):
    Abbundance=0
    if not n.children:
        #print n.name
        a=getAbbundance(n.name, table)
        n.add_feature('a',a)
        return a
    else:
        for N in n.children:
            Abbundance+=traversingAbbundance(N,table)
        n.add_feature('a',Abbundance)
        return Abbundance

    

def cleanTree(t):
    t.add_feature('T','NA')
    for n in t.traverse():
        n.add_feature('ag','NA')
        n.add_feature('a','NA')
        n.add_feature('atot','NA')
        n.add_feature('permbetaA',[])
        n.add_feature('permbetaCA',[])
        n.add_feature('betaCA',0)
        n.add_feature('betaA',0)

def getSignNode(t, alpha=0.05, rewrite=True):
    Sign=[]
    for n in t.traverse():
        try:
            if rewrite:
                raise AttributeError
            if n.pvalue<alpha:
                Sign.append(n)
        except AttributeError:
            L=len(n.permbetaA)
            p=sum([x>=n.betaA for x in n.permbetaA])/float(L)
            n.add_feature('pvalue',p)
            if n.pvalue<alpha:
                Sign.append(n)
    
    res=sorted(Sign, key=lambda x: x.betaA, reverse=True)
    return res

def getSignNodeMult(t, alpha=0.05, rewrite=True):
    Sign=[]
    L=float(len(t.permbetaA))
    for n in t.traverse():
        try:
            if rewrite:
                raise AttributeError
            n.pvalue
            Sign.append(n)
        except AttributeError:
            p=sum([x>=n.betaA for x in n.permbetaA])/L
            n.add_feature('pvalue',p)
            p=sum([x>=n.betaCA for x in n.permbetaCA])/L
            n.add_feature('pvalueCA',p)
            Sign.append(n)
    
    res=sorted(Sign, key=lambda x: (x.pvalue,x.betaA), reverse=True)
    m=float(len(res))
    Sign=[y for x,y in enumerate(res) if y.pvalue<(alpha/float((x+1)))]
    if len(Sign)==0:
        return Sign
    if alpha*len(Sign)/m < (1/L):
        print str(L)+' replicates are too low to evaluate correctly significance, you should increase to at least '+ str(1+m/(alpha*len(Sign)))
    return Sign

def getSignalNode(t):
    Sign=[]
    for n in t.traverse():
        if n.betaA>n.betaCA:
            Sign.append(n)
    return Sign
    
def getFPFN(t,simNodes, Sign, alpha=0.05):
    Vero=set(sum([x.get_leaf_names() for x in simNodes],[]))
    #Sign=getSignNode(t, alpha)
    Forse=[set()]
    for n in Sign:
        S=set(n.get_leaf_names())
        #S=Forse[-1].union(S)
        Forse.append(S)
    
    Forse.pop(0)
    FP=[len(x.difference(Vero)) for x in Forse]
    FN=[len(Vero.difference(x)) for x in Forse]
    TP=[len(Vero.intersection(x)) for x in Forse]
    L=len(t)
    TN=[L-sum(x) for x in zip(FP,FN,TP)]
    return {"FP":FP, "FN":FN, "TP":TP,"TN":TN}
    
def getFPFN(t,simNodes, Sign):
    Vero=set(sum([x.get_leaf_names() for x in simNodes],[]))
    Forse=[]
    for x in Sign:
        Forse.append(x.get_leaf_names())
    
    Forse=set(sum(Forse,[]))
    Res={}
    Res["FP"]=FP=len(Forse.difference(Vero)) 
    Res["FN"]=FN=len(Vero.difference(Forse))
    Res["TP"]=TP=len(Vero.intersection(Forse))
    L=len(t)
    Res["TN"]=TN=L-sum(Res.values())
    
    Res["accuracy"]=(TP+TN)/float(sum([FP,FN,TP,TN]))
    try:
        Res["precision"]=TP/float(TP+FP)
    except ZeroDivisionError:
        Res["precision"]=None
    try:
        Res["sensitivity"]=TP/float(TP+FN)
    except ZeroDivisionError:
        Res["sensitivity"]=None
    try:
        Res["specificity"]=TN/float(TN+FP)
    except ZeroDivisionError:
        Res["specificity"]=None
    try:
        Res["MCC"]=((TP*TN)-(FP*FN))/((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))**0.5
    except ZeroDivisionError:
        Res["MCC"]=None    
    return Res
    
     
        
              


def TreeKL(t, table, groups, q=1,record=True, Perm=False):
    """
    table contiene come chiavi tutti i nomi di specie con value un array dei conteggi nelle varie subpopolazione, in piu ha due chiavi speciali Ni con array dei conteggi totali per campione e Ns per specie
    groups dizionario di come i campioni si compongono in gruppi
    """
    print "sono Vecchio", __name__,table
    #minifunction that implement basic KL function
    def KLi(L,T,p,ptot):
        temp=(L/T)*p*log(p/ptot)
        temp[p==0]=0
        return temp 
    #minifunction that implent the atomic "base sum" of hill numbers
    def HillAtomic(L,p,q,T):
        #print L,p,T
        q=float(q)
        if q==1:
            H=-(L*p*log(p))/T
            if H.__class__==numpy.array([]).__class__:
                H[p==0]=0
            elif p==0:
                H=0
        else:
            H=L*(p/T)**q
        return H
    #getting total count/abbundance for sample (ta), environment (tag) and total (tatot)
    #multiplication per 1.0 makes array from integer to float
    ta=traversingAbbundance(t.get_tree_root(), table)*1.0
    
    tag=[]
    for g in groups.values():
            tag.append(sum(ta[g]))
    tag=numpy.array(tag)
    tatot=sum(ta)
    #print ta, tag, tatot
    for n in t.traverse():
        #n.a=n.a/totalcount
        n.add_feature('atot',sum(n.a))
        temp=[]
        for g in groups.values():
            temp.append(sum(n.a[g]))
        n.add_feature('ag',numpy.array(temp))    
    
    #from count to relative abbundance for each nodes
    #notice that total count and root count are the same
    for n in t.traverse():
        n.a=n.a/ta
        n.ag=n.ag/tag
        n.atot=n.atot/tatot

    #Calculate the T factor of correction of Chao formula for each node for sample, environment and total
    T={'Ta':0,'Tag':0,'Ttot':0}
    for n in t.traverse():
        T['Ta']+=n.dist*n.a
        T['Ttot']+= n.dist*n.atot
        T['Tag']+= n.dist*n.ag
    t.add_feature('T',T)
    #print T, t.atot,t.a,t.ag
    #Calculate indices 
    Hgamma=0
    HalphaC=0
    HalphaA=0
    HalphaAA=0
    betaCA=0
    betaA=0
    betaAkl=0
    for n in t.traverse():
        if sum(n.a)==0: 
            if record and not Perm:
                n.add_feature('betaA',0)
                n.add_feature('betaCA',0)
                n.add_feature('betaAkl',0)
            if Perm:
                try:
                    n.permbetaA.append(0)
                    n.permbetaCA.append(0)
                except AttributeError:
                    n.add_feature('permbetaA',[0])
                    n.add_feature('permbetaCA',[0])
            continue
        n.betaAkl=(tag/tatot)*KLi(L=n.dist,T=t.T['Tag'],p=n.ag,ptot=n.atot)
        betaAkl+=n.betaAkl
        #print n.betaAkl, t.T['Tag']
        #print tag/tatot
        HgammaAtomic=HillAtomic(L=n.dist,p=n.atot,q=q, T=t.T['Ttot'])
        Hgamma+=HgammaAtomic
        
        HalphaCAtomic=HillAtomic(L=n.dist,p=n.a,q=q,T=t.T['Ta'])
        HalphaC+=sum((ta/tatot)*HalphaCAtomic)
        
        HalphaAAtomic=HillAtomic(L=n.dist,p=n.ag,q=q,T=t.T['Tag'])
        HalphaAA+=HalphaAAtomic
        HalphaA+=sum((tag/tatot)*HalphaAAtomic)
        
        nodebetaA=HgammaAtomic-sum((tag/tatot)*HalphaAAtomic)
        nodebetaCA=sum((tag/tatot)*HalphaAAtomic)-sum((ta/tatot)*HalphaCAtomic)
        
        if record and not Perm:
            n.add_feature('betaA',nodebetaA)
            n.add_feature('betaCA',nodebetaCA)
        if Perm:
            try:
                n.permbetaA.append(nodebetaA)
                n.permbetaCA.append(nodebetaCA)
            except AttributeError:
                n.add_feature('permbetaA',[nodebetaA])
                n.add_feature('permbetaCA',[nodebetaCA])
        
        betaCA+=nodebetaCA
        betaA+=nodebetaA
        #print 'res', -n.dist*t.T['Ttot']*apa*log(apa)/(1-(1-apa)**sum(t.a)), sum(ni[n.a!=0]*-(n.dist*t.T['T'][n.a!=0]*pa[n.a!=0]*log(pa[n.a!=0]))/(1-(1-pa[n.a!=0])**t.a[n.a!=0]))
    return {"Hgamma":Hgamma, "HalphaBySamples":HalphaC, "MI_treeAndEnvironment":betaA, "MI_treeAndSampleGivenEnvironment":betaCA, "MI_KL":betaAkl,"HalphaByEnvironment":HalphaA, "HgammaEachEnvironment":HalphaAA}
if __name__=='__main__':
    t=Tree(sys.argv[1])
    R=TreeDiversity(t,q=1)
    for stat in R:
        print((stat,R[stat]))
    
