"""
perform partitioning of diversity and defining significance of grouping based on permutation
"""
import sys
#import scipy
import numpy
import math
from numpy import array, random, dtype
from utility import Table
#from DiversityTreeKLFabris import cleanTree, TreeKL
from DiversityMatrixTreeKL import TreeExtractInfo, MatrixKL
from ete2 import Tree
from time import time
from copy import deepcopy
from Bio import AlignIO, Phylo
#from ete2 import faces, NodeStyle
import os
import numpy
from pandas import DataFrame,MultiIndex, crosstab


def Diversity2Perc(D, Cardinality,q=1):
    res=(((1/D)**(q-1.000001)-(1/Cardinality)**(q-1.000001))/(1-(1/Cardinality)**(q-1.000001)))
    return res

def flatTree(db):
    for i in db.tree.traverse():
        if i.is_leaf():
            i.dist=1
        else:
            i.dist=0
def HPD(lista,alpha):
    n=len(lista)
    m=max(1,int(math.ceil(alpha*n)))
    lista.sort()
    a=lista[:m]
    b=lista[(n - m ):]
    return min([x for x in zip(*[a,b])],key=lambda x: x[1]-x[0])
    
    
def Hill(lista,q):
    H=0
    q=float(q)
    if q==1:
        for n in lista:
            if n==0: continue
            H+=-n*numpy.log(n)
        
        D=numpy.exp(H)
    else:
        for n in lista:
            H+=n**q
        
        D=H**(1/(1-q))
    return H,D



def entropy(table,q=1):
    freq=table.values()
    tot=float(sum(freq))
    if tot>1:
        freq=[x/tot for x in freq]
    H,D=Hill(freq,q)
    return H

def entropySh(table):
    freq=table.values()
    tot=float(sum(freq))
    Freq=freq/tot
    E=0
    for f in Freq:
        E+=f*numpy.log(f)
    return E


    
def alignmentEntropy(alignment, table,q=1):
    A=deepcopy(alignement_records)
    A=[[list(x.seq),x.description] for x in A]
    N,S=zip(*A)
    S=zip(*S)
    freqSeq=[table[n] for n in N]
    Compress=Table(S)
    H=0
    for sitepattern in Compress:
        Freq={'A':0 ,'C':0,'T':0,'G':0}
        OldfreqSeq=[]
        sitepattern=sitepattern.upper()
        for b in sitepattern:
            f=freqSeq.pop(0)
            try:
                Freq[b]+=f
            except KeyError:
                pass
        
        H+=Compress[sitepattern]*Hill(Freq.values())
    
    H=H/float(sum(Compress.values()))
    return H

def chaoshenCorr(table):
    n=0
    f1=0
    for c in table.values():
        n+=c
        if c==1:
            f1+=1
    if f1==n:
        f1=n-1
    
    C=1-f1/float(n)
    t={}
    for c in table:
        t[c]=(1-(1-C*table[c])**n)
    return C,t

def Dirichlet(table,a=1):
    for c in table:
        table[c]+=a
    return table

class DBdata:
    def __init__(self):
        self.countTable=[]
        self.samplesNames=[]
        self.SeqName=[]
        self.groups={}
        self.ChaoShenCorr=False
        self.withTree=True
        self.withGroups=True
    def noGroup(self):
        sample=dict([[y,[x]] for x, y in enumerate(self.samplesNames)])
        self.groups=sample
        self.expandTable()
    def readGroupTable(self,filename):
        """
               A given Sample belong always to the same environment, then each sample is cited only once in the file
               the file has two columns, 

               sample1  environment1
               sample2  environment1
               sample3  environment2
               sample4  environment2
        """
        sample=dict([[y,x] for x, y in enumerate(self.samplesNames)])
        assert len(sample)>1
        s=open(filename,'r').readlines()
        S=[x.split() for x in s]
        #assert len(set(S[0]))==len(S[0])
        for s, e in S:
            try:
                self.groups[e].append(sample[s])
            except KeyError:
                self.groups[e]=[sample[s]]
        self.expandTable()
    def readSampleTable(self,filename):
        """
        Sample file is a tabbed or space delimited table with 3 columns
        In first case the columns are the following:
        sample count Treeleaf
        
        combination sample and treeleaf not represented are considered count zero.
 
        there should be only one instance of combination of  Treeleaf and samplenumber
        A contingency table would be built based with one row for each Treeleaf and one columns for each sample, while several columns will belong to a given environment
        """
        D=DataFrame.from_csv(filename, header=None, index_col=None,sep="\t")
        Dc=crosstab(D.iloc[:,0],D.iloc[:,2],values=D.iloc[:,1], aggfunc=sum)
        Dc=Dc.transpose()
        Dc.fillna(0,inplace=True)
        self.samplesNames=map(str,Dc.columns.tolist())
        self.SeqName=map(str,Dc.index.tolist())
        self.countTable=Dc.values
    def OLDreadSampleTable(self,filename):
        """
        Sample file is a tabbed or space delimited table with 3 columns
        In first case the columns are the following:
        sample count Treeleaf
        
        combination sample and treeleaf not represented are considered count zero.
 
        there should be only one instance of combination of  Treeleaf and samplenumber
        A contingency table would be built based with one row for each Treeleaf and one columns for each sample, while several columns will belong to a given environment
        """
        try:
            s=filename.readlines()
        except AttributeError:
            s=open(filename,'r').readlines()
        S=[x.strip().split("\t") for x in s]
        SS=zip(*S)
        assert len(SS)==3, SS
        names=list(set(SS[2]))
        leaves=dict([[y,x] for x, y in enumerate(names)])
        self.samplesNames=list(set(SS[0]))
        self.SeqName=names
        sample=dict([[y,x] for x, y in enumerate(self.samplesNames)])
        
        row=len(leaves)
        col=len(sample)
        Table=array(row*col*[0])
        Table.shape=row,col
        for i in S:
            Table[leaves[i[2]],sample[i[0]]]=int(i[1])
        self.countTable=Table
    def readTree(self, filename, format=1):
        self.tree=Tree(filename, format)
    def readTreesPandasMultiple(self, filename):
        self.trees=list(Phylo.parse(filename, "newick"))
    def readTreePandas(self, filename):
        self.tree=Phylo.read(filename, "newick")
        self.TreeSummary=TreeExtractInfo(self.tree)
        self.TreeStat={"ITEi":DataFrame(index=self.TreeSummary[-1].index, columns=[]),
                       "ITSgivenEi":DataFrame(index=self.TreeSummary[-1].index, columns=[])}
    def readAlignment(self, filename):
        self.alignment=AlignIO.read(filename, 'fasta')
    def fromTree(self, filename,names,table, ambienti=2,replicates=2):
        #names,DB,Ns, table, t=FakeCommunity(Tree(filename),10000, ambienti=2,replicates=2, bonus=bonus)
        self.countTable=table
        self.SeqName=names
        self.samplesNames=range(table.shape[1])
        self.nogroup={0:self.samplesNames}
        counter=0
        for g in range(ambienti):
            self.groups[g]=range(counter,counter+replicates)
            counter+=replicates
        self.expandTable()
        if filename.__class__==''.__class__:
            self.readTree(filename)
        elif filename.__class__==Tree().__class__:
            self.tree=filename
        else:
            raise
        #return names,DB,Ns, table, t
        
    def fromSimulator(self, filename, N=10000, ambienti=2,replicates=2, bonus=0.1,Xs=[], prob=[]):
        t=Tree(filename)
        if not Xs:
            for a in range((ambienti-1)):
                X=SelectNode(list(t.traverse()))
                Xs.append(X)
        names,DB,Ns, table, t=FakeCommunity(t,10000, ambienti=ambienti,replicates=replicates, bonus=bonus, Xs=Xs, prob=prob)
        self.countTable=DB
        self.SeqName=names
        self.samplesNames=range(DB.shape[1])
        self.nogroup={0:self.samplesNames}
        counter=0
        for g in range(ambienti):
            self.groups[g]=range(counter,counter+replicates)
            counter+=replicates
        self.expandTable()
        self.readTree(filename)
        return names,DB,Ns, table, t
        
    def readTable(self,filename):
        """
        table file is a tabbed file with on each row a unique sequence, first columns its ID, then counts for each sample
        first row is special because have one column left and had the name of the groups on the same order of the label
        ex.
            group1    group1    group2    group3    group2
        id1    count11    count12    count13    count14    count15
        id1    count21    count22    count23    count24    count25
        ...
        idn    countn1    countn2    countn3    countn4    countn5
        """
        handle=open(filename,'r')
        groups=handle.readline()
        self.nogroup={0:db.samplesNames}
        self.groups
        for pos, group in enumerate(groups):
            try:
                groups[group].append(pos)
            except KeyError:
                groups[group]=[pos]
        countTable=zip([x.split('\t') for x in handle.readlines()])
        self.SeqName=countTable[0]
        self.countTable=numpy.array(zip(countTable[1:]))
    def expandTable(self, Groups=None):
        if not Groups:
            Groups=self.groups.items()
        self.itemTable=[]
        self.indexItemTable={'id':0,'sample':1,'group':2}
        for i in range(self.countTable.shape[0]):
            for g,columns in Groups:
                for c in columns:
                    self.itemTable+=int(self.countTable[i,c])*[[self.SeqName[i],c,g]]
        
        self.itemTable=numpy.array(self.itemTable)
        self.indexID={}
        D=zip(self.SeqName, range(self.countTable.shape[0]))
        print D[0]
        self.indexID.update(D)
        #self.Nread=self.itemTable.shape[0]
    def compressTable(self):        
        r, c=self.countTable.shape
        T=numpy.array(r*[c*[0]])
        ID=self.indexItemTable['id']
        SAMPLE=self.indexItemTable['sample']
        for i in self.itemTable:
            T[self.indexID[i[ID]],i[SAMPLE]  ]+=1
        
        #table={Ni:numpy.sum(T,0), Ns:numpy.sum(T,1)}
        table={}
        for i,name in zip(T,self.SeqName):
            table[name]=i
        return table
    def GetEntropiesPandas(self, record=True, Perm=False, q=1, branchScore=True, Pairwise=False, EqualEffort=False):
        #passare dall'itemtable alla tabella panda larga.
        A=DataFrame(["\t".join(x) for x in self.itemTable])
        AA=A.iloc[:,0].value_counts()
        AAindex=[x.split("\t") for x in AA.index]
        z=MultiIndex.from_tuples(AAindex, names=["Taxon","Sample","Group"])
        AA.index=z
        D=AA.unstack(level=0).transpose()
        ch={}
        ch.update(zip(map(str,range(len(self.samplesNames))),self.samplesNames))
        D.rename(columns=ch,inplace=True)
        #T=self.compressTable()
        Depths, desc, L=self.TreeSummary
        #print "UP",EqualEffort
        H=MatrixKL(Depths, desc,L, D, Perm=Perm, Pairwise=Pairwise, EqualEffort=EqualEffort)
        if Perm:
            ITEi=H["ITEi"]
            ITEi.name=len(self.TreeStat["ITEi"].columns)
            ITSgivenEi=H["ITSgivenEi"]
            ITSgivenEi.name=len(self.TreeStat["ITSgivenEi"].columns)
            self.TreeStat["ITEi"]=self.TreeStat["ITEi"].join(DataFrame(ITEi))
            self.TreeStat["ITSgivenEi"]=self.TreeStat["ITSgivenEi"].join(DataFrame(ITSgivenEi))
        return H
    
    def GetEntropies(self, record=True, Perm=False, q=1, branchScore=True):
        T=self.compressTable()
        groups={1:range(self.countTable.shape[1])}
        if self.withGroups:
            groups=self.groups
        try:
            H=TreeKL(t=self.tree, table=T, groups=groups,q=q, record=record, Perm=Perm,branchScore=branchScore)
        except KeyError:
            Diff=list(set(self.tree.get_tree_root().get_leaf_names()).symmetric_difference(T.keys()))
            Col=self.countTable.shape[1]
            T.update([[i,array(Col*[0])] for i in Diff])
            H=TreeKL(t=self.tree, table=T, groups=groups,q=q, record=record, Perm=Perm, branchScore=branchScore)
        return H
        #return TreeKLCS(t=self.tree, table=T, groups=groups,record=record, Chao=self.ChaoShenCorr, Perm=Perm)
    def RarefactionOfRead(self, nrarefazioni, repXrar=1, Perm=0, record=True):
        """
        nrarefazioni e' il numero di prove tra zero e tutti i read da fare
        oppure e' la lista in cui le dimensioni di campionamento da usare sono elecante  
        """
        self.expandTable()
        Nread=self.itemTable.shape[0]
        res=[]
        if nrarefazioni.__class__==[].__class__:
            Ns=nrarefazioni
        else:
            Ns=range(Nread,0, -Nread/nrarefazioni)
        for n in Ns:
            for r in range(repXrar):
                self.expandTable()
                sys.stdout.write(str(r)+' '+str(n)+'\n')
                selected=random.permutation(range((self.itemTable.shape[0])))[:n]
                #selected=random.randint(self.itemTable.shape[0], size=n)
                self.itemTable=self.itemTable[selected,:]
                tempCS=self.GetEntropies(record=record).items()
                #cleanTree(self.tree)
                #temp=self.GetEntropies(groups=groups, Tree=True, record=False, Chao=False)
                #cleanTree(self.tree)
                #tempH=self.PartionEntropy(groups=groups)
                #tempP=self.PermutationTest(reps=Perm)
                res+=[[r,self.itemTable.shape[0]]+tempCS]
                
        
        return res
#    def getCounts(self,type):
#        g=self.indexItemTable[type]
#        return Table(self.itemTable[:,g])
    def getEntropy(self, freq, entropyfunc=entropy, chaoshen=False, dirichlet=False, a=0.5,q=1):
        CorrH=1
        if chaoshen & q==1:
            freq, CorrH=chaoshenCorr(freq)
        if dirichlet:
            freq=Dirichlet(freq,a)
        if entropyfunc==TreeDiversity:
            Ent=entropyfunc(t=self.tree, q=q, table=freq)['Hp']
        else:
            Ent=entropyfunc(freq)
        
        return Ent
#    def PartionEntropies(self, groups=False, S=None, entropyfuncs=[TreeDiversity,entropy]):
#        res=[]
#        for entropyfunc in entropyfuncs:
#            res+=self.PartionEntropy(groups, S, entropyfunc)
#        return res
    def PermutationTest(self, reps=100, recordPerm=True,q=1):
        res=[]
        t0=time()
        OLD=deepcopy(self.itemTable)
        t1=0
        t2=0
        for i in range(reps):
            if i/100==i/100.0:
                sys.stdout.write('iteration '+str(i)+'\n')
            S=self.itemTable[:,self.indexItemTable['id']]
            t2a=time()
            numpy.random.shuffle(S)
            t2+=time()-t2a
            self.itemTable[:,self.indexItemTable['id']]=S
            t1a=time()
            res+=[ self.GetEntropies(Perm=recordPerm,q=q)]
            t1+=time()-t1a
        
        self.itemTable=OLD
        sys.stdout.write('total time:'+str(time()-t0)+' time calculating:'+str(t1)+' time shuffling:'+str(t2)+'\n')
        return res
    def PermutationTestPandas(self, reps=100, recordPerm=True,q=1):
        res=[]
        t0=time()
        OLD=deepcopy(self.itemTable)
        t1=0
        t2=0
        for i in range(reps):
            if i/100==i/100.0:
                sys.stdout.write('iteration '+str(i)+'\n')
            S=self.itemTable[:,self.indexItemTable['id']]
            t2a=time()
            numpy.random.shuffle(S)
            t2+=time()-t2a
            self.itemTable[:,self.indexItemTable['id']]=S
            t1a=time()
            res+=[ self.GetEntropiesPandas(Perm=recordPerm,q=q)]
            t1+=time()-t1a
        
        self.itemTable=OLD
        sys.stdout.write('total time:'+str(time()-t0)+' time calculating:'+str(t1)+' time shuffling:'+str(t2)+'\n')
        return res
    def BootstrapPandas(self, reps=100):
        res=[]
        t0=time()
        OLD=deepcopy(self.itemTable)
        t1=0
        t2=0
        for i in range(reps ):
            if i/100==i/100.0:
                sys.stdout.write('iteration '+str(i)+'\n')
            t2a=time()
            INDEX=numpy.random.random_integers(self.itemTable.shape[0]-1,size=self.itemTable.shape[0])
            self.itemTable=OLD[INDEX,:]
            t2+=time()-t2a
            t1a=time()
            temp=self.GetEntropiesPandas(Perm=False)
            res.append([temp[x] for x in temp if x in ["Hgamma","HalphaByEnvironment","MI_treeAndEnvironment"] ])
            t1+=time()-t1a
        
        self.itemTable=OLD
        Keys=["Hgamma","HalphaByEnvironment","MI_treeAndEnvironment"] 
        sys.stdout.write('total time:'+str(time()-t0)+' time calculating:'+str(t1)+' time resampling:'+str(t2)+'\n')
        return DataFrame(numpy.array(res), columns=Keys)
        
##    def PartionEntropy(self,S=None, entropyfunc=TreeKLCS ):
##        groups=True
##        if not self.withGroups:
##            groups=False
##        if S==None:
##            S=self.itemTable[:,self.indexItemTable['id']]
##        Hs=self.getEntropy(Table(S, cases=self.SeqName),entropyfunc=entropyfunc)
##        Hs_c=0
##        for s in self.samplesNames:
##            COUNTS=Table(S[self.itemTable[:,self.indexItemTable['sample']]==str(s)], cases=self.SeqName)
##            Hs_ci=self.getEntropy(COUNTS,entropyfunc=entropyfunc)
##            #sys.stdout.write(str(Hs_ci)+'\n')
##            print s,sum(COUNTS.values()), Hs_ci, float(self.itemTable.shape[0])
##            Hs_c+=sum(COUNTS.values())*Hs_ci/float(self.itemTable.shape[0])
##        Halpha=Hs_c
##        Hgamma=Hs
##        Hs_a=0
##        if groups:
##            for g in self.groups:
##                COUNTS=Table(S[self.itemTable[:,self.indexItemTable['group']]==str(g)], cases=self.SeqName)
##                print 'cc',g
##                Hs_a+=sum(COUNTS.values())*self.getEntropy(COUNTS,entropyfunc=entropyfunc)/float(self.itemTable.shape[0])
##            
##            HbetaC_A=Hs_a-Hs_c
##            HbetaA=Hs-Hs_a
##            return [Hgamma,Halpha, HbetaA,HbetaC_A ]
##        else:
##            Hbeta=Hs-Hs_c
##            return [Hgamma, Halpha, Hbeta]
#          
#    def PermutationTest(self, groups=False, reps=1000):
#        res=[]
#        t0=time()
#        S=deepcopy(self.itemTable[:,self.indexItemTable['id']])
#        for i in range(reps):
#            if i/100==i/100.0:
#                sys.stdout.write(str(i)+'\n')
#            numpy.random.shuffle(S)
#            res+=[ self.PartionEntropies(S=S,groups=groups)]
#        sys.stdout.write(str(time()-t0))
#        return res
#    def RarefactionOfRead(self, nrarefazioni, groups, repXrar=1, Perm=0):
#        """
#        nrarefazioni e' il numero di prove tra zero e tutti i read da fare
#        oppure e' la lista in cui le dimensioni di campionamento da usare sono elecante  
#        """
#        Nread=self.itemTable.shape[0]
#        res=[]
#        if nrarefazioni.__class__==[].__class__:
#            Ns=nrarefazioni
#        else:
#            Ns=range(Nread,0, -Nread/nrarefazioni)
#        for n in Ns:
#            for r in range(repXrar):
#                sys.stdout.write(str(r)+' '+str(n)+'\n')
#                selected=random.permutation(range((self.itemTable.shape[0])))[:n]
#                #selected=random.randint(self.itemTable.shape[0], size=n)
#                self.itemTable=self.itemTable[selected,:]
#                temp=self.PartionEntropies(groups=groups)
#                #tempP=self.PermutationTest(reps=Perm)
#                res+=[[r,self.itemTable.shape[0]]+temp]
#        
#        return res
def SelectNode(nodes):
    #selezioni un nodo proporzionale alla lunghezza del ramo su cui insiste
    l=sum([x.dist for x in nodes])
    S=0
    r=numpy.random.random_sample()
    for i in nodes:
            S+=i.dist/l
            if S>r:
                break
    return i


def FakeCommunity(t, N, replicates=2, ambienti=2, bonus=10, Xs=[], prob=None):
    from copy import copy
    names=t.get_leaf_names()
    Community=[]
    TheoreticalCommunity=[]
    if len(prob)==0:
        prob=numpy.exp(numpy.random.standard_normal(len(names)))
        prob=prob+min(prob)
        prob=prob/sum(prob)
    namedict={}
    d=zip(names,range(len(names)))
    namedict.update(d)
    for a in range(ambienti):
        TheoreticalCommunity.append(copy(prob))
        if (a!=0):
            Leaves2change=Xs[-a].get_leaf_names()
            for l in Leaves2change:
                TheoreticalCommunity[-1][namedict[l]]*=bonus
            
            TheoreticalCommunity[-1]=numpy.array(TheoreticalCommunity[-1])/float(sum(TheoreticalCommunity[-1]))
        for r in range(replicates):
            if r>0:
                TheoreticalCommunity.append(copy(TheoreticalCommunity[-1]))
            Community.append(numpy.random.multinomial(N,TheoreticalCommunity[-1]))
    
    TheoreticalCommunity=zip(*TheoreticalCommunity)
    TheoreticalCommunity=map(numpy.array,TheoreticalCommunity )
    table={}
    table.update(zip(names,TheoreticalCommunity))
    return names, numpy.array(Community).T,Xs,table, t



def layout2NHXw(node):
    M=max([x.betaA for x in node.get_tree_root().traverse()])
    V=int(round(10*(0.1+node.betaA/M)))
    node.add_feature('W',V)

def layout2NHXcol(node):
    node.add_feature('C',"255.0.0")
    

def layout_branch(node):
    style = NodeStyle()
    M=max([x.betaA for x in node.get_tree_root().traverse()])
    G=grad(11)
    V=int(round(10*(0.1+node.betaA/M)))
    
    style["hz_line_width"]=V
    style['hz_line_color']=G[V-1]
    node.img_style=style

def layout_Node(node):
    nstyle = NodeStyle()
    nstyle["shape"] = "sphere"
    nstyle["size"] = 10
    node.set_style(nstyle)

def tohex(r,g,b):
    hexchars = "0123456789ABCDEF"
    return "#" + hexchars[r / 16] + hexchars[r % 16] + hexchars[g / 16] + hexchars[g % 16] + hexchars[b / 16] + hexchars[b % 16]

def grad(cat):
    r=map(int,numpy.cumsum(numpy.array((cat-1)*[255.0/(cat-1)])))
    r=[0]+r
    b=[255-x for x in r]
    g=cat*[64]
    return [tohex(r,g,b) for r,g,b in zip(r,g,b)]

def parseTree(t, replicates=1):
    A=Tree(t)
    regions =set([x.split('_')[-3] for x in A.get_leaf_names() if len(x.split('_'))==4 ])
    res=[]
    #regions=['NI','SI','CI']
    #regions=['MPE4','MPE5']
    names=[]
    for n in A.get_leaf_names():
        temp=[]
        names.append(n)
        try:
            prefix, sample, cluster, count=n.split('_')
            temp=map(int,[[0,str(count)][sample==n] for n in regions])
        except ValueError:
            temp=map(int,[0]*len(regions))
        
        res.append(temp)
    
    return names, numpy.array(res), len(regions), replicates

def parseTreeWithDB(t,Ordername, method='PitTrap'):
    A=Tree(t)
    A.get_leaf_names()
    con=MySQLdb.connect(host='cerbero.ba.itb.cnr.it', user='prin2007',passwd='PyroNoise', db='TAXONOMYdb')
    cur=con.cursor()
    sql="""SELECT 454Reads.ReadAccno, 454Reads.Region, 454Reads.Rich
    FROM 454Reads
    INNER JOIN BestHit_Order ON BestHit_Order.QueryName = 454Reads.ReadAccno
    AND 454Reads.Run = '"""+method+"""'
    AND BestHit_Order.order_name = '"""+Ordername+"""'"""
    cur.execute(sql)
    results=cur.fetchall()
    DB={}
    DB.update([[x,[y,z]] for x in results])
    res=[]
    names=[]
    for n in A.get_leaf_names():
        temp=[]
        names.append(n)
        try:
            region,abb=DB[n]
            if region==regions[0]:
               temp=[abb,0]
            else:
                temp=[0,abb]
        except KeyError:
            temp=[0,0]
        res.append(temp)
    
    return names, numpy.array(res).T






    
                
        
