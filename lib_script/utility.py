import sys

def concatGenerators(*n):
    while 1:
        res=[]
        for a in n:
            res.append(a.next())
        yield res

def Table(x,cases=[]):
    table={}
    for case in cases:
        table[case]=0
    if len(cases)!=0:
        table['NA']=0
    for obs in x:
        if table.has_key(obs):
            table[obs]+=1
        else:
            if not cases:
                table[obs]=1
            else:
                table['NA']+=1
    return table
def output_dict(dict, prefix='',file=''):
    if file=='':
        out=sys.stdout
    else:
        out=open(file,'w')
    for key in dict:
        out.write(prefix+str(key)+'\t')
        out.write(str(dict[key])+'\n')
        
def table(x,cases=[],table={}):
    if table=={}:
        for case in cases:
            table[case]=0
        if cases:
            table['NA']=0
    for obs in x:
        try:
            table[obs]+=1
        except KeyError:
            if not cases:
                table[obs]=1
            else:
                table['NA']+=1
    return table

class contatore:
    def __init__(self, casi):
        self.data={}
        for caso in casi:
            self.data[caso]=[0]
        if casi:
            self.data['NA']=[0]
    def conta(self,osservazioni):
        for obs in osservazioni:
            try:
                self.data[obs][-1]+=1
            except KeyError:
                table['NA'][-1]+=1
    def NuovaSerie(self):
        for k in self.data:
            self.data[k].append(0)
    

AmbDict={}
AmbDict['A']=set('A')
AmbDict['C']=set('C')
AmbDict['G']=set('G')
AmbDict['T']=set(['T'])
AmbDict['U']=set(['T','U'])
AmbDict['M']=set(['A','C'])
AmbDict['R']=set(['A','G'])
AmbDict['W']=set(['T','A'])
AmbDict['S']=set(['C','G'])
AmbDict['Y']=set(['T','C'])
AmbDict['K']=set(['T','G'])
AmbDict['V']=set(['A','C','G'])
AmbDict['H']=set(['A','C','T'])
AmbDict['D']=set(['A','T','G'])
AmbDict['B']=set(['T','C','G'])
AmbDict['N']=set(['A','C','G','T'])
AmbDict['?']=set(['A','C','G','T'])
AmbDict['-']=set(['-'])

def transpose(m):
    if isinstance(m, list):
        if isinstance(m[0], list):
            return map(list, zip(*m))
        else:
            return zip(*m) # faster
    else:
        if isinstance(m[0], list):
            return tuple(map(list, zip(*m)))
        else:
            return tuple( zip(*m) )
