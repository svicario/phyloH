from ete2 import Tree, nexml
from lib_script.entropyPerm5 import *
from lib_script.DiversityTreeKLFabris import *
from lib_script import Itol,ItolExport

import sys



def Diversity2Perc(D, Cardinality,q=1):
    res=(((1/D)**(q-1.000001)-(1/Cardinality)**(q-1.000001))/(1-(1/Cardinality)**(q-1.000001)))
    return res

def hsv2rgb (h, s,v):
    from math import floor
    if (h == 360):h=0     
    if (h == -1):h=0 
    h =h/60.0 
    i =floor(h)
    f =h - i
    p1 =v*(1-s)
    p2 =v*(1-(s*f))
    p3 =v*(1-(s*(1-f)))
    CMD={ 0:lambda v,p1,p2,p3: [int(255*i) for i in [v,p3,p1]],
          1:lambda v,p1,p2,p3: [int(255*i) for i in[p2,v,p1]],
          2:lambda v,p1,p2,p3: [int(255*i) for i in[p1,v,p3]],
          3:lambda v,p1,p2,p3: [int(255*i) for i in[p1,p2,v]],
          4:lambda v,p1,p2,p3: [int(255*i) for i in[p3,p1,v]],
          5:lambda v,p1,p2,p3: [int(255*i) for i in[v,p1,p2]]}
    return CMD[i](v,p1,p2,p3)

def spacedColors(cat):
    step=360.0/cat
    h=[step*x for x in range(cat)]
    colors=[tohex(*hsv2rgb(i,s=1,v=1)) for i in h]
    return colors


def tohex(r,g,b):
    hexchars = "0123456789ABCDEF"
    return "#" + hexchars[r / 16] + hexchars[r % 16] + hexchars[g / 16] + hexchars[g % 16] + hexchars[b / 16] + hexchars[b % 16]

def fromhex(string):
    hexchars = "0123456789ABCDEF"
    string=string[1:]
    A=[hexchars.find(x)+1 for x in string]
    r=A[0]*A[1]
    g=A[2]*A[3]
    b=A[4]*A[5]
    return {"red":r-1,"green":g-1,"blue":b-1}
##def grad(cat, switch=True):
##    r=map(int,numpy.cumsum(numpy.array((cat-1)*[255.0/(cat-1)])))
##    r=[0]+r
##    b=[255-x for x in r]
##    g=cat*[64]
##    if not switch:
##        temp=r
##        r=b
##        b=g
##        g=temp
##    return [tohex(r,g,b) for r,g,b in zip(r,g,b)]

def grad(cat, switch=True):
    g=[0]+map(int,numpy.cumsum(numpy.array((cat-1)*[255.0/(cat-1)])))
    g.reverse()
    r=cat*[255]
    b=cat*[0]
    return [tohex(r,g,b) for r,g,b in zip(r,g,b)]

def layout_branch(node,f):
    try:
        style = NodeStyle()
    except NameError:
        style={}
    
    M=max([x.__getattribute__(f) for x in node.get_tree_root().traverse()])
    G=grad(11)
    V=int(round(10*(0.1+node.__getattribute__(f)/M)))
    if V<=0:V=1
    style["hz_line_width"]=V
    style['hz_line_color']=G[V-1]
    style['size']=0
    style["bgcolor"]='#FFFFFF'
    node.img_style=style

def layout_branchB(node):
    DIF=node.betaA-node.betaCA    
    M=max([abs(x.betaA-x.betaCA) for x in node.get_tree_root().traverse()])
    G=grad(11)
    V=int(round(10*(0.1+(DIF)/M)))
    try:
        style = NodeStyle()
    except NameError:
        style={}
    style["hz_line_width"]=abs(V)
    style['hz_line_color']=G[V-1]
    style['size']=0
    style["bgcolor"]='#FFFFFF'
    node.img_style=style

def BuildNEXML(string_tree, otusName="myTaxa"):
    nexml_project = nexml.Nexml()
    tree_collection = nexml.Trees()
    nexml_tree = nexml.NexmlTree(string_tree)
    
    #N=nexml_tree.get_leaf_names()
    #Ntaxon=[nexml.Taxon(label=x,id=x) for x in N]
    #nexml_Taxa=nexml.Taxa(otu=Ntaxon)
    
    # Jaime: OTUs must be defined within the NexML project, so they are
    # not automatically create when a newick tree is processed. However,
    # you can create OTUs and bind them to the tree very easly:
    nexml_Taxa=nexml.Taxa(id=otusName)
    counter=0
    for leaf in nexml_tree.iter_leaves():
        # create a new Taxon (otu) nexml object and set its ID to leaf
        # name
        idotu="otu"+str(counter)
        counter+=1
        Ntaxon= nexml.Taxon(label=leaf.name,id=idotu)
        # Bind nexml_node of this leaf to the same OTU id
        leaf.nexml_node.otu = idotu
        # And add OTU to our Taxa group object
        nexml_Taxa.add_otu(Ntaxon)
    
    nexml_project.add_otus(nexml_Taxa)
    tree_collection.add_tree(nexml_tree)
    tree_collection.set_otus(otusName)
    nexml_project.add_trees(tree_collection)
    return nexml_project, nexml_tree

def Features2NeXML(node,features,datatype="double"):
    for f in features:
        value=node.__getattribute__(f)
        nexml_meta = nexml.LiteralMeta(datatype=datatype, property=f, content=str(value))
        node.nexml_edge.add_meta(nexml_meta)
        return value



def GraphFeatures2Itollabel(tree_obj,SN):
    #out=["\t".join(["label","hz_line_color"])]
    out=[]
    m=[x.__getattribute__("betaA") for x in tree_obj.traverse()]
    M=max(m)
    m=min(m)
    A=[[x,[str(round(M*((x/10.0)-0.15),3)),str(round(M*((x/10.0)-0.05),3))]] for x in range(1,12,1)]
    D={}
    D.update(A)
    D[1][0]=str(round(m,3))
    for n in tree_obj.traverse():
        out.append("\t".join(map(str,[n.name,"range",n.img_style["hz_line_color"]," - ".join(D[n.img_style["hz_line_width"]])])))
        if n in SN:
            out.append("\t".join([n.name,"clade","#00FFFF","Signif"]))
        else:
            out.append("\t".join([n.name,"clade","#000000","NotSignif"]))
    
    S="\n".join(out)
    return S

def HistFeatures2Itol(db):
    tree_obj=db.tree
    table=db.compressTable()
    abbundances=[]
    for leaf in table.keys():
        a=[leaf]
        for g in db.groups:
            a.append(sum(table[leaf][db.groups[g]]))
        abbundances.append(a)
    names=db.groups.keys()
    out=[",".join(["LABELS"]+names), ",".join(["COLORS"]+spacedColors(len(names)))]
    for a in abbundances:
        out.append(",".join(map(str,a)))
    return "\n".join(out)

def SignFeature2Itol(SN):
    out=[]
    for i in SN:
        out.append("\t".join([i.name,"Signif","#00FFFF"]))
    return "\n".join(out)

def mapWidthANDColor2Phyloxml(node):
    color=phyloxml.BranchColor(**fromhex(node.img_style["hz_line_color"]))
    node.phyloxml_clade.set_color(color)
    node.phyloxml_clade.set_width(node.img_style["hz_line_width"])




com={"-q":1, "-qt":None, "-x":"nexml","-R":"/opt/exp_soft/uniba/R/bin/","None":1}
count=1
key=None
for i in sys.argv:
    if i[0]=="-":
        key=i
        com[key]=None
    elif key!=None:
    	com[key]=i
    	key=None

spiegazione=""" 
 -f filename	   Use this file as the phylogeny file [phylo].
 -o filename       Use this file to record output
 -s filename       Use this file as the sample file [sample].
 -g filename       Use this file as the group file [group]
 -r INT		   Number of randomizations to use [999]
 -q float          q parameter in the hill series (q=1 index is beta is Chao phylogenetic entropy, q=2 Rao phylogenetic diversity,q  zero is faith phylogenetic diversity)
 -x string         two possible strings :"nexml" or "phyloxml" to select the xml output of the results
 -h 0 or 1         boolean to check if you want html output     
"""
if (not com.has_key('-s')) | (not com.has_key('-f')):
    print spiegazione
    raise

#Check if pyQT4 is there
try:
    from ete2 import NodeStyle, TreeStyle, faces
except ImportError:
    del com["-qt"]
#Loading the data
db=DBdata()
db.readSampleTable(com["-s"])
db.readTree(com['-f'])
if com.has_key("-g"):
    db.readGroupTable(com["-g"])

if com["-x"]=="nexml":
    nexml_project, db.tree=BuildNEXML(db.tree.write())
elif com["-x"]=="phyloxml":
    from ete2 import Phyloxml, phyloxml
    phylo = phyloxml.PhyloxmlTree(newick=db.tree.write())
    
    project = Phyloxml()
    project.add_phylogeny(phylo)
    db.tree=phylo

#add node label
for i,j in enumerate(db.tree.traverse()):
    if not j.is_leaf():
        j.name="L"+str(i)
    



#Getting Observed Values
R=db.GetEntropies(q=com["-q"])
print R
#Getting Confindence Interval
#bot=db.Bootstrap(reps=com["-r"])
#cleanTree(db.tree)

#Getting entropy evaluation of non Phylogenetic part of the data. S stand for sample and E for environment where sampled

countbyS=sum(db.countTable)
HS=Hill(countbyS/float(sum(countbyS)),com["-q"])[0]
countbyGroup=[]
for g in db.groups:
    countbyGroup.append(sum(countbyS[db.groups[g]]))
HE=Hill([x/float(sum(countbyGroup)) for x in countbyGroup],com["-q"])[0]

#Getting Significance
Per=db.PermutationTest(reps=int(com["-r"]),q=com["-q"])
SN=getSignNodeMult(db.tree, alpha=0.05)

db.tree.get_children()[0]
#Formatting Output
##B=[]
##for i in bot:
##    b=i.pop("HgammaEachEnvironment")
##    a=i.items()
##    a.sort()
##    a=zip(*a)[1]
##    B.append(list(a)+list(b))
##
##B=zip(*B)

P=[]
for i in Per:
    b=i.pop("HgammaEachEnvironment")
    a=i.items()
    a.sort()
    a=zip(*a)[1]
    P.append(list(a)+list(b))

P=zip(*P)
Names=i.keys()
Names.sort()
Names+=["HgammaFor_"+g for g in db.groups.keys()]


Venn="""
<?xml version="1.0" encoding="windows-1252" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 20001102//EN" "http://www.w3.org/TR/2000/CR-SVG-20001102/DTD/svg-20001102.dtd" [
]>
<svg width="600" height="500" x="0" y="0" xmlns="http://www.w3.org/2000/svg" onmouseup="add(evt)" onmousedown="grab(evt)" onmousemove="null">
	<title>Report of PhyloH with use  of BioVenn and Itol graphics</title>
	<desc></desc>
	<script>
	<![CDATA[
	var Root=document.documentElement
	standardize(Root)
	function standardize(R){
		var Attr={
			'onmouseup':'add(evt)',
			'onmousedown':'grab(evt)',
			'onmousemove':null
		}
		assignAttr(R,Attr)
	}
	function grab(evt){
		var O=evt.target
		var Attr={
			'onmousemove':'slide(evt,"'+O.id+'")',
			'onmouseup':'standardize(Root)'
		}
		assignAttr(Root,Attr)
	}
	function slide(evt,id){
		if(id!='rect'&&id!='circlex'&&id!='circley'&&id!='circlez'){
			var o=document.getElementById(id)
			o.setAttributeNS(null, 'x', evt.clientX)
			o.setAttributeNS(null, 'y', evt.clientY)
		}
	}
	function assignAttr(O,A){
		for (i in A) O.setAttributeNS(null,i, A[i])
	}
	]]>
	</script>
	<rect id="rect" x="0" y="0" width="500" height="500" style="fill:#FFFFFF"/>
	<circle id="circlex" r="130.35407334401" cx="315.30560225805" cy="205.48651044486" style="fill:#FF0000;opacity:.5"/>
	<circle id="circley" r="140.7984913873" cx="255.17612189703" cy="349.76907151185" style="fill:#00FF00;opacity:.5"/>
	<circle id="circlez" r="92.174249216841" cx="303.89253831671" cy="349.76907151185" style="fill:#0000FF;opacity:.5"/>
	<text id="titlex" x="170" y="125" font-size="12" font-family="Courier New bold" font-weight="bold" style="fill:#000000">Haplotypes on Tree(T)</text>
	<text id="titley" x="27" y="306" font-size="12" font-family="Courier New bold" font-weight="bold" style="fill:#000000">Sample names (S)</text>
	<text id="titlez" x="401" y="358" font-size="12" font-family="Courier New bold" font-weight="bold" style="fill:#000000">Environment/Group names (E)</text>
	<text id="titleMITSgivenE" x="203" y="240" font-size="12" font-family="Courier New bold" font-weight="bold" style="fill:#000000">I(T,S|E)="""+str(round(R["MI_treeAndSampleGivenEnvironment"],3))+"""</text>
	<text id="titleMITE" x="267" y="301" font-size="12" font-family="Courier New bold" font-weight="bold" style="fill:#000000">I(T,E)="""+str(round(R["MI_treeAndEnvironment"],3))+"""</text>
	<text id="titleHTE" x="323" y="246" font-size="12" font-family="Courier New bold" font-weight="bold" style="fill:#000000">H(T|E)="""+str(round(R["HalphaByEnvironment"],3))+"""</text>
	<text id="titleHTS" x="302" y="176" font-size="12" font-family="Courier New bold" font-weight="bold" style="fill:#000000">H(T|S)="""+str(round(R["HalphaBySamples"],3))+"""</text>
	<text id="titleHS" x="134" y="363" font-size="12" font-family="Courier New bold" font-weight="bold" style="fill:#000000">H(S|E,T)="""+str(round(HS-HE-R["MI_treeAndSampleGivenEnvironment"],3))+"""</text>
	<text id="titleHE" x="292" y="367" font-size="12" font-family="Courier New bold" font-weight="bold" style="fill:#000000">H(E|T)="""+str(round(HE-R["MI_treeAndEnvironment"],3))+"""</text>
</svg>
"""

OUT="<res>"
OUT+="<HSample><Observed>"+str(HS)+"</Observed><Diversity>"+str(exp(HS))+"</Diversity></HSample>"
OUT+="<HEnvironment><Observed>"+str(HE)+"</Observed><Diversity>"+str(exp(HE))+"</Diversity></HEnvironment>"
outdict={"HEnvironment":{"Observed":HE,"Diversity":exp(HE),"MaxDiversity":len(db.groups.keys())},"HSample":{"Observed":HS,"Diversity":exp(HS),"MaxDiversity":len(countbyS)}}
count=0
for label, Pvalues in zip(Names,P):
    OUT+="<"+label+">\n"
    outdict[label]={}
    if label.find("HgammaFor")>-1:
        Nlabel="HgammaEachEnvironment"
        OUT+="<Observed>"+str(R[Nlabel][count])+"</Observed>\n"
        outdict[label]["Observed"]=R[Nlabel][count]
        D=exp(R[Nlabel][count])
        if label.find("alpha")>-1: name="alpha"
        elif label.find("gamma")>-1: name="gamma"
        OUT+="<"+name+"_Diversity>"+str(D)+"</"+name+"_Diversity>\n"
        outdict[label][name+"_Diversity"]=D
        count+=1
    else:
        OUT+="<Observed>"+str(R[label])+"</Observed>\n"
        outdict[label]["Observed"]=R[label]
        if label.find("MI")>-1:
            outdict[label]["Pvalue"]=sum(R[label]<=Pvalues)/float(len(Pvalues))
            OUT+="<Pvalue>"+str(sum(R[label]<=Pvalues)/float(len(Pvalues)))+"</Pvalue>\n"
            D=exp(R[label])
            OUT+="<beta_Diversity>"+str(D)+"</beta_Diversity>\n"
            outdict[label]["beta_Diversity"]=D
            Cardinality=exp(HE)
            S=Diversity2Perc(D,Cardinality)
            OUT+="<Percentage_Overlap>"+str(S)+"</Percentage_Overlap>\n"
            outdict[label]["Percentage_Overlap"]=S
        else:
            D=exp(R[label])
            OUT+="<alpha_Diversity>"+str(D)+"</alpha_Diversity>\n"
            outdict[label]["alpha_Diversity"]=D

    
    OUT+="</"+label+">\n"

R["MI_KL"]=list(R["MI_KL"])
for temp, g in zip(R["MI_KL"],db.groups.keys()):
    outdict["KL_of_"+g]={"Observed":temp}
    OUT+="<HbetaFrom_"+g+">"+str(temp)+"</HbetaFrom_"+g+">"

OUT+="<Counts><Total>"+str(sum(countbyS))+"</Total></Counts>"


OUT+="</res>"
handle=open(com["-o"]+".xml","w")
handle.write(OUT)
handle.close()
if bool(int(com["-h"])):
    CountHTML='<table border=1>'
    CountHTML+='<tr><th colspan='+str(len(countbyS))+'>total</th></tr>'
    CountHTML+='<tr><td align="center" colspan='+str(len(countbyS))+'>'+str(sum(countbyS))+'</td></tr>'
    CountHTML+="</tr>"
    for k in db.groups.keys():
            CountHTML+='<th colspan='+str(len(countbyS[db.groups[k]]))+'>'+k+'</th>'
    CountHTML+="</tr>"
    for k in db.groups.keys():
            CountHTML+='<td align="center" colspan='+str(len(countbyS[db.groups[k]]))+'>'+str(sum(countbyS[db.groups[k]]))+'</td>'
    CountHTML+="</tr>"
    for k in db.groups:
        for s in db.groups[k]:
            CountHTML+='<td align="center">'+str(countbyS[s])+"</td>"
    CountHTML+="</tr></table>"                                                                         

    translator={"HEnvironment":"H(E)","HSample":"H(S)","Hgamma":"H(T)","HalphaByEnvironment":"H(T|E)","HalphaBySamples":"H(T|S)","MI_treeAndEnvironment":"I(T,E)","MI_treeAndSampleGivenEnvironment":"I(T,S|E)"}
    for g in db.groups.keys():
        translator["HgammaFor_"+g]="H(T|E="+g+")"

    def makeHTMLTable(outdict,pattern, translator=None):
        from numpy import floor, log10
        def signif(x,ndigits=1):
            coef=1
            if x==0:
                return 0.0
            if x<0:
                coef=-1
            return round(x, -1+ndigits-int(floor(log10(coef*x))))
        outhtml=""
        KEYs=[k for k in outdict.keys() if k.find(pattern)>-1]
        KEYs.sort()
        outhtml+="""<table border="1"><tr>"""
        for k in KEYs:
            outhtml+="<th colspan="+str(len(outdict[k]))+">"+str(k)+"</th>"

        outhtml+="</tr>\n<tr>"
        for k in KEYs:
            for kk in outdict[k].keys():
                if (kk=="Observed")&bool(translator):
                    outhtml+="<th>"+translator[k]+"</th>"
                else:
                    outhtml+="<th>"+str(kk)+"</th>"

        outhtml+="</tr>\n"
        for k in KEYs:
            for kk in outdict[k].keys():
                try:
                    outhtml+='<td align="center">'+str(signif(outdict[k][kk],3))+"</td>"
                except TypeError:
                    outdict[k][kk]

        outhtml+="</tr></table>"
        return outhtml




    outhtml="<!DOCTYPE html><html><body><H1>Partitioning information in  "+com["-o"] +"</H1><H2>Global Statistics</H2>"

    outhtml+="<H3>Experimental Design:Counts of observations across groups and samples within groups</H3>"
    outhtml+= CountHTML                                                                           
    outhtml+="""<H3>Experimental Design Diversity: entropy and diversity of observation  in the different groups</H3>"""
    outhtml+=makeHTMLTable(outdict,"HEn",translator)
             
    outhtml+="""<H3>Experimental Design Diversity: entropy and diversity of observation  in the sample within the groups</H3>"""
    outhtml+=makeHTMLTable(outdict,"HSa",translator)
             
    outhtml+="""<H3>Gamma Diversity: diversity using all data and in each group</H3>"""
    outhtml+=makeHTMLTable(outdict,"gamma",translator)
             
    outhtml+="""<H3>Alpha Diversity: mean within group diversity</H3>"""
    outhtml+=makeHTMLTable(outdict,"alpha",translator)

    outhtml+="""<H3>Beta Diversity or Mutual Information between the phylogeny and a given grouping: diversity across group  and across sample within same group</H3>"""
    outhtml+=makeHTMLTable(outdict,"MI_t",translator)
    outhtml+="<H3>Difference of each group from total: phylogenetic Kullback-Leiber distance between each group and the overall sample</H3>"
    outhtml+=makeHTMLTable(outdict,"KL_of_")
    outhtml+="<H2>Venn diagram of the partitioning of the information across the three attributes present on each observation</H2>"
    outhtml+="N.B. 1)the diagram is not area-proportional. </p><p>2) E is always within S, given that each sample belong to only one environment type or sample group.</p>"
    outhtml+="<p>3)Hgamma=HalphabyEnvironment+Hbeta = H(T)=H(T|E)+I(T,E)</p>"
    outhtml+="<p>while taking in account sample info Hgamma=HalphaBySample+HbetabySamplegivenEnvironment+Hbeta = H(T)=H(T|S)+I(T,S|E)+I(T,E)</p>"
    outhtml+=Venn
    outhtml+="<H2>Per Node Statistics mapped on the phylogeny</H2>"
    outhtml+="""
                <p>Three types of data are shown on the tree.</p>
                <p>1)The color of the branches cyan indicates a contribution to I(T,E) higher than the null distribution,
                while branches are black otherwise.</p>
                <p>2)The background of each branch is a gradient from yellow to red for increased contribution to I(T,E). For details look at the legend on the side</p>
                <p>3)Bar plot on each tips indicates the number of count in each group</p>
                <p>Look at the tree find an relevant branches and text search the label of the branch to access the correct row on the by node statistics table.</p>
                <p>Go itol using the link to modify the tree, or use the itol table and the labelled tree to add further data set (i.e. taxonomic name)</p>
                """
    outhtml+='<img src="'+com["-o"]+'BetaEtreeR.pdf.svg" alt="some_text">'

    outhtml+="<H2>Per Node Statistics in tabular format</H2>"
    nodeHeader=["betaA","pvalue","Significant","ag","agCount","betaCA","pvalueCA"]
    outhtml+='<table border="1"><tr><th>NodeName</th>'
    for h in nodeHeader:
        if h=="ag":
            h="ag"+str(db.groups.keys())
        outhtml+='<th>'+h+'</th>'

    outhtml+='</tr>'


#Put the correct proportion on each branches
R=db.GetEntropies(q=com["-q"])
def seebeta(n):
    return n.betaA
N=list(db.tree.traverse())
N.sort(key=seebeta, reverse=True)

if bool(int(com["-h"])):
    for n in N:
        outhtml+='<tr><td><a name="'+n.name+'">'+n.name+'</a></td>'
        for h in nodeHeader:
            if h=="Significant":
                outhtml+='<td>'+str(n in SN)+'</td>'
            elif h=="agCount":
                outhtml+='<td>'
                for tag,ag in zip(countbyGroup,n.__getattribute__("ag")):
                    #print tag, ag
                    outhtml+=str(int(tag*ag))+" "
                
                outhtml+='</td>'
            elif h=="ag":
                outhtml+='<td>'+", ".join([str(round(x,3)) for x in n.__getattribute__(h)])+'</td>'
            else:
                outhtml+='<td>'+str(n.__getattribute__(h))+'</td>'
        
        outhtml+='</tr>\n'

    outhtml+='</table>'

    ##for n in N:
    ##    print "["+",".join(map(str,[n.__getattribute__("betaA"),sum(n.__getattribute__("betaAkl"))]))+"],"

    outhtml+='</body>\n'
    handle=open(com["-o"]+".html","w")
    handle.write(outhtml)
    handle.close()

[layout_branch(x,"betaA") for x in db.tree.traverse()]
for i in SN:
    pass
    #i.img_style["bgcolor"]='#FFFF00'

if com["-x"]=="nexml":
    [Features2NeXML(x,features=["pvalue"]) for x in SN]
    [Features2NeXML(x,features=["betaA","betaCA"]) for x in db.tree.traverse()]
    nexml_project.export(open(com["-o"]+".nexml","w"))

if bool(int(com["-h"])):
    handle=open(com["-o"]+".TreeLabeled","w")
    handle.write(db.tree.write(format=1))
    handle.close()
    S=GraphFeatures2Itollabel(db.tree,SN)
    handle=open(com["-o"]+"_tableXitol.txt","w")
    handle.write(S)
    handle.close()
    SS=HistFeatures2Itol(db)
    handle=open(com["-o"]+"_tableHistXitol.txt","w")
    handle.write(SS)
    handle.close()
    #SSS=SignFeature2Itol(SN)
    #handle=open(com["-o"]+"_tableSignNodeXitol.txt","w")
    #handle.write(SSS)
    #handle.close()
    test = Itol.Itol()
    test.add_variable('treeFile',com["-o"]+".TreeLabeled")
    test.add_variable('treeName',com["-o"])
    test.add_variable('treeFormat',"newick")
    test.add_variable('showInternalIDs','1')
    test.add_variable('colorDefinitionFile',com["-o"]+"_tableXitol.txt")
    test.add_variable('dataset1File',com["-o"]+"_tableHistXitol.txt")
    test.add_variable('dataset1Label','Counts')
    test.add_variable('dataset1Separator','comma')
    test.add_variable('dataset1Type','multibar')
    #test.add_variable('branchLabelsFile',com["-o"]+"_tableSignNodeXitol.txt")
    test.print_variables()
    good_upload = test.upload()
    if good_upload == False:
        print 'There was an error:'+test.comm.upload_output
    link=test.get_webpage()
    itol_exporter = test.get_itol_export()
    itol_exporter.set_export_param_value('format', 'svg')
    itol_exporter.set_export_param_value('rangesCover','clades')
    itol_exporter.set_export_param_value('showInternalLabels','1')
    itol_exporter.set_export_param_value('colorBranches','1')
    itol_exporter.set_export_param_value('datasetList','dataset1')
    itol_exporter.export(com["-o"]+"BetaEtree.svg")
    import re
    import sys
    print "try to add link from tree to table and then add to html file"
    def addHref(matchobj):
        name=matchobj.group()[1:-1]
        return '><a href="#'+name+'">'+name+'</a><'
    f=open(com["-o"]+"BetaEtree.svg", "r").read()
    A=re.sub(">L[0-9]*<",addHref,f)
    link="<a href='"+link+"'>Click here to modify image</a>"
    A=link+A
    F=open(com["-o"]+".html", "r").read()
    B=re.sub('<img src=".+" alt="some_text">',A,F)
    F=open(com["-o"]+".html", "w")
    F.write(B)
    F.close()




if com["-x"]=="phyloxml":
    [mapWidthANDColor2Phyloxml(x) for x in db.tree.traverse()]
    project.export(outfile=open(com["-o"]+"BetaE_phyloxml.xml","w"))
    #following line make xml readable to archeopeteryx viewer
    def correctPhyloXML(name):
        handle=open(name,"r")
        A=handle.read()
        handle.close()
        import re
        B=re.sub('<phy:Phyloxml xmlns:phy="http://www.phyloxml.org/1.10/phyloxml.xsd">','<phyloxml xmlns:phy="http://www.phyloxml.org/1.10/phyloxml.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"  xmlns="http://www.phyloxml.org" >',A)
        C=re.sub('<phy:phylogeny>','<phy:phylogeny rooted="true">',B)
        D=re.sub('/phy:Phyloxml','/phyloxml',C)
        handle=open(name,"w")
        handle.write(D)
        handle.close()
    
    correctPhyloXML(com["-o"]+"BetaE_phyloxml.xml")




[layout_branch(x,"betaCA") for x in db.tree.traverse()]
if com.has_key("-qt"):
    db.tree.render(com["-o"]+"BetaSGivenEtree.pdf")

if com["-x"]=="phyloxml":
    [mapWidthANDColor2Phyloxml(x) for x in db.tree.traverse()]
    project.export(outfile=open(com["-o"]+"BetaSGivenE_phyloxml.xml","w"))
    correctPhyloXML(com["-o"]+"BetaSGivenE_phyloxml.xml")


[layout_branchB(x) for x in db.tree.traverse()]
if com.has_key("-qt"):
    db.tree.render(com["-o"]+"Difftree.pdf")


if com["-x"]=="phyloxml":
    [mapWidthANDColor2Phyloxml(x) for x in db.tree.traverse()]
    project.export(outfile=open(com["-o"]+"Diff_phyloxml.xml","w"))
    correctPhyloXML(com["-o"]+"Diff_phyloxml.xml")









