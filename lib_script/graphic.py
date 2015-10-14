import Itol,ItolExport
import numpy
from numpy import exp
from entropyPerm5 import Diversity2Perc


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
    node._img_style=style

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
    #del node.img_style
    node._img_style=style


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
        out.append("\t".join(map(str,[n.name,"range",n._img_style["hz_line_color"]," - ".join(D[n._img_style["hz_line_width"]])])))
        if n in SN:
            out.append("\t".join([n.name,"clade","#00FFFF","Signif"]))
        else:
            out.append("\t".join([n.name,"clade","#000000","NotSignif"]))
    
    S="\n".join(out)
    return S

def HistFeatures2Itol(db):
    tree_obj=db.tree
    table=db.compressTable()
    tag=[numpy.sum(db.countTable[:,db.groups[x]]) for x in db.groups]
    abbundances=[]
    for i in tree_obj.get_leaves():
        if numpy.sum(i.ag)>0:
            a=[i.name]
            a+=list(map(int,i.ag*tag))
            abbundances.append(a)
    # abbundances=[]
    # for leaf in table.keys():
    #     a=[leaf]
    #     for g in db.groups:
    #         a.append(sum(table[leaf][db.groups[g]]))
    #     abbundances.append(a)
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

def makeHMTLDistMatrix(outdict,pattern, db):
    Gdict=outdict[pattern]
    l=len(db.groups)
    outhtml="<table border='1'><tr><th></th>"
    for gg in db.groups:
        outhtml+="<th>"+gg+"</th>"
    outhtml+="</tr>\n"
    for g in db.groups:
        outhtml+="<tr><th>"+g+"</th>"
        temp=[x for x in Gdict if x[0]==g]
        for gg in db.groups:
            if (g,gg) in temp:
                outhtml+="<td>"+str(Gdict[(g,gg)])+"</td>"
            else:
                outhtml+="<td></td>"
        outhtml+="</tr>"
    outhtml+="</table border='1'>"
    return outhtml
        
    
    
    
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

def makeXMLoutput(HS,HE,db,Names,P,countbyS, R, com, Gdict, HSgivenE):
    OUT="<res>\n"
    OUT+="<HSample>\n\t<Observed>"+str(HS)+"</Observed>\n\t<Diversity>"+str(exp(HS))+"</Diversity>\n</HSample>\n"
    OUT+="<HEnvironment>\n\t<Observed>"+str(HE)+"</Observed>\n\t<Diversity>"+str(exp(HE))+"</Diversity>\n</HEnvironment>\n"
    outdict={"HEnvironment":{"Observed":HE,"Diversity":exp(HE),"MaxDiversity":len(db.groups.keys())},"HSample":{"Observed":HS,"Diversity":exp(HS),"MaxDiversity":len(countbyS)}}
    count=0
    for label, Pvalues in zip(Names,P):
        outdict[label]={}
        if label.find("HgammaFor_")>-1:
            Group=label.split("For_")[1]
            OUT+="<Hgamma group='"+Group+"'>\n"
            Nlabel="HgammaEachEnvironment"
            OUT+="\t<Observed>"+str(R[Nlabel][count])+"</Observed>\n"
            outdict[label]["Observed"]=R[Nlabel][count]
            D=exp(R[Nlabel][count])
            if label.find("alpha")>-1: name="alpha"
            elif label.find("gamma")>-1: name="gamma"
            OUT+="\t<Diversity>"+str(D)+"</Diversity>\n"
            outdict[label][name+"_Diversity"]=D
            count+=1
            OUT+="</Hgamma>\n"
        elif label not in ["MI_KL"]:
            OUT+="<"+label+">\n"
            OUT+="\t<Observed>"+str(R[label])+"</Observed>\n"
            outdict[label]["Observed"]=R[label]
            if label.find("MI_tree")>-1:
                outdict[label]["Pvalue"]=sum(R[label]<=Pvalues)/float(len(Pvalues))
                OUT+="\t<Pvalue>"+str(sum(R[label]<=Pvalues)/float(len(Pvalues)))+"</Pvalue>\n"
                D=exp(R[label])
                OUT+="\t<Diversity>"+str(D)+"</Diversity>\n"
                outdict[label]["beta_Diversity"]=D
                if label.find("treeAndEnvironment")>-1:
                    Cardinality=exp(HE)
                elif label.find("SampleGivenEnvironment")>-1:
                    Cardinality=exp(HSgivenE)
                S=Diversity2Perc(D,Cardinality)
                OUT+="\t<Percentage_Overlap>"+str(S)+"</Percentage_Overlap>\n"
                outdict[label]["Percentage_Overlap"]=S
            elif (label not in ["tag","tot"]) :
                D=exp(R[label])
                OUT+="\t<Diversity>"+str(D)+"</Diversity>\n"
                outdict[label]["alpha_Diversity"]=D
            OUT+="</"+label+">\n"

    R["MI_KL"]=list(R["MI_KL"])
    for temp, g in zip(R["MI_KL"],db.groups.keys()):
        outdict["KL_of_"+g]={"Observed":temp}
        OUT+="<MI_KL group='"+g+"'>\n"
        OUT+="\t<Observed>"+str(temp)+"</Observed>\n"
        D=exp(temp)
        OUT+="\t<Diversity>"+str(D)+"</Diversity>\n"
        Cardinality=exp(HE)
        S=Diversity2Perc(D,Cardinality)
        OUT+="\t<Percentage_Overlap>"+str(S)+"</Percentage_Overlap>\n</MI_KL>\n"
        
    OUT+="<PairwiseTurnover>\n"
    for g, gg in Gdict:
        OUT+="<Turnover pair='"+g+"_"+gg+"'>"+str(Gdict[(g,gg)])+"</Turnover>"
    OUT+="<Counts><Total>"+str(sum(countbyS))+"</Total></Counts>\n"
    outdict["PairwiseTurnover"]=Gdict
    nomi=numpy.array(db.samplesNames)
    for g in db.groups:
        OUT+="<Samples group='"+g+"'>"+",".join(list(nomi[db.groups[g]]))+"</Samples>\n"
    
    OUT+="</res>"
    handle=open(com["-o"]+".xml","w")
    handle.write(OUT)
    handle.close()
    return outdict

def makeHMTLoutput(db,countbyS, HS, HE, R, com, outdict,call):
    Venn="""
    <svg width="600" height="500" x="0" y="0" xmlns="http://www.w3.org/2000/svg">
        <title>Report of PhyloH with use  of BioVenn and Itol graphics</title>
        <desc></desc>
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
    CountHTML='<table border=1>'
    CountHTML+='<tr><th colspan='+str(len(countbyS))+'>total</th></tr>'
    CountHTML+='<tr><td align="center" colspan='+str(len(countbyS))+'>'+str(sum(countbyS))+'</td></tr>'
    CountHTML+="</tr>"
    for k in db.groups.keys():
            CountHTML+='<th colspan='+str(len(countbyS[db.groups[k]]))+'>'+k+'</th>'
    CountHTML+="</tr><tr>"
    for k in db.groups.keys():
            CountHTML+='<td align="center" colspan='+str(len(countbyS[db.groups[k]]))+'>'+str(sum(countbyS[db.groups[k]]))+'</td>'
    CountHTML+="</tr><tr>"
    for k in db.groups:
        for s in db.groups[k]:
            CountHTML+='<td align="center">'+str(countbyS[s])+"</td>"
    CountHTML+="</tr></table>"                                                                         

    translator={"HEnvironment":"H(E)","HSample":"H(S)","Hgamma":"H(T)","HalphaByEnvironment":"H(T|E)","HalphaBySamples":"H(T|S)","MI_treeAndEnvironment":"I(T,E)","MI_treeAndSampleGivenEnvironment":"I(T,S|E)"}
    for g in db.groups.keys():
        translator["HgammaFor_"+g]="H(T|E="+g+")"
    
    outhtml="<!DOCTYPE html><html><body><H1>Partitioning information in  "+com["-o"] +"</H1><p>analysis called:"+call+"</p><H2>Global Statistics</H2>"

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
    outhtml+="""<H3>Beta Diversity: Pairwise TurnOver </H3>"""
    outhtml+=makeHMTLDistMatrix(outdict, "PairwiseTurnover",db)
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
    outhtml+='<img src="'+com["-o"]+'BetaEtree.svg" alt="some_text">'

    return outhtml

def makeHMTLTree(N,outhtml, NodeTaxonDB,db, SN, countbyGroup, com, tabPrune):
    headnode=[x.split("\t")[0] for x in tabPrune.split("\n")]
    outhtml+="<H2>Per Node Statistics in tabular format</H2>"
    nodeHeader=["betaA","pvalue","Significant","ag","agCount","betaCA","pvalueCA"]
    outhtml+='<table border="1"><tr><th>TaxonName</th><th>NodeName</th>'
    for h in nodeHeader:
        if h in ["ag", "agCount"]:
            h+=str(db.groups.keys())
        outhtml+='<th>'+h+'</th>'
    outhtml+='</tr>'
    for n in N:
        try:
            tx=NodeTaxonDB[n.name]
        except KeyError:
            tx=""
        if n.name in headnode:
            outhtml+='<tr><td><a href="#'+n.name+'Col" id="'+n.name+'">'+tx+'</a></td><td>'+n.name+'</td>'
        else:
            outhtml+='<tr><td><a id="'+n.name+'">'+tx+'</a></td><td>'+n.name+'</td>'
        for h in nodeHeader:
            if h=="Significant":
                outhtml+='<td>'+str(n in SN)+'</td>'
            elif h=="agCount":
                outhtml+='<td>'
                #outhtml+=",".join(zip(countbyGroup,n.__getattribute__("ag")))
                countLinneage=[]
                for tag,ag in zip(countbyGroup,n.__getattribute__("ag")):
                    countLinneage.append(str(int(tag*ag)))
                    #outhtml+=str(int(tag*ag))+", "
                outhtml+=",".join(countLinneage)
                outhtml+='</td>'
            elif h=="ag":
                outhtml+='<td>'+", ".join([str(round(x,3)) for x in n.__getattribute__(h)])+'</td>'
            else:
                outhtml+='<td>'+str(n.__getattribute__(h))+'</td>'
        outhtml+='</tr>\n'
    outhtml+='</table>'
    outhtml+="<H2>Collapsed sequences</H2>"
    outhtml+='<table border="1"><tr><th>NodeName</th><th>Collapsed Sequences</th>'
    for i in tabPrune.split("\n")[1:]:
        node,sequences=i.split("\t")
        outhtml+='<tr><td><a id="'+node+'Col">'+node+'</a></td><td>'+sequences+'</td>'
    ##for n in N:
    ##    print "["+",".join(map(str,[n.__getattribute__("betaA"),sum(n.__getattribute__("betaAkl"))]))+"],"
    outhtml+='</body>\n'
    handle=open(com["-o"]+".html","w")
    handle.write(outhtml)
    handle.close()

def makeITOLcall(NodeTaxonDB,com,SN, db):
    if NodeTaxonDB:
        taxo="\n".join(["\t".join([x,x]) for x in NodeTaxonDB.keys()])
        handle=open(com["-o"]+"_taxonLabel.txt","w")
        handle.write(taxo)
        handle.close()
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
    #handle=open(com["-o"]+"_CollapseList.txt","w")
    #handle.write(Collapse)
    #handle.close()
    #SSS=SignFeature2Itol(SN)
    #handle=open(com["-o"]+"_tableSignNodeXitol.txt","w")
    #handle.write(SSS)
    #handle.close()
    test = Itol.Itol()
    test.add_variable('treeFile',com["-o"]+".TreeLabeled")
    test.add_variable('treeName',com["-o"])
    if NodeTaxonDB:
        test.add_variable('branchLabelsFile',com["-o"]+"_taxonLabel.txt")
    test.add_variable('treeFormat',"newick")
    #test.add_variable('preCollapsedFile',com["-o"]+"_CollapseList.txt")
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
    SCRIPT="""
            <style>
                text:hover
                {
                    opacity: 1;
                }
            </style>
            """
    def addHref(matchobj):
        name=matchobj.group()[1:-1]
        return '><a xlink:href="#'+name+'">'+name+'</a><'
    def addOpacity(matchobj):
        name=matchobj.group()[:]
        return "opacity=0 "+name
    def addStyle(matchobj):
        name=matchobj.group()[:]
        return name+SCRIPT
    f=open(com["-o"]+"BetaEtree.svg", "r").read()
    f=re.sub("</style>",addStyle,f)
    f=re.sub(">L[0-9]*<",addHref,f)
    f=re.sub('fill="black">[0-9_. A-z-]+</text>',addOpacity,f)
    f=re.sub('<path fill="white" .+\r',"",f)
    link="<a href='"+link+"'>Click here to modify image</a>"
    B=link+f
    F=open(com["-o"]+".html", "r").read()
    B=re.sub('<img src=".+" alt="some_text">',B,F)
    F=open(com["-o"]+".html", "w")
    F.write(B)
    F.close()