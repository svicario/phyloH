import Itol,ItolExport
import numpy
from numpy import exp,arange,log
from pandas import DataFrame, Series,MultiIndex, cut
import os
import sys
if len(sys.argv)>1:
    sep=os.sep
    pathScript=os.path.dirname(os.path.abspath(sys.argv[0]))
    pathAncillary=os.path.dirname(os.path.abspath(sys.argv[0]))+sep+"AncillaryFiles"+sep


def AddTaxonomy(H, db, taxonomyfile):
    #Define reference taxonomy
    handle=open(taxonomyfile,"r")
    TAX=handle.readlines()
    handle.close()
    TAX=[x.strip().split("|") for x in TAX]
    taxDB={}
    for i in TAX:
        speciesName=i[0]
        # species names in taxonomy have this modification.
        #speciesName=speciesName.strip().replace("sp.","sp").replace("_"," ").replace("-"," ")
        taxDB[speciesName]=i[1:]
    
    if taxDB.has_key(" "):
        del taxDB['']
    if taxDB.has_key("_"):
        del taxDB['_']
    #build node translation map
    NodeTaxonDB={}
    bad=[]
    check=[]
    def NodeName(node):
        #taxroot=taxDB[" ".join(l[0].split("_")[1:])]
        l=[x.name for x in node.get_terminals()]
        taxroot=None
        for i in l:
            try:
                X=taxDB[i]
                #X=taxDB[" ".join(i.split("_")[1:])]
                if not taxroot:
                    taxroot=X
                else:
                    taxroot=[y for x,y in zip(taxroot,X) if x==y ]
            except KeyError:
                bad.append(i)
        if not taxroot:
            return None
        # if len(taxroot)==0: 
        #   print l, taxroot, X, node,i.split("_")[1:]
        #   check.append(node)
        return taxroot[-1]
    
    def Traverse(node):
        NAME=NodeName(node)
        if NAME:
            NodeTaxonDB[node.name]=NAME
        else:
            NodeTaxonDB[node.name]=node.name
        for c in node.clades:
            #if not c.is_terminal():
            Traverse(c)
                #if NodeTaxonDB[c.name]==NodeTaxonDB[node.name]:
                #    NodeTaxonDB[c.name]=c.name
    
    Traverse(db.tree.clade)
    index,values=zip(*NodeTaxonDB.items())
    NodeTaxonDB=Series(values,index=index, name="Taxonomy")
    #NodeTaxonDB.index=MultiIndex.from_tuples(zip(index,len(index)*[False]))
    #NodeTaxonDB.index.names=["Name","Is_Leaf"]
    return NodeTaxonDB

def getSignNodeMult(H,TreeStat, alpha=0.05):
    Sign=[]
    Tempo=[]
    for n,k in enumerate(["ITEi", "ITSgivenEi"]):
        L=float(len(TreeStat[k].columns))
        #if H[k].__class__!=DataFrame().__class__:
        #    H[k]=DataFrame({"Bits":H[k]})
        temp=(TreeStat[k].values>=DataFrame(H[k].nats).values).sum(axis=1)
        temp=temp/L
        H[k]["pvalue"]=temp
        SGN=H[k].pvalue<alpha/(H[k].pvalue.rank(method="first")+1)
        H[k]["MultTest"]=SGN
        m=H[k].shape[0]
        if alpha*SGN.sum()/m < (1/L):
            if SGN.sum()>0:
                print str(L)+' replicates are too low to evaluate correctly significance, you should increase to at least '+ str(1+m/(alpha*SGN.sum()))
        #If you change here, remember to change also in QRtree
        Tempo+=[(H["MI"].iloc[n,0]<=TreeStat[k].sum(axis=0).values).sum()/L]
    H["MI"]["pvalue"]=Tempo

def DecorateH(H, db, alpha=0.05, taxonomy=None):
    #Calculating Pvalue and Turnover from MI
    H["HSgivenE"]=H["HS"]-H["HE"]
    H["MI_treeAndEnvironment"]=DataFrame([H["MI_treeAndEnvironment"]],columns=["nats"], index=["I(T,G)"])
    H["MI_treeAndSampleGivenEnvironment"]=DataFrame([H["MI_treeAndSampleGivenEnvironment"]],columns=["nats"], index=["I(T,S|G)"])
    H["ITEi"]=DataFrame(H["ITEi"],columns=["nats"])
    H["ITSgivenEi"]=DataFrame(H["ITSgivenEi"],columns=["nats"])
    namesTE=["MI_treeAndEnvironment","ITEi"]
    namesTSgivenE=["MI_treeAndSampleGivenEnvironment", "ITSgivenEi"]
    for n in namesTE+namesTSgivenE:
        if n in namesTE:
            H[n]["TurnOver"]=H[n].nats/H["HE"]
        
        if n in namesTSgivenE:
            H[n]["TurnOver"]=H[n].nats/H["HSgivenE"]
    #print H["MI_treeAndEnvironment"]
    H["MI"]=H["MI_treeAndEnvironment"].append(H["MI_treeAndSampleGivenEnvironment"])
    del H["MI_treeAndEnvironment"]
    del H["MI_treeAndSampleGivenEnvironment"]
    if db.TreeStat["ITEi"].values.shape[1]:
        getSignNodeMult(H, db.TreeStat, alpha=0.05)
    #Several line of makeup to make nice table
    #Formatting Counts
    Counts=H["counts"]
    Counts.index.name=""
    Counts.index=[""]
    #Counts=Counts.swaplevel("Group","Sample",axis=1)
    NEWColumns=Counts.columns.sort_values()
    Counts=Counts[NEWColumns]
    Levels=Counts.columns.levels
    Labels=Counts.columns.labels
    #Levels=[[H["tot"]]]+[list(Levels[0])]+[H["tag"].values[0]]+[list(Levels[1])]+[list(Counts.values[0,Labels[1]])]
    Levels=[[H["tot"]]]+[list(Levels[0])]+[H["tag"].values[0]]+[list(Levels[1])]+[list(Counts.values)]
    #print Labels
    #Labels=[len(Labels[0])*[0]]+[list(Labels[0])]+[list(Labels[0])]+[list(Labels[1])]+[list(Labels[1])]
    Labels=[len(Labels[0])*[0]]+[list(Labels[0])]+[list(Labels[0])]+[list(Labels[1])]+[range(Counts.shape[1])]
    #print Labels,Levels
    COL=MultiIndex(levels=Levels, labels=Labels, names=["Total Counts","Group Name","Group Counts","Sample Name","Sample Counts"])
    CCounts=DataFrame([Counts.shape[0]*[""]], index=COL,columns=[""])
    H["counts"]=CCounts
    glevel=H["Pie"].shape[1]
    slevel=CCounts.shape[0]
    del H["tag"]
    del H["tot"]
    #Formatting Gammas
    H['HgammaEachEnvironment']["Overall"]=H["Hgamma"]
    temp=DataFrame(H['HgammaEachEnvironment'], columns=["nats"])
    temp["Diversity"]=exp(temp)
    H["Gammas"]=temp
    del H['HgammaEachEnvironment']
    del H["Hgamma"]
    #Formatting Alphas
    H["Alphas"]=DataFrame([H['HalphaBySamples'], H['HalphaByEnvironment']], index=["H(T|S)", "H(T|G)"], columns=["nats"])
    H["Alphas"]["Diversity"]=exp(H["Alphas"].nats)
    del H['HalphaBySamples']
    del H['HalphaByEnvironment']
    #Formatting Experimental Design
    temp=DataFrame([H["HE"],H["HS"],H["HSgivenE"]],index=["H(G)","H(S)","H(S|G)"], columns=["nats"])
    temp["Diversity"]=exp(temp)
    temp["MaxDiversity"]=[glevel,slevel,slevel/float(glevel)]
    H["ExperimentalDesign"]=temp
    del H["HE"]
    del H["HS"]
    del H["HSgivenE"]
    #Formatting By Branch result
    H["ITEi"].columns=MultiIndex.from_tuples([("I(Ti,G)",x) for x in list(H["ITEi"].columns)])
    H["ITEi"].columns.names=["Metric","Stat"]
    H["ITSgivenEi"].columns=MultiIndex.from_tuples([("I(Ti,S|G)",x) for x in list(H["ITSgivenEi"].columns)])
    H["ITSgivenEi"].columns.names=["Metric","Stat"]
    H["Pie"]=H["Pie"].fillna(0)
    temp=zip(["By Group Relative Frequency"]*len(H["Pie"].columns),list(H["Pie"].columns))
    H["Pie"].columns=MultiIndex.from_tuples(temp)
    H["Pie"].columns.names=["Metric","Stat"]
    temp=H["ITEi"].join(H["Pie"])
    temp=temp.join( H["ITSgivenEi"])
    H["MIByBranch"]=temp.fillna(0)
    del H["ITSgivenEi"]
    del H["ITEi"]
    del H["Pie"]
    #adding taxonomy annotation
    if taxonomy:
        NodeTaxonDB=AddTaxonomy(H, db,taxonomy)
        NodeTaxonDB=NodeTaxonDB[H["MIByBranch"].index.get_level_values("Name")]
        #print "Adding Taxonomy"
        #print H["MIByBranch"].index
        #print "sum"
        #print H["MIByBranch"].sum()
        H["MIByBranch"].set_index(keys=NodeTaxonDB, append=True, inplace=True)
        #print H["MIByBranch"].sum()
        #print H["MIByBranch"].index
        H["MIByBranch"].reorder_levels(["Taxonomy","Name","Is_Leaf"], axis=0)
        #print H["MIByBranch"]
    else:
        H["MIByBranch"].set_index(keys=Series(["Unknown"]*H["MIByBranch"].shape[0],name="Taxonomy"), append=True, inplace=True)
        #print H["MIByBranch"]
    H["MI_KL"]=DataFrame(H["MI_KL"], columns=["KullBack-Lieber(PG(i)||Ptot(i)"])
    H["MI_KL"]["pvalue"]=numpy.sum(H["KL_perm"].values.transpose()>H["MI_KL"].values,axis=1)/float(H["KL_perm"].shape[0])
    SGN=H["MI_KL"].pvalue<alpha/(H["MI_KL"].pvalue.rank(method="first")+1)
    H["MI_KL"]["Seq_Bonferroni"]=SGN
    H["MIByBranch"].sort_values(by=[("I(Ti,G)","TurnOver")], axis="index",inplace=True,ascending=False)

def spacedColors(cat):
    step=360.0/cat
    h=[step*x for x in range(cat)]
    colors=[tohex(*hsv2rgb(i,s=1,v=1)) for i in h]
    return colors

def tohex(r,g,b):
    hexchars = "0123456789ABCDEF"
    return "#" + hexchars[r / 16] + hexchars[r % 16] + hexchars[g / 16] + hexchars[g % 16] + hexchars[b / 16] + hexchars[b % 16]

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

def ForITOL(H):
    NBINS=10
    values, bins=cut(H["MIByBranch"]["I(Ti,G)"].TurnOver,bins=NBINS, retbins=True)
    try:
        from matplotlib import pyplot as plt
        import matplotlib
    except ImportError:
        if NBINS==10:
            zz=Series(["#FFF1A9", "#FEE187", "#FECA66", "#FEAB49", "#FD8C3C","#FC5B2E","#ED2E21", "#D41020", "#B00026", "#800026"], index=values.cat.categories)
        else: raise ImportError
    else:
        cm = plt.get_cmap('YlOrRd')
        z=arange(1,(NBINS+1),1)/float(NBINS)
        zz=Series([matplotlib.colors.rgb2hex(x).upper() for x in cm(z)], index=values.cat.categories)
    XITOL=DataFrame({"branch name":list(values.index.get_level_values("Name")), "mode":"range","label":list(values.values), "color":zz[values]})
    H["MIByBranch"].loc[:,("I(Ti,G)","Color")]=zz[values].values
    #print "wwww"
    #print H["MIByBranch"]["I(Ti,G)"]
    #print H["MIByBranch"].columns
    #H["MIByBranch"].reindex(H["MIByBranch"].index)
    L=len(H["MIByBranch"].columns)
    H["MIByBranch"]=H["MIByBranch"].iloc[:,sum([range(4),[L-1],range(4,L-1)],[])]
    #H["MIByBranch"]=H["MIByBranch"].iloc[:,[0,1,2,3,11,4,5,6,7,8,9,10]]
    #print "CIAO"
    #print H["MIByBranch"].columns
    XITOL.set_index("branch name",inplace=True)
    XITOL["label"]=["_to_".join(x.split(", "))[1:-1] for x in XITOL["label"]]
    #print XITOL.iloc[0:3,:]
    label=numpy.array(["NotSignificant","Significant"])[(H["MIByBranch"]["I(Ti,G)"].MultTest*1).values]
    color=numpy.array(["#000000","#00FFFF"])[(H["MIByBranch"]["I(Ti,G)"].MultTest*1).values]
    XITOLbis=DataFrame({"mode":"clade","label":label, "color":color}, index=list(values.index.get_level_values("Name")) )
    XITOLbis.name="branch name"
    XITOL=XITOL.append(XITOLbis)
    XITOL=XITOL[[  "mode", "color","label"]]
    Pie=H["MIByBranch"]["By Group Relative Frequency"].query("Is_Leaf==True")
    color=spacedColors(Pie.shape[1])
    Pie.index=Pie.index.get_level_values("Name")
    Pie.columns=MultiIndex(levels=[[Pie.columns],[color]],labels=[range(Pie.shape[1])]*2, names=["LABELS","COLORS"])
    Pie.index.name=""
    #Transform in integer to do not upset ITOL
    HIST=(Pie*H["counts"].index.get_level_values("Total Counts")[0]).astype(int)
    return XITOL, HIST,H,values.cat.categories.tolist()
def secondaryOutput(H, db, com):
    DESC="\n".join([x[0]+"\t"+" ".join(x[1]) for x in db.TreeSummary[1]])
    handle=open(com["-o"]+"_descedant.csv","w")
    handle.write(DESC)
    handle.close()
    for k in H:
        temp=H[k].to_csv(sep="\t")
        if (k=="MIByBranch"):
            handle=open(com["-o"]+".mibybranch","w")
        else:
            handle=open(com["-o"]+"_"+k+".csv","w")
        handle.write(temp)
        handle.close()
def MakeHTML(H,com, errormessage=None):
    Titles={
        "counts":["Experimental Design:","Counts of observations across groups and samples within groups"],
        "ExperimentalDesign":["Entropy across samples, groups and samples within groups",
                              "MaxDiversity gives the values for a maximally balanced experimental design"],
        "Gammas":["Gamma diversities:",
                  "Total entropy and diversity within each group and overall data",
                  "Unit measure for Diversity is equivalent number of independent equi-abundant linneages"],
        'Alphas':["Alpha diversities:",
                  "Mean entropy and diversity within sample or group",
                  "Unit measure for Diversity is equivalent number of independent equi-abundant linneages"],
        "MI":["Beta diversity:",
              "Information shared across Tree and Sample or Group vector expressed  as nats and turnover of linneage.",
              "Turnover is the percentage of observations not shared across groups",
              "Pvalue is computed with Permutation procedure"],
        'MI_KL':["Difference of each group from total:",
                 "phylogenetic Kullback-Leiber distance between each group and the overall data", "P Value and Seq_Bonferroni are the probability that the difference is caused by sampling fluctuation and significance after Sequential Bonferroni correction and alpha=0.05, respectevely"],
        'DistTurnover':["Pairwise TurnOver between groups"]    
    }
    def floater(x, perc=False):
        from math import floor, log10
        if numpy.abs(x)==numpy.inf:
			return "Inf"
        if perc:
            x=100*x
        if (x/2 == float(x)/2) and x>0:
            temp= str(round(x, 2-int(floor(log10(x)))))
        else:
            temp= str(x)
        if perc:
            temp+="%"
        return temp
    def floaterPerc(x):
        return floater(x,perc=True)
    def Tagify(text, tag):
        return "<"+tag+">"+text+"</"+tag+">"
    htmltext=open(pathAncillary+sep+"Template.html","r").read()
    GlobalStatistics=""
    GlobalStatistics+=Tagify("The run call was:","H3")+Tagify(com["call"],"p")+"\n"
    if errormessage:
        GlobalStatistics+=Tagify("Errors", "H3")
        for e in errormessage:
            GlobalStatistics+=Tagify(e, "p")
    for k in [x for x in ['counts', 'ExperimentalDesign', 'Gammas','Alphas', 'MI', 'MI_KL','DistTurnover'] if x in H]:
        GlobalStatistics+=Tagify(Titles[k].pop(0),"H2")+"\n"
        for title in Titles[k]:
            GlobalStatistics+=Tagify(title,"p")+"\n"
        print k
        if k!="DistTurnover":
            GlobalStatistics+=H[k].to_html(float_format=floater,na_rep="",formatters={"TurnOver": lambda x: floater(x,perc=True)})
        else:
            GlobalStatistics+=H[k].to_html(float_format=floaterPerc,na_rep="")
    GlobalStatistics+="\n"
    temp=zip(*H["MIByBranch"].index.tolist())
    temp[0]=['<a name="'+x+'">'+x+'</a>' for x in temp[0]]
    tempname=H["MIByBranch"].index.names
    OUT=H["MIByBranch"]+0
    OUT.index=MultiIndex.from_tuples(zip(*temp))
    OUT.index.names=tempname
    BranchStatistics=OUT.to_html(float_format=floater, formatters={("I(Ti,G)","TurnOver"): lambda x: floater(x,perc=True),
        ("I(Ti,S|G)","TurnOver"): lambda x: floater(x,perc=True)}, escape=False)
    
    htmltext=htmltext.format(BranchStatistics=BranchStatistics, GlobalStatistics=GlobalStatistics,nameRun=com["-o"] )
    handle=open(com["-o"]+".html","w")
    handle.write(htmltext)
    handle.close()
    return "bau"
    

def makeITOL3output(tree, XITOL, HIST, bins, com):
    handle=open(com["-o"]+".TreeLabeled","w")
    handle.write(tree.format("newick"))
    handle.close()
    histtxt=open(pathAncillary+"ITOLStackedHistogram.txt","r").read()
    formatdict={
        "Data":HIST.to_csv(header=False,sep=","),
        "NomeRun":com["-o"],
        "NumeroUnoRipetutoQuantiICampioni":",".join(["1"]*HIST.shape[1]),
        "ColoriDelleCategorie": ",".join(HIST.columns.get_level_values("COLORS").tolist()).strip(),
        "NomiDelleCategorie": ",".join(HIST.columns.get_level_values("LABELS").tolist())
    }
    histtxt=histtxt.format(**formatdict)
    handle=open(com["-o"]+"_tableHistXitol.txt","w")
    handle.write(histtxt)
    handle.close()
    Size={}
    bins=["_to_".join(x.split(", "))[1:-1] for x in bins]
    Size.update([(y,str(log(x+1)*51/log(10))) for x,y in enumerate(bins)])
    Circletxt=open(pathAncillary+"ITOLInternalNodeCircle.txt","r").read()
    assert len(bins)==10, "Bins are expected to be 10, have you changed the number of the colors?"
    temp=XITOL[XITOL.loc[:,"mode"]=="range"]
    SS=zip(temp.index.tolist(),["2"]*temp.shape[0],[Size[x] for x in temp.label], temp.color, ["1"]*temp.shape[0],["1"]*temp.shape[0],temp.label)
    formatdict={
        "Data":"\n".join([",".join(x) for x in SS]),
        "CUTBins":",".join(bins)
    }
    Circletxt=Circletxt.format(**formatdict)
    handle=open(com["-o"]+"_tableCircleXitol.txt","w")
    handle.write(Circletxt)
    handle.close()
    import errno
    import shutil
    from shutil import ignore_patterns
    def copyD(src, dest):
        try:
            shutil.copytree(src, dest, ignore=ignore_patterns('.py', 'specificfile.file'))
        except OSError as e:
            # If the error was caused because the source wasn't a directory
            if e.errno == errno.ENOTDIR:
                shutil.copy(src, dest)
            else:
                print('Directory not copied. Error: %s' % e)
    copyD(pathAncillary+sep+"external",os.getcwd()+sep+"external")
    copyD(pathAncillary+sep+"PhyloCommunity.js",os.getcwd()+sep+"PhyloCommunity.js")
     
    
    
def makeITOLcall(tree,XITOL, HIST, com):
    handle=open(com["-o"]+".TreeLabeled","w")
    handle.write(tree.format("newick"))
    handle.close()
    buffITOL=XITOL.to_csv(header=False, sep="\t")
    handle=open(com["-o"]+"_tableXitol.txt","w")
    handle.write(buffITOL)
    handle.close()
    buffHIST=HIST.to_csv()
    handle=open(com["-o"]+"_tableHistXitol.txt","w")
    handle.write(buffHIST)
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
    test.add_variable('treeFormat',"newick")
    #test.add_variable('preCollapsedFile',com["-o"]+"_CollapseList.txt")
    test.add_variable('showInternalIDs','1')
    test.add_variable('colorDefinitionFile',com["-o"]+"_tableXitol.txt")
    test.add_variable('dataset1File',com["-o"]+"_tableHistXitol.txt")
    test.add_variable('dataset1Label','Counts')
    test.add_variable('dataset1Separator','comma')
    test.add_variable('dataset1Type','multibar')
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
    #f=re.sub('fill="black">[0-9_. A-z-]+</text>',addOpacity,f)
    f=re.sub('<path fill="white" .+\r',"",f)
    link="<a href='"+link+"'>Click here to modify image</a>"
    B=link+f
    F=open(com["-o"]+".html", "r").read()
    B=re.sub('<img src=".+" alt="some_text">',B,F)
    F=open(com["-o"]+".html", "w")
    F.write(B)
    F.close()
