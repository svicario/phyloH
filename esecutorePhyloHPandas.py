from lib_script.entropyPerm5 import *
from lib_script.Decorators import DecorateH, ForITOL, MakeHTML,makeITOLcall,secondaryOutput, makeITOL3output
#from lib_script.graphic import *
from StringIO import StringIO
import copy

def BlackLister(db, toTrash="SampleAndTree", Suffix="QUERY__"):
    """
    make a blacklist
    """
    NameInTree=set([x.name for x in db.tree.get_terminals()])
    #to be changed when all will be pandas
    NameInSample=set(db.indexID.keys())
    if toTrash=="TreeOnly":
        blacklist=NameInTree.difference(NameInSample)
    elif toTrash=="SampleAndTree":
        blacklist=NameInSample.difference(NameInTree)
    elif toTrash=="Suffix":
        blacklist=[x for x in NameInTree+NameInSample if x.find(Suffix)==0]
    else:
        raise ValueError, "toTrash value not among the 3 possible: TreeOnly, SampleAndTree,Suffix"
    return blacklist

def collapseQR(node, blacklist,recipient=None):
    """
    Collapse branch with only leaf with in blacklist
    """
    for i in node.clades:
        if i.is_terminal():
            if i.name in blacklist:
                #recipient.children.append(queryOne)
                if noSingleton:
                    pass
        #if all([x.find(prefix)==0 for x in i.get_leaf_names()]):
        termi=[x.name for x in i.get_terminals()]
        if set(termi).difference(blacklist)==set([]):
            if recipient:
                    recipient.clades.append((i.name,termi))
            for c in i.clades:
                i.collapse(c)
        else:
                collapseQR(i,blacklist,recipient)

def MakeQRdata(db, blacklist):
    """
    function that modify sample an and group so to observe difference
    between observed sample ( with name that start with "QUERY") and
    reference one used only to build the tree
    Reference sequence 
    """
    totNames=set([x.name for x in db.tree.get_terminals()])
    Qlist=blacklist
    Rlist=totNames.difference(blacklist)
    OUT=""
    for r in Rlist:
        OUT+="\t".join(["R","1",r])+"\n"
    for q in Qlist:
        OUT+="\t".join(["Q","1",q])+"\n"
    #handle=open("QvsR.sample","w")
    handle=StringIO()
    handle.write(OUT)
    handle.pos=0
    return handle



def QRtree(db, H):
    """
    function that remove the contribution to mutual information caused
    by an excess of Reference observation to Query observation.
    So that result are only caused by query that do not find match in
    reference.
    Input are the tree and R the dictionary that collect result from
    tree calculation.
    Tree and R are modified and no explicit result is given.
    """
    betaAplus=0
    Pie=H["MIByBranch"]["By Group Relative Frequency"]
    Pie.Q
    H["MIByBranch"]=temp=H["MIByBranch"][Pie.Q<Pie.R]
    H["MI"].loc["I(T,E)"]["nats","TurnOver"]= temp["I(Ti,E)"]["nats","TurnOver"].sum()
    H["MI"].loc["I(T,S|E)"]["nats","TurnOver"]= temp["I(Ti,S|E)"]["nats","TurnOver"].sum()

if __name__=="__main__":
    com={"-x":"nexml","None":1, "--QR":"0","--QRC":"0","-k":0, "-e":0,"-q":1,"-G":0,"-H":"0","-g":"locationID","-M":"0"}
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
    spiegazione=""" 
     -f filename        Use this file as the phylogeny file [phylo].
     -o filename        Use this file to record output
     -t filename        Use this file to get taxonomic information on some tips and map them in internal node
     -s filename        Use this file as the sample file [sample].If "-G 1" the file is a CSV that follow darwin archive standard and heading following TWDG
     -g filename/string Use this file as the group file [group]/ if "-G 1" is a string that indicates the column in the CSV where to retrive groupping.
     -r INT             Number of randomizations to use [999]
     -x string          Two possible strings :"nexml" or "phyloxml" to select the xml output of the results
     -h 0 or 1          Boolean to check if you want html output
     --QR  0 or 1       Identify linneage present in the Query but not in the Reference (Need found observation to be tagged with "Query" prefix)
     --QRC 0 or 1       Collapse branch with only query before analysis
     -k    0 or 1       Perform pairwise comparison among all sample to explore variation across all data
     -e    0 or 1       Assume that each sampling site have the same sampling effort. This option equalize the weigth of the different sample
     --treesimplify     Collapse after analysis all descedant nodes that have weighted length to tips less that given [0.01]
     -G    0 or 1       Geographic analysis mode: group and sample and potentially taxonomy  and tree are all embedded in a CSV that includes also geographical location of observations.
     -H    INT          If data are Geographic, group observation by hexagonal grid of a given span in meters, if zero no gridding is applied and locationID parameters is expected
     -M    string       Cross obseration with a shapecollection  -G need to be 1 and -H need to be 0
    """
    if (( (com["--QR"]=="0") and (not '-s' in com) ) | ((not '-f' in com ) and  (com["-G"]=="0"))):
        print spiegazione
        raise ImportError
    com["call"]=" ".join(sys.argv)
    print "sono qui"
    if com["-G"]=="1":
        from lib_script.geoAddOn import *
        if int(com["-H"])==0:
            print "H0"
            if com["-M"]=="TNC":
                if com["-g"].lower() in ["biomes", "major_habitat_type"]: com["-g"]="WWF_MHTNAM"
                if com["-g"].lower() in ["biogeographic_region", "realm"]: com["-g"]="WWF_REALM2"
                ShapeLoader(shapeName="terr-ecoregions-TNC/tnc_terr_ecoregions",samplefilename=com["-s"],pathPythonScript=com["call"].split(" ")[0], group=com["-g"])
                IN=os.path.dirname(os.path.abspath(com["-s"]))+os.path.sep+"grid_"+os.path.basename(com["-s"])
                sample,groupBy=makePhyloHInput(csv=IN,sep=",", groupBy=com["-g"],  shape=False, makePhylo=("-f" not in com),  makeTaxo=("-t" not in com))
                print ("Sample", sample, groupBy)
                com["-g"]="Group_"+com["-g"]
                com["-s"]="Sample"
                com["sample"]=sample
                com["groupBy"]=groupBy
                print "option TNC with GRASS"
            else:
                sample,groupBy=makePhyloHInput(csv=com["-s"],sep=",", groupBy=com["-g"],  shape=False, makePhylo=("-f" not in com),  makeTaxo=("-t" not in com))
                com["-g"]="Group_"+ groupBy
                com["sample"]=sample
                com["groupBy"]=groupBy
                IN=com["-s"]
                print "option CSV"
        else:
            GridMaker(com["-s"],com["-H"], com["call"].split(" ")[0])
            IN=os.path.dirname(os.path.abspath(com["-s"]))+os.path.sep+"grid_"+os.path.basename(com["-s"])
            sample,groupBy=makePhyloHInput(csv=IN,sep=",", groupBy="cellid",  shape=False, makePhylo=("-f" not in com),  makeTaxo=("-t" not in com))
            com["-g"]="Group_cellid"
            com["sample"]=sample
            com["groupBy"]=groupBy
            print "option GRID with GRASS"
        namefile=".".join(IN.split(".")[:-1])
        com["-s"]="Sample"
        if ("-t" not in com):
            com["-t"]=namefile+".taxonomy"
        if not com.has_key("-f"):
            com["-f"]=namefile+".tree"
    print com
    #print sys.argv[0]
    #print os.path.abspath(sys.argv[0])
    db=DBdata()
    db.readTreePandas(com['-f'])
    if com["--QR"]=="0":
        db.readSampleTable(com["-s"])
        if com.has_key("-g"):
            db.readGroupTable(com["-g"])
        else:
            db.noGroup()
    else:
        if com["--QRC"]=="1":
            blacklist=BlackLister(toTrash="SampleAndTree")
            QueryTreeSegment=[]
            collapse(db.tree,blacklist,QueryTreeSegment)
            QueryTreeSegment="\n".join([x[0]+"\t"+" ".join(x[1]) for x in QueryTreeSegment])
            if QueryTreeSegment:
                QueryTreeSegment="Kept\tCollapsed\n"+QueryTreeSegment
                handle=open("PreliminaryCol.tab","w")
                handle.write(QueryTreeSegment)
                handle.close
        sample=QRdata(db)
        db.readSampleTable(sample)
        db.noGroup()
    
    Per=db.PermutationTestPandas(reps=int(com["-r"]),q=com["-q"])
    H=db.GetEntropiesPandas(q=com["-q"], Pairwise=int(com["-k"]), EqualEffort=int(com["-e"]))
    H["KL_perm"]=DataFrame([x["MI_KL"] for x in Per])
    DecorateH(H, db, taxonomy=com.setdefault("-t",None))
    if com["--QR"]==1:
        blacklist=BlackLister(toTrash="SampleAndTree")
        QueryTreeSegment=[]
        collapse(db.tree,blacklist,QueryTreeSegment)
        QueryTreeSegment="\n".join([x[0]+"\t"+" ".join(x[1]) for x in QueryTreeSegment])
        if QueryTreeSegment:
            QueryTreeSegment="Kept\tCollapsed\n"+QueryTreeSegment
            handle=open("PreliminaryCol.tab","w")
            handle.write(QueryTreeSegment)
            handle.close
    #Start ouput
    for k in ['counts', 'ExperimentalDesign', 'Gammas','Alphas', 'MI', 'MI_KL','DistTurnover']:
        pass
        #print H[k]
    errormessage=[]
    diffName=set(db.SeqName).difference([x.name for x in db.tree.get_terminals()])
    if diffName:
        errormessage.append("\n".join(["Following OTUs are not present in the Tree:",",".join(list(diffName))]))
        print "Errore"
    HTMLout=MakeHTML(H,com,errormessage=errormessage)
    buffITOL, buffHIST,H, bins=ForITOL(H)
    #print H["MIByBranch"].loc[:,("I(Ti,G)","Color")]
    secondaryOutput(H,db,com)
    #print H["MIByBranch"].loc[:,("I(Ti,G)","Color")]
    makeITOL3output(db.tree,buffITOL, buffHIST, bins, com)
    #print H["MIByBranch"].loc[:,("I(Ti,G)","Color")]
    
    #CallITOL for graph
    #makeITOLcall(db.tree,buffITOL, buffHIST, com)
    if com["-G"]=="1":
        if com["-H"]=="0":
            if com["-M"]=="TNC":
                makePhyloHOutput(path="./", Z="maximumDepthInMeters", GeoJson=False,prefix=com["-o"],shape="subBiome/subBiome.shp", Sample=com["sample"], groupBy=com["groupBy"])
            else:
                makePhyloHOutput(path="./", Z="maximumDepthInMeters", GeoJson=False,prefix=com["-o"], Sample=com["sample"], groupBy=com["groupBy"])
        else:
            makePhyloHOutput(path="./", Z="maximumDepthInMeters", GeoJson=False,prefix=com["-o"],shape="subgrid_ease/subgrid.shp", Sample=com["sample"], groupBy=com["groupBy"])
    print (os.listdir("./"))
