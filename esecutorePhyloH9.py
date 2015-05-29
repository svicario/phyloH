from ete2 import Tree, nexml
from lib_script.entropyPerm5 import *
from lib_script.DiversityTreeKLFabris import *
from lib_script.graphic import *
from StringIO import StringIO
import copy


import sys


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

def mapWidthANDColor2Phyloxml(node):
    color=phyloxml.BranchColor(**fromhex(node.img_style["hz_line_color"]))
    node.phyloxml_clade.set_color(color)
    node.phyloxml_clade.set_width(node.img_style["hz_line_width"])

def collapse(node, recipient=None,prefix="QUERY_", noSingleton=False):
    """
    Collapse branch with only leaf with   "Query" prefix on name
    """
    for i in node.get_children():
        if i.is_leaf():
            if i.name.find(prefix)==0:
                #recipient.children.append(queryOne)
                if noSingleton:
                    pass
        if all([x.find(prefix)==0 for x in i.get_leaf_names()]):
            newi=Tree()
            newi.name=i.name
            newi.dist=i.dist
            if recipient:
                    recipient.children.append(newi)
            for c in i.get_children():
                i.remove_child(c)
                newi.children.append(c)
        else:
                collapse(i,recipient)

def tabulatePrunedBranches(recipient):
    """
    make a tabbed file of kept ancestor and removed children
    """
    header="Kept\tCollapsed"
    Out=header
    for child in recipient.children:
        Out+="\n"+child.name+"\t"+",".join(child.get_leaf_names())
    if Out==header:
        Out=""
    return Out
def QRdata(db):
    """
    function that modify sample an and group so to observe difference
    between observed sample ( with name that start with "QUERY") and
    reference one used only to build the tree
    Reference sequence 
    """
    totNames=db.tree.get_leaf_names()
    L=len(totNames)
    Qlist=[y for y in totNames if "QUERY" in y]
    if Qlist == []: raise InputError, "no Query"
    Rlist=[y for y in totNames if "QUERY" not in y]
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

def QRtree(db, R):
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
    for n in db.tree.traverse():
        if n.a[db.groups["Q"]]<n.a[db.groups["R"]]:
            n.betaA=0
        else:
            betaAplus+=n.betaA
    R['MI_treeAndEnvironment']= betaAplus
    R['tag']= array([ 0,  0])
    R['MI_treeAndSampleGivenEnvironment']= 0.0
    R['Hgamma']= 0
    R['HalphaByEnvironment'] = 0
    R['MI_KL']= numpy.array([ 0,  0])
    R['HgammaEachEnvironment'] = numpy.array([ 0,  0])
    R['tot'] = 0
    R['HalphaBySamples']= 0
    
if __name__=="__main__":
    com={"-q":1, "-qt":None, "-x":"nexml","-R":"/opt/exp_soft/uniba/R/bin/","None":1, "--QR":"0","--QRC":"1"}
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
     -f filename       Use this file as the phylogeny file [phylo].
     -o filename       Use this file to record output
     -t filename       Use this file to get taxonomic information on some tips and map them in internal node
     -s filename       Use this file as the sample file [sample].
     -g filename       Use this file as the group file [group]
     -r INT            Number of randomizations to use [999]
     -q float          q parameter in the hill series (q=1 index is beta is Chao phylogenetic entropy, q=2 Rao phylogenetic diversity,q  zero is faith phylogenetic diversity)
     -x string         two possible strings :"nexml" or "phyloxml" to select the xml output of the results
     -h 0 or 1         boolean to check if you want html output
     --QR  0 or 1             identify linneage present in the Query but not in the Reference (Need found observation to be tagged with "Query" prefix)
     --QRC 0 or 1            collapse branch with only query before analysis
    """
    if (( (com["--QR"]=="0") and (not '-s' in com) ) | (not '-f' in com )):
        print spiegazione
        raise ImportError
    call=" ".join(sys.argv)
    #Check if pyQT4 is there
    try:
        from ete2 import NodeStyle, TreeStyle, faces
    except ImportError:
        del com["-qt"]
    #Loading the data
    db=DBdata()
    db.readTree(com['-f'])
    if com["--QR"]=="0":
        db.readSampleTable(com["-s"])
        if com.has_key("-g"):
            db.readGroupTable(com["-g"])
        else:
            db.noGroup()
    else:
        if com["--QRC"]=="1":
            QueryTreeSegment=Tree()
            collapse(db.tree,QueryTreeSegment)
            tabPrune=tabulatePrunedBranches(QueryTreeSegment)
            if tabPrune:
                handle=open("PreliminaryCol.tab","w")
                handle.write(tabPrune)
                handle.close
        sample=QRdata(db)
        db.readSampleTable(sample)
        db.noGroup()
    
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

    outdict=makeXMLoutput(HS,HE,db,Names,P,countbyS, R, com)
    
    if bool(int(com["-h"])):
        outhtml=makeHMTLoutput(db,countbyS, HS, HE, R, com, outdict,call)
    
    #Define reference taxonomy
    NodeTaxonDB={}
    if "-t" in com:
        handle=open(com["-t"],"r")
        TAX=handle.readlines()
        handle.close()
        TAX=[x.split("|") for x in TAX]
        taxDB={}
        for i in TAX:
            speciesName=i[0]
            # species names in taxonomy have this modification.
            #speciesName=speciesName.strip().replace("sp.","sp").replace("_"," ").replace("-"," ")
            taxDB[speciesName]=i[1:]
        
        if taxDB.has_key(" "):
            del taxDB['']
        #build node translation map
        NodeTaxonDB={}
        bad=[]
        check=[]
        def NodeName(node):
            #taxroot=taxDB[" ".join(l[0].split("_")[1:])]
            l=node.get_leaf_names()
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
            for c in node.get_children():
                if not c.is_leaf():
                    Traverse(c)
                    if NodeTaxonDB[c.name]==NodeTaxonDB[node.name]:
                        NodeTaxonDB[c.name]=c.name
        
        Traverse(db.tree)
        #print taxDB["FJ402946.1.1211_U"]
        #Use node translation map
        #for i in db.tree.traverse():
        #    if not i.is_leaf():
        #        i.name=NodeTaxonDB[i.name]
    #Put the correct proportion on each branches
    R=db.GetEntropies(q=com["-q"])
    if com["--QR"]=="1":
        QRtree(db,R)
    #remove node Assigned


    QueryTreeSegment=Tree()
    collapse(db.tree,QueryTreeSegment)
    tabPrune=tabulatePrunedBranches(QueryTreeSegment)
    if tabPrune:
        handle=open(com["-o"]+"_Collapsed.tab","w")
        handle.write(tabPrune)
        handle.close
    def seebeta(n):
        return n.betaA
    N=list(db.tree.traverse())
    N.sort(key=seebeta, reverse=True)

    if bool(int(com["-h"])):
        makeHMTLTree(N,outhtml, NodeTaxonDB,db, SN, countbyGroup, com,tabPrune)

    [layout_branch(x,"betaA") for x in db.tree.traverse()]
    for i in SN:
        pass
        #i.img_style["bgcolor"]='#FFFF00'

    if com["-x"]=="nexml":
        [Features2NeXML(x,features=["pvalue"]) for x in SN]
        [Features2NeXML(x,features=["betaA","betaCA"]) for x in db.tree.traverse()]
        nexml_project.export(open(com["-o"]+".nexml","w"))

    if bool(int(com["-h"])):
        makeITOLcall(NodeTaxonDB,com,SN, db)

    try:
        [layout_branch(x,"betaCA") for x in db.tree.traverse()]
    except ValueError:
        print "Experimental design does not seem to include replicates and  I(T,S|E)should be equal to zero"
    else:
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









