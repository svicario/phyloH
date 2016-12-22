from pandas import DataFrame, Series,to_datetime, read_csv
import shapefile
from geojson import Feature, Point, FeatureCollection, GeometryCollection
import geojson
from ete2 import PhyloTree
import os
import shutil

def makePhyloHInput(csv="allPMK1970-1996II.tab",sep=",", specialCleaning=False, index_col=None,groupBy="locationID", organismQuantity="individualCount", Z="maximumDepthInMeters", shape=False, makePhylo=False, makeTaxo=False):
    """
    Extract from a csv/tab file a shape file and produce sample and group file for PhyloH input
    The input file need to have a row for each group of organism of the same species observed in a given event.
    Convention for column names are taken from http://rs.tdwg.org/dwc/terms.
    The index_col need to give the number of the column that has the eventID information (the id that identify the event of observation)
    In order to build the Shape file user need to include columns with name decimalLongitude and decimalLatitude.
    The optional time information should be stored in the column name eventDate. pandas tools would enforce shape file convention on that columns
    The output format of eventDate is 8 characters with the first 4 indicating the year, the second two the month, and last two the day
    #Two mode of use: Planned and Random Survey.
    1)In Planned Survey each sample is composed by the observations that have same locationID and same eventID. Grouping is by default locationId but could be also something else
    2)In Random Survey sample and group are all identical to the cellid variable.
    If the variable cellid is present the random survey mode is on, if not it is assumed the presence of locationID and eventID
    """
    #Definitions
    #Z coordinate
    Zs=["minimumElevationInMeters","maximumElevationInMeters", "minimumDepthInMeters","maximumDepthInMeters","minimumDistanceAboveSurfaceInMeters","maximumDistanceAboveSurfaceInMeters"]
    #print sep, "ssssciao"
    D=DataFrame.from_csv(csv, sep=sep,index_col=index_col)
    #Lower Taxon name
    taxon=["scientificName","nameComplete"]
    try:
      Taxon=[x for x in taxon if x in D.columns][0]
    except IndexError:
      raise BaseException("no species columns")
    #Abbundances
    abbundance=["individualCount","organismQuantity"]
    D[Taxon]=D[Taxon].str.replace(" ", "_")
    if organismQuantity==None:
      for a in abbundance:
        if a in D.columns:
          organismQuantity=a
          if organismQuantity==None:
            organismQuantity="individualCount"
            D["individualCount"]=1
    #Defining Grouping
    Location="locationID"
    if "cellid" in D.columns:
        Location="cellid"
        #if "locationID" in D.columns:
            #D["OLDlocationID"]=D.locationID
        #D["locationID"]=D.cellid
    elif "locationID" not in D.columns:
        raise ValueError("No Sampling LocationID Specified!")
    D[Location]=D[Location].astype("str")
    if "eventDate" in D.columns:
        Time=Series([x.isoformat().split("T")[0].replace("-","") for x in to_datetime(D.eventDate,errors="coerce")], index=D.index)
        D.loc[:,"eventDate"]=Time
        #print D.loc[:,"eventDate"]
    Sample="eventID"
    if "eventID" not in D.columns:
                if (Location=="cellid")&(Location!=groupBy):
                    Sample=Location
                elif "eventDate" in D.columns:
                                D.loc[:,"eventID"]=D[Location]+D.loc[:,"eventDate"]
                else:
                    Sample=groupBy
    if specialCleaning:
        D[Taxon]=Series([str(x).replace(" ","_") for x in list(D[Taxon].values)], index=D.index)
        D["maximumDepthInMeters"]=-1*D["maximumDepthInMeters"]
        D["minimumDepthInMeters"]=-1*D["minimumDepthInMeters"]
        D["Deep"]=(D.maximumDepthInMeters< -80)
        D["Deep"].replace({True:"Shallow",False:"Deep"})
        groupBy.append("Deep")
        D=D[D.Higher_taxon=="Echinodermata"]
        DD=D.dropna(subset=["decimalLongitude","decimalLatitude","maximumDepthInMeters","minimumDepthInMeters",organismQuantity,Taxon])
        DD[[organismQuantity]]=DD[[organismQuantity]].astype(int)
        DD["eventID"]=DD.locationID+"-"+DD.eventDate+"-"+Series(map(str,map(int,DD.NumberSample)), index=DD.index)
    else:
        DD=D
    DDD=DD.groupby(Sample)
    proGeoColumns=["decimalLongitude","decimalLatitude"]
    if Z:
        proGeoColumns.append(Z)
    
    namefile=".".join(csv.split("/")[-1].split(".")[:-1])
    #print namefile
    if shape:
        ProSHAPE=DDD[proGeoColumns].first()
        w = shapefile.Writer(shapefile.POINT)
        w.field('locationID')
        #w.field("Date","D",8)
        for station,ProGeo in zip(ProSHAPE.index, ProSHAPE.values):
            #if the ProGeo is of two elements only Lat, Long, while if 3 then Lat, Long, Z
            w.point(*ProGeo)
            w.record(station)
        w.save(namefile)
    
    if makePhylo or makeTaxo:
        #taxo
        Classification=[ "kingdom","phylum","class","order","family","genus","subgenus"]+[Taxon]
        DT=DD[[x for x in Classification if x in DD.columns]]
        ActualClassification=list(DT.columns)
        DT["Accession"]=DT[Taxon]
        DT=DT[["Accession"]+ActualClassification]
        #DT=DT.drop_duplicates()
        if makeTaxo:
            DT.to_csv(namefile+".taxonomy",na_rep="NA",sep="|", header=False, index=False)
        if makePhylo:
            TT=Taxo()
            TT.build(namefile+".taxonomy", accession=True)
            #print makePhyloTaxo, "ciao"
            tree=TT.tree.write(namefile+".tree", format=3)
            #print tree
            handle=open(namefile+".tree","w")
            handle.write(tree)
            handle.close()
            #print namefile
    SampleData=DD[[Sample,organismQuantity, Taxon]]
    SampleData.to_csv("Sample", sep="\t", header=False, index=False)
    
    Group=DDD[[groupBy]].first()
    Group.to_csv("Group_"+groupBy, sep="\t",header=False)
    DD.to_csv(namefile+"_Clean.tab",sep="\t")
    return Sample, groupBy

def makePhyloHOutput(path="./", Z="maximumDepthInMeters", GeoJson=True, prefix="ciccio", shape=None, Sample=None, groupBy=None):
    """
    take the _Clean.tab file, the xml and optinal shapefile file to build a ShapeFile and GeoJson file 
    """
    #Working assumptions
    ##In general I assume that long and lat are in decimal and their name match standard one of TWDG
    ##if no shapefile is given I assume that projection is always WGS 1984 and the csv has always a locationID to group observation  and eventID to identify single sampling event
    ##if a shapefile is given I assume that grouping is done with cellid column assuming a single sample for each group, and do not consider locationID and eventID
    ##If a date it assume eventDate
    ##Z coordinate is built using the different component described in TWDG, see Elevation function for details
    
    #Definitions
    ##Allowed Z components
    Zs=["minimumElevationInMeters","maximumElevationInMeters", "minimumDepthInMeters","maximumDepthInMeters","minimumDistanceAboveSurfaceInMeters","maximumDistanceAboveSurfaceInMeters"]
    ## Default projection
    prj='GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]'
    ## Default coordinates
    proGeoColumns=["decimalLongitude","decimalLatitude"]
    ##Default attributes
    proValueColumns=[Sample,groupBy,"Diversity","Beta_Diversity", "Beta_Pvalue", "SGN_SeqBonferroni"]
    if Sample==groupBy:
        proValueColumns=[Sample,"Diversity","Beta_Diversity", "Beta_Pvalue", "SGN_SeqBonferroni"]
    
    #Getting original observation
    filenametab= [ x for x in os.listdir(path) if x.split("_")[-1]=="Clean.tab"][0]
    D=DataFrame.from_csv(path+"/"+filenametab,sep="\t")
    
    #Searching Z components and Building a single z value
    proZ=[x for x in D.columns if x in Zs]
    if proZ:
      DictZ=Elevation(D[proZ])
      if DictZ.has_key("PointEstimateZ"):
        Z=DictZ["PointEstimateZ"]
      else:
        Z=DictZ["MaxZ"]
      D["Z"]=Z
      #proGeoColumns.append("Z")
      #proValueColumns+=proZ
    
    #filename=[ x for x in os.listdir(path) if x.find("Group")==0][0]
    #z=DataFrame.from_csv(path+filename,sep="\t", header=None, index_col=0)
    #Groups=z.to_dict()[1]
    #print Groups
    if "eventDate" in D.columns:
        #Formatting should be already done in Input
        #Time=Series([x.isoformat().split("T")[0].replace("-","") for x in to_datetime(D.eventDate)], index=D.index)
        #D["eventDate"]=Time
        #proGeoColumns.append("eventDate")
        #proValueColumns=["eventDate"]+proValueColumns
        pass
    #print D.columns
    
    #D["Grouping"]=D["locationID"].replace(Groups)
    D[Sample]=D[Sample].astype("str")
    D[groupBy]=D[groupBy].astype("str")
    DDD=D.groupby(Sample).first()
    #This is necessery because otherwise group columns disapear from data set
    DDD.reset_index(inplace=True)
    print "grouping"
    print Sample
    print DDD.head()
    filename=[ x for x in os.listdir(path) if x.find(prefix+"_MI_KL.csv")>-1][0]
    z=DataFrame.from_csv(path+filename,sep="\t", header=0)
    z.index=z.index.astype("str")
    z.columns=["Beta_Diversity","Beta_Pvalue", "SGN_SeqBonferroni"]
    print( z.index,DDD.columns)
    DDD=DDD.join(z,on=groupBy)
    

    filename=[ x for x in os.listdir(path) if x.find(prefix+"_Gammas.csv")>-1][0]
    z=DataFrame.from_csv(path+filename,sep="\t", header=0)
    z.index=z.index.astype("str")
    DDD=DDD.join(z.Diversity, on=groupBy)
    
    #DDD.reset_index(inplace=True)
    #DDD["locationID"]=DDD["locationID"].astype("str")
    #print DDD
    #print DDD[["Diversity","Beta_Diversity" ]].first()
    print proValueColumns
    ProValues=DDD[proValueColumns].round(4)
    #print ProValues.iloc[0,:]
    ProGeos=DDD[proGeoColumns]
    ProValues.set_index(Sample, inplace=True, drop=False)
    print "perovalues"
    print ProValues
    namefile=filenametab.split("_Clean.tab")[0]
    if not shape:
        w = shapefile.Writer(shapefile.POINTZ)
        #if "eventDate" in D.columns:
            #w.field('eventDate',"D",8)
        if Sample!=groupBy:
            w.field("EventID","C",40)
        w.field('Grouping',"C",40)
        w.field('Diversity',"N",12,7)
        w.field('Beta_Diversity',"N",12,7)
        w.field('Beta_P_value',"N",12,7)
        w.field('Beta_SignSeqBonferroni',"C",5)
        if "Z" in D.columns:
            w.field("Z","N",8,0)
        for ProGeo, ProValue in zip(ProGeos.values,ProValues.values):
            #if the ProGeo is of two elements only Long Lat, while if 3 then Long, Lat, Z, if 4 also eventID in format 8 numbers
            #print ProGeo
            #print [len(str(x)) for x in ProValue]
            w.point(*ProGeo)
            w.record(*ProValue)
            #print w.records
        
        w.save(namefile)
        #I assume that projection is always what written in prj 
        handle=open(namefile+".prj","w")
        handle.write(prj)
        handle.close()
    else:
        #return ProValues
        GridDecorator(shape, ProValues)
    if GeoJson:
        geoPoints=[]
        for i in range(len(ProGeos)):
            temp={}
            temp.update(zip(proValueColumns,list(ProValues.iloc[i,:].values)))
            posvalue=ProGeos.shape[1]
            if ("eventDate" in D.columns) & posvalue>3:
              posvalue-=1
            Position=list(ProGeos.iloc[i,:posvalue].values)
            #print Position
            #print temp
            ID=ProValues.iloc[i]["eventID"]
            geoPoints.append(Feature(geometry=Point(Position), id=ID, properties=temp))
        outjson=FeatureCollection(geoPoints)
        handle=open(prefix+".geojson","w")
        handle.write(geojson.dumps(outjson, indent=4, separators=(',', ': ')))
        handle.close()
    if "eventDate" in ProValues.columns:
        ProValues.drop(["eventDate"], inplace=True, axis=1)
    if "Z" in ProValues.columns:
        ProValues.drop(["Z"], inplace=True, axis=1)
    #print (ProValues.columns,ProGeos.columns)
    OUTD=ProGeos.join(ProValues)
    ZZ=[x for x in D.columns if (x in Zs)&(x not in OUTD.columns)]
    OUTD=OUTD.join(DDD[ZZ])
    OUTD.to_csv(namefile+"_ByeventID.tab",sep="\t", index=False)
    return OUTD

def Elevation(DictZ):
  #Zs=["minimumElevationInMeters","maximumElevationInMeters", "minimumDepthInMeters","maximumDepthInMeters","minimumDistanceAboveSurfaceInMeters","maximumDistanceAboveSurfaceInMeters"]
  #Mins={"ElevationInMeters":0,"DepthInMeters":0,"DistanceAboveSurfaceInMeters":0}
  default={"ElevationInMeters":0,"DepthInMeters":0,"DistanceAboveSurfaceInMeters":0}
  Mins={}
  Maxs={}
  for i in DictZ:
    dif=i.find("minimum")
    if dif>-1:
      Mins[i[7:]]=DictZ[i] 
    else:
      Maxs[i[7:]]=DictZ[i]
  lackMins=set(Maxs.keys()).difference(Mins.keys())
  for l in lackMins:
    Mins[l]=Maxs[l]
  lackMaxs=set(Mins.keys()).difference(Maxs.keys())
  for l in lackMaxs:
    Maxs[l]=Mins[l]
  lackBoth=set(default.keys()).difference(Mins.keys())
  for k in lackBoth:
    Maxs[k]=default[k]
    Mins[k]=default[k]
  Max=CalcElevation(Maxs)
  Min=CalcElevation(Mins)
  if all(Max==Min):
    return {"PointEstimateZ":Max}
  else:
    return {"MaxZ":Max, "MinZ":Min}
  

def CalcElevation(DictZ):
  return DictZ["ElevationInMeters"]-DictZ["DepthInMeters"]+DictZ["DistanceAboveSurfaceInMeters"]

def GridMaker(samplefilename, sizeg, pathPythonScript):
    from GrassLand import *
    dfOr=DataFrame.from_csv(samplefilename,sep=",")
    #dfOr.to_csv("pipe_"+samplefilename,sep="|")
    gisbase,gisdb, mapset=SetUp(pathPythonScript=pathPythonScript)
    LoadData(gisbase,gisdb, mapset, samplefile=samplefilename)
    hexagrid(gisbase,gisdb, mapset, radius=sizeg)
    #result.csv need to changed to match input defined in com[-s]
    #rowname in sample.csv become first column, so need to fix again them as rownames, so index_col=1
    df=DataFrame.from_csv("result.csv",sep="|", index_col=1)
    #dfOr=DataFrame.from_csv(samplefilename,sep="|"
    dfOr.loc[:,"cellid"]="N"+df.cellid.astype("str")
    dfOr.to_csv(os.path.dirname(os.path.abspath(samplefilename))+os.path.sep+"grid_"+os.path.basename(samplefilename),sep=",")

def ShapeLoader(shapeName,samplefilename ,pathPythonScript, group="WWF_MHTNAM"):
    from GrassLand import *
    dfOr=read_csv(samplefilename,sep=",")
    gisbase,gisdb, mapset=SetUp(pathPythonScript=pathPythonScript)
    LoadData(gisbase,gisdb, mapset, samplefile=samplefilename)
    addVec(gisbase, gisdb,mapset,pathPythonScript=pathPythonScript,dirName=os.path.dirname(shapeName), fileName=os.path.basename(shapeName))
    Cross(gisbase,gisdb, mapset, group=group,location="ll", vector_name="Biome", point_name="sample")
    df=DataFrame.from_csv("result.csv",sep="|", index_col=None)
    dfOr.loc[:,"cellid"]=df.cellid
    prefix=""
    if df.loc[:,group].dtype!=Series(["ww"]).dtype:
        prefix="N"
    df.loc[:,group]=[x.replace(" ", "_") for x in df[group].astype("str").tolist()]
    print dfOr.shape, df.shape
    dfOr.loc[:,group]=prefix+df.loc[:,group]
    #dfOr.loc[:,"cellid"]=df.cellid
    print(dfOr.head())
    dfOr=dfOr[dfOr.cellid.notnull()]
    dfOr.loc[:,"cellid"]="N"+dfOr.cellid.astype("int").astype("str")
    dfOr.to_csv(os.path.dirname(os.path.abspath(samplefilename))+os.path.sep+"grid_"+os.path.basename(samplefilename),sep=",")
#def GridMaker(samplefilename, sizeg):
    #"""
    #SetUp a grid in projection Equal Area Ease 2.0 and add a column in observation CSV to indicates to what cell belong each observation
    #"""
    #import os
    ##load samplefile and format ready for grass with newfilename sample.csv
    #sf=DataFrame.from_csv(samplefilename, sep=",")
    #sf.to_csv("sample.csv", index=False, sep="|")
    #os.system(".\load.py")
    #os.system(".\hgrid.py "+str(sizeg))
    ##result.csv need to changed to match input defined in com[-s]
    ##rowname in sample.csv become first column, so need to fix again them as rownames, so index_col=1
    #df=DataFrame.from_csv("result.csv",sep="|", index_col=1)
    #dfOr=DataFrame.from_csv("sample.csv",sep="|")
    #dfOr.loc[:,"cellid"]=df.cellid
    #dfOr.to_csv("sampleRich.csv",sep=",")

def GridDecorator(filename, dictRES, addN=True):
    """
    Open the grid previously defined and write a new shape file in which this grid is decorated with attribute of diversity
    The link variable between filename and dictRES is the variable Grouping 
    """
    print("GRIDDECORATOR")
    print dictRES.head()
    from pandas import isnull
    pathname=os.path.dirname(filename)
    basename=os.path.basename(filename)
    R=shapefile.Reader(filename)
    w = shapefile.Writer(R.shapeType)
    if dictRES.shape[1]==6:
        w.field('EventID',"C",40)
    w.field('Grouping',"C",40)
    w.field('Diversity',"N",12,7)
    w.field('Beta_Diversity',"N",12,7)
    w.field('Beta_P_value',"N",12,7)
    w.field('Beta_SignSeqBonferroni',"C",6)
    print dictRES.head()
    for shapeRecord in R.shapeRecords():
        Key=str(shapeRecord.record[0])
        if addN:
            Key="N"+Key
        if Key in dictRES.index:
            if not any(isnull(dictRES.loc[Key,:].iloc[-2:])):
                if len(shapeRecord.shape.parts)==1:
                    w.poly(parts=[shapeRecord.shape.points])
                else:
                    Start=[None]+list(shapeRecord.shape.parts)
                    End=list(shapeRecord.shape.parts)+[None]
                    Pair=zip(Start,End)[1:]
                    WKT=[shapeRecord.shape.points[x[0]:x[1]] for x in Pair]
                    w.poly(parts=WKT)
                #print [shapeRecord.record[0]]+list(dictRES.loc[str(shapeRecord.record[0]),[groupBy, "Diversity", "Beta_Diversity"]])
                w.record(*list(dictRES.loc[Key,:]))
            else:
                print(dictRES.loc[Key,:])
        else:
            print(Key)
    w.save(os.path.join(pathname,"decorated_"+basename))
    shutil.copy(os.path.join(pathname,"".join(basename.split(".")[:-1])+".prj"),os.path.join(pathname,"".join(("decorated_"+basename).split(".")[:-1])+".prj"))
    
class Taxo:
    """
    Handle taxonomic information and sequence redundancy and possible conflict with taxonomy.
    Taxonomy is stored as a phylogenetic tree with branch length equal to 1.
    """
    def __init__(self):
        self.taxoDB={}
        self.tree=PhyloTree()
        self.tree.name="NoName"
    def CleanTaxonName(self,List):
        self.CleanedTaxon={}
        for s in List:
            #s=s.replace("_"," ")
            sout=s.strip().split(" sp")[0]
            sout=sout.split(" strain")[0]
            sout=sout.split("_strain")[0]
            self.CleanedTaxon[s]=sout
    def MatchTaxo2Seq(self,filenametax,filenameseq):
        pass
    def growingTree(self,node, taxoPath):
        """
        Using a taxonomic path, list of taxonomic name ordered from more general to more specific,
        build a tree with name matching this taxonomy.
        It assume that base node is the correct one
        """
        if not taxoPath:
            return "END"
        NEW=True
        for n in node.children:
            if n.name==taxoPath[0]:
                try:
                    taxoPath.pop(0)
                except IndexError:
                    print n.name 
                self.growingTree(n, taxoPath)
                NEW=False
                break
        if NEW:
            name=taxoPath.pop(0)
            n=node.add_child(name=name)
            self.growingTree(n, taxoPath)
        
    def build(self,filenametax, accession=False):
        """
        process a taxonomy file and build a taxonomic tree
        """
        self.filenametax=filenametax
        handle=open(filenametax,"r")
        TAX=handle.readlines()
        handle.close()
        self.TAX=[]
        TAX=[x.split("|") for x in TAX]
        for tax in TAX:
            #reformat species name
            #I assume that last name of the path is the species name
            #reformatting are not fullproof.
            if len(tax)<2: print tax
            speciesName=tax[-1].strip()
            #speciesName=speciesName.strip().replace("sp.","sp").replace("_"," ").replace("-"," ")
            tax[-1]=speciesName
            if accession:
                Naccession=tax.pop(0)
                self.taxoDB[Naccession]=speciesName
            basename=tax.pop(0)
            if self.tree.name=="NoName":
                self.tree.name=basename
            assert self.tree.name==basename, (tax,basename,self.tree.name)
            self.TAX.append([basename]+tax)
            self.growingTree(self.tree, tax)
        self.TAX=set(["|".join(x) for x in self.TAX])
        self.TAX=[x.split("|") for x in self.TAX]
