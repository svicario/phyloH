from pandas import DataFrame, Series,to_datetime
import shapefile
from geojson import Feature, Point, FeatureCollection, GeometryCollection
import geojson
from ete2 import PhyloTree
import os

def makePhyloHInput(csv="allPMK1970-1996II.tab",sep=",", specialCleaning=False, index_col=None,groupBy=["locationID"], organismQuantity="individualCount", Z="maximumDepthInMeters", shape=False, makePhyloTaxo=False):
    """
    Extract from a csv/tab file a shape file and produce sample and group file for PhyloH input
    The input file need to have a row for each group of organism of the same species observed in a given event.
    Convention for column names are taken from http://rs.tdwg.org/dwc/terms.
    The index_col need to give the number of the column that has the eventID information (the id that identify the event of observation)
    In order to build the Shape file user need to include columns with name decimalLongitude and decimalLatitude.
    The optional time information should be stored in the column name eventDate. pandas tools would enforce shape file convention on that columns
    The output ormat of eventDate is 8 characters with the first 4 indicating the year, the second two the month, and last two the day
    
    """
    Zs=["minimumElevationInMeters","maximumElevationInMeters", "minimumDepthInMeters","maximumDepthInMeters","minimumDistanceAboveSurfaceInMeters","maximumDistanceAboveSurfaceInMeters"]
    #print sep, "ssssciao"
    D=DataFrame.from_csv(csv, sep=sep,index_col=index_col)
    taxon=["scientificName","nameComplete"]
    #print D.columns
    try:
      Taxon=[x for x in taxon if x in D.columns][0]
    except IndexError:
      raise BaseException("no species columns")
    abbundance=["individualCount","organismQuantity"]
    D[Taxon]=D[Taxon].str.replace(" ", "_")
    if organismQuantity==None:
      for a in abbundance:
        if a in D.columns:
          organismQuantity=a
          if organismQuantity==None:
            organismQuantity="individualCount"
            D["individualCount"]=1
    Location="locationID"
    if "locationID" not in D.columns:
    	if "location" in  D.columns:
    		Location="location"
    	else:
    		raise ValueError("No Sampling Location Specified!")
    D[Location]=D[Location].astype("str")
    if "eventDate" in D.columns:
        Time=Series([x.isoformat().split("T")[0].replace("-","") for x in to_datetime(D.eventDate,errors="coerce")], index=D.index)
        D.loc[:,"eventDate"]=Time
        #print D.loc[:,"eventDate"]
    if "eventID" not in D.columns:
                if "eventDate" in D.columns:
                                D.loc[:,"eventID"]=D[Location]+D.loc[:,"eventDate"]
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
    DDD=DD.groupby("eventID")
    proGeoColumns=["decimalLongitude","decimalLatitude"]
    if Z:
        proGeoColumns.append(Z)
    namefile=".".join(csv.split("/")[-1].split(".")[:-1])
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
    
    Classification=[ "kingdom","phylum","class","order","family","genus","subgenus"]+[Taxon]
    DT=DD[[x for x in Classification if x in DD.columns]]
    ActualClassification=list(DT.columns)
    DT["Accession"]=DT[Taxon]
    DT=DT[["Accession"]+ActualClassification]
    #DT=DT.drop_duplicates()
    DT.to_csv(namefile+".taxonomy",na_rep="NA",sep="|", header=False, index=False)
    TT=Taxo()
    TT.build(namefile+".taxonomy", accession=True)
    if makePhyloTaxo:
        #print makePhyloTaxo, "ciao"
        tree=TT.tree.write(namefile+".tree", format=3)
        #print tree
        handle=open(namefile+".tree","w")
        handle.write(tree)
        handle.close()
        #print namefile
    Sample=DD[["eventID",organismQuantity, Taxon]]
    Sample.to_csv("Sample", sep="\t", header=False, index=False)
    for group in groupBy:
        Group=DDD[[group]].first()
        Group.to_csv("Group_"+group, sep="\t",header=False)
    DD.to_csv(namefile+"_Clean.tab",sep="\t")
    return DD

def makePhyloHOutput(path="./", Z="maximumDepthInMeters", GeoJson=True, prefix="ciccio"):
    """
    take the _Clean.tab file and the xml file from JST to build a ShapeFile and GeoJson file 
    """

    #I assume that projection is always WGS 1984 and long and lat are in decimal and the csv has always a locationID and eventID
    #If a date it assume eventDate
    Zs=["minimumElevationInMeters","maximumElevationInMeters", "minimumDepthInMeters","maximumDepthInMeters","minimumDistanceAboveSurfaceInMeters","maximumDistanceAboveSurfaceInMeters"]
    
    prj='GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]'
    proGeoColumns=["decimalLongitude","decimalLatitude"]
    proValueColumns=["locationID","eventID","Grouping","Diversity","Beta_Diversity"]
    filenametab= [ x for x in os.listdir(path) if x.split("_")[-1]=="Clean.tab"][0]
    D=DataFrame.from_csv(path+"/"+filenametab,sep="\t")
    proZ=[x for x in D.columns if x in Zs]
    if proZ:
      DictZ=Elevation(D[proZ])
      if DictZ.has_key("PointEstimateZ"):
        Z=DictZ["PointEstimateZ"]
      else:
        Z=DictZ["MaxZ"]
      D["Z"]=Z
      proGeoColumns.append("Z")
      proValueColumns+=proZ
    
    filename=[ x for x in os.listdir(path) if x.find("Group")==0][0]
    z=DataFrame.from_csv(path+filename,sep="\t", header=None)
    Groups=z.to_dict()[1]
    if "eventDate" in D.columns:
        #Formatting should be already done in Input
        #Time=Series([x.isoformat().split("T")[0].replace("-","") for x in to_datetime(D.eventDate)], index=D.index)
        #D["eventDate"]=Time
        #proGeoColumns.append("eventDate")
        proValueColumns=["eventDate"]+proValueColumns
    #print D.columns
    
    D["Grouping"]=D["locationID"].replace(Groups)
    D["Grouping"]=D["Grouping"].astype("str")
    
    DDD=D.groupby(["Grouping"]).first()
    filename=[ x for x in os.listdir(path) if x.find(prefix+"_MI_KL.csv")>-1][0]
    z=Series.from_csv(path+filename,sep="\t", header=0)
    z.index=z.index.astype("str")
    z.name="Beta_Diversity"
    DDD=DDD.join(z)
    

    filename=[ x for x in os.listdir(path) if x.find(prefix+"_Gammas.csv")>-1][0]
    z=DataFrame.from_csv(path+filename,sep="\t", header=0)
    z.index=z.index.astype("str")
    DDD=DDD.join(z.Diversity)
    
    DDD.reset_index(inplace=True)
    DDD["locationID"]=DDD["locationID"].astype("str")
    #print DDD
    #print DDD[["Diversity","Beta_Diversity" ]].first()
    ProValues=DDD[proValueColumns].round(4)
    #print ProValues.iloc[0,:]
    ProGeos=DDD[proGeoColumns]
    print ProGeos
    w = shapefile.Writer(shapefile.POINTZ)
    if "eventDate" in D.columns:
        w.field('eventDate',"D",8)
    w.field('locationID',"C",40)
    w.field('EventID',"C",40)
    w.field('Grouping',"C",40)
    w.field('Diversity',"N",8,7)
    w.field('Beta_Diversity',"N",8,7)
    if "Z" in D.columns:
        w.field("Z","N",8,0)
    for ProGeo, ProValue in zip(ProGeos.values,ProValues.values):
        #if the ProGeo is of two elements only Long Lat, while if 3 then Long, Lat, Z, if 4 also eventID in format 8 numbers
        #print ProGeo
        #print [len(str(x)) for x in ProValue]
        w.point(*ProGeo)
        w.record(*ProValue)
        #print w.records
    namefile=filenametab.split("_Clean.tab")[0]
    w.save(namefile)
    #I assume that projection is always what written in prj 
    handle=open(namefile+".prj","w")
    handle.write(prj)
    handle.close()
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
        #print outjson
        handle=open(namefile+".geojson","w")
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