import os
import sys
import subprocess
import tarfile
import shutil
from pandas import read_csv

Locations="AncillaryFiles/locations.tgz"
addon="AncillaryFiles"
def SetUp(pathPythonScript):
    """
    Change Enviromental Variables to work with GRASS and setup working DB
    """
    #Set up locations
    wdIN=os.path.dirname(os.path.abspath(pathPythonScript))
    wdout=os.getcwd()
    
    os.makedirs(os.path.join(wdout,"grass/geodb7"))
    
    grass7bin = 'grass74'
    startcmd = [grass7bin, "-c", "EPSG:4326", "-e", "grass/geodb7/ll"]
    
    p = subprocess.Popen(startcmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
        print >>sys.stderr, "ERROR: Cannot find GRASS GIS 7 start script (%s)" % startcmd
        print err
        sys.exit(-1)
        
    startcmd = [grass7bin, "-c", "EPSG:6933", "-e", "grass/geodb7/cea"]
    
    p = subprocess.Popen(startcmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
        print >>sys.stderr, "ERROR: Cannot find GRASS GIS 7 start script (%s)" % startcmd
        sys.exit(-1)
    
    
    grass7bin = 'grass74'
    #gisdb = os.path.join(os.path.expanduser("~"), "grass", "geodb7")
    gisdb=os.path.join(os.path.realpath(wdout), "grass", "geodb7")
    mapset   = "PERMANENT"
    
    startcmd = [grass7bin, '--config', 'path']
    
    p = subprocess.Popen(startcmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
        print >>sys.stderr, "ERROR: Cannot find GRASS GIS 7 start script (%s)" % startcmd
        sys.exit(-1)
    gisbase = out.strip('\n\r')
    
    os.environ['GISBASE'] = gisbase
    # the following not needed with trunk
    os.environ['PATH'] += os.pathsep + os.path.join(gisbase, 'extrabin')
    # add path to GRASS addons
    home = os.path.expanduser("~")
    os.environ['PATH'] += os.pathsep + os.path.join(home, '.grass7', 'addons', 'scripts')
    
    # define GRASS-Python environment
    gpydir = os.path.join(gisbase, "etc", "python")
    sys.path.append(gpydir)
    
    ########### DATA
    # Set GISDBASE environment variable
    os.environ['GISDBASE'] = gisdb
    return gisbase,gisdb, mapset

def LoadData(gisbase,gisdb, mapset, location="ll", samplefile="sample.csv", sep=","):
    # find decimalLatitude and Longitude in table
    DF=read_csv(samplefile, header=0,sep=sep)
    print DF.columns
    Lat=[number for number,value in enumerate(DF.columns) if value=="decimalLatitude"][0]+1
    Long=[number for number,value in enumerate(DF.columns) if value=="decimalLongitude"][0]+1
    def MakeUp(label,size): return " ".join(["'"+label.upper()+"'","varchar({size})".format(size=size)])
    form=[ MakeUp(x, DF.loc[:,x].astype("str").apply(len).max()) for x in DF.columns]
    form[Lat-1]="decimalLatitude double precision"
    form[Long-1]="decimalLongitude double precision"
    # import GRASS Python bindings (see also pygrass)
    import grass
    import grass.script as gscript
    import grass.script.setup as gsetup
    ###########
    # launch session
    gsetup.init(gisbase,gisdb, location, mapset)
    gscript.run_command('v.in.ascii',input=samplefile,output='sample',x=Long,y=Lat,skip=1,overwrite=True,quiet=True, separator=sep)

def MakeHexagrid(gisbase,gisdb, mapset, location="cea",radius=15000, point_name="sample", point_loc="ll"):
    import grass
    reload(grass)
    import grass.script as gscript
    import grass.script.setup as gsetup
    gsetup.init(gisbase,
            gisdb, location, mapset)
    radius=int(radius)
    if radius==0:
        radius = 10000
    radius2 = radius*2
    gscript.run_command('v.proj',location=point_loc,input=point_name,output='sample',overwrite=True,quiet=True)
    gscript.run_command('g.region',vect='sample',quiet=True)
    gscript.run_command('g.region',n='n+%d' % radius2 ,e='e+%d' % radius2, \
            w='w-%d' % radius2,s='s-%d' % radius2,quiet=True)
    gscript.run_command('v.mkgrid',flags='h',map='grid',box=(radius,radius),overwrite=True,quiet=True)

def addVec(gisbase,gisdb, mapset, pathPythonScript, dirName='terr-ecoregions-TNC/' ,fileName='tnc_terr_ecoregions',location="ll", name="Biome"):
    import grass
    reload(grass)
    import grass.script as gscript
    import grass.script.setup as gsetup
    gsetup.init(gisbase,gisdb, location, mapset)
    wdIN=os.path.dirname(os.path.abspath(pathPythonScript))
    #temp=os.path.join(wdIN,addon,dirName)
    gscript.run_command('v.in.ogr',input=dirName,layer=fileName,output=name, overwrite=True)

def Cross(gisbase,gisdb, mapset, location="cea", vector_name="grid", point_name="sample", group=None):
    import grass
    reload(grass)
    import grass.script as gscript
    import grass.script.setup as gsetup
    gsetup.init(gisbase,
            gisdb, location, mapset)
    suffix=""
    if location=="cea":
        suffix="_ease"
    gscript.run_command('v.out.ogr',input=point_name,output=point_name+suffix,overwrite=True,quiet=True) 
    gscript.run_command('v.out.ogr',input=vector_name,output=vector_name+suffix,overwrite=True,quiet=True,flags='m')
    gscript.run_command('v.db.addcolumn',map='sample',column='cellid integer',quiet=True)
    gscript.run_command('v.db.addcolumn',map='sample',column=group+ " VARCHAR(10)",quiet=True)
    gscript.run_command('v.what.vect',map='sample',column='cellid',query_map=vector_name,query_column='cat',quiet=True)
    if group:
        gscript.run_command('v.what.vect',map='sample',column=group,query_map=vector_name,query_column=group,quiet=True)
    result = gscript.read_command('db.select',sql='select * from '+point_name) 
    with open('result.csv','w') as t:
        t.write(result)
        t.close()
    s = gscript.read_command('db.select', flags='c', \
            sql='select cellid from sample where cellid is not null')
    a = set(s.split('\n'))
    a.discard('')
    b = str(','.join(str(e) for e in a))
    subname='sub'+vector_name
    gscript.run_command('v.extract',input=vector_name,output=subname, where="cat in (%s)" % b,overwrite=True,quiet=True)
    gscript.run_command('v.out.ogr',input=subname,output=subname+suffix,overwrite=True,quiet=True,flags="m",format="ESRI_Shapefile")


def hexagrid(gisbase,gisdb, mapset, location="cea",radius=15000):
    import grass
    reload(grass)
    import grass.script as gscript
    import grass.script.setup as gsetup
    gsetup.init(gisbase,
            gisdb, location, mapset)
    radius=int(radius)
    if radius==0:
        radius = 10000
    radius2 = radius*2
    gscript.run_command('v.proj',location='ll',input='sample',output='sample',overwrite=True,quiet=True)
    gscript.run_command('g.region',vect='sample',quiet=True)
    gscript.run_command('g.region',n='n+%d' % radius2 ,e='e+%d' % radius2, \
            w='w-%d' % radius2,s='s-%d' % radius2,quiet=True)
    gscript.run_command('v.mkgrid',flags='h',map='grid',box=(radius,radius),overwrite=True,quiet=True)
    gscript.run_command('v.out.ogr',input='sample',output='sample_ease',overwrite=True,quiet=True) 
    gscript.run_command('v.out.ogr',input='grid',output='grid_ease',overwrite=True,quiet=True)
    gscript.run_command('v.db.addcolumn',map='sample',column='cellid integer',quiet=True)
    gscript.run_command('v.what.vect',map='sample',column='cellid',query_map='grid',query_column='cat',quiet=True)
    result = gscript.read_command('db.select',sql='select * from sample') 
    with open('result.csv','w') as t:
        t.write(result)
        t.close()
    s = gscript.read_command('db.select', flags='c', \
            sql='select cellid from sample where cellid is not null')
    a = set(s.split('\n'))
    a.discard('')
    b = str(','.join(str(e) for e in a))
    gscript.run_command('v.extract',input='grid',output='subgrid', where="cat in (%s)" % b,overwrite=True,quiet=True)
    gscript.run_command('v.out.ogr',input='subgrid',output='subgrid_ease',overwrite=True,quiet=True,format="ESRI_Shapefile")
