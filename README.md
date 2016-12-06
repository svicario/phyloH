phyloH
======

The script estimate the phylogenetic entropy and the phylogenetic mutual information across sample given a reference tree and a table of observations

The script user the API https://github.com/albertyw/itol-api to load the reference 
tree and the estimate of entropy and mutual information on itol web site (http://itol.embl.de/) and produce an output.

With Complements to: [ETE2 python library](http://pythonhosted.org/ete2/),[iTOL (Interactive Tree of Life)](http://itol.embl.de/), [Python API for the Interactive Tree of Life ](https://github.com/albertyw/itol-api), [numpy](http://www.numpy.org/),[scipy](http://www.scipy.org/), [Biopython](http://biopython.org/wiki/Main_Page)
###RUNNING FROM COMMAND LINE
python esecutorePhyloHPandas.py [-option] [value]

option -s, -r  are mandatory  
if -G 0 then also -f -g are mandatory  
-H is read only if -G 1  
-M is read only if -G 1 and -H different from zero  

to obtain the help call
"python esecutorePhyloHPandas.py"


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
 
  
 run example using following call 
 ../esecutorePhyloHPandas.py -f test/Echinodermata.tree -s test/sampleTest  -g test/GroupTest -q 1 -r 2 -o xxx -e 0 -k 0
 
 The script is exposed behind a webservice. For details visit https://www.biodiversitycatalogue.org/rest_methods/143

Send e-mails to Saverio Vicario (the repo owner) for questions and comments