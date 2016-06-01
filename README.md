phyloH
======

The script estimate the phylogenetic entropy and the phylogenetic mutual information across sample given a reference tree and a table of observations

The script user the API https://github.com/albertyw/itol-api to load the reference 
tree and the estimate of entropy and mutual information on itol web site (http://itol.embl.de/) and produce an output.

With Complements to: [ETE2 python library](http://pythonhosted.org/ete2/),[iTOL (Interactive Tree of Life)](http://itol.embl.de/), [Python API for the Interactive Tree of Life ](https://github.com/albertyw/itol-api), [numpy](http://www.numpy.org/),[scipy](http://www.scipy.org/), [Biopython](http://biopython.org/wiki/Main_Page)
###RUNNING FROM COMMAND LINE
python esecutorePhyloH5.py [-option] [value]

option -f  -s, -g, -r are mandatory

to obtain the help call
"python esecutorePhyloH9.py"


     -f filename       Use this file as the phylogeny file [phylo].
     -o filename       Use this file to record output
     -t filename       Use this file to get taxonomic information on some tips and map them in internal node
     -s filename       Use this file as the sample file [sample].
     -g filename       Use this file as the group file [group]
     -r INT            Number of randomizations to use [999]
     -q float          q parameter in the hill series (q=1 index is beta is Chao phylogenetic entropy, q=2 Rao phylogenetic diversity,q  zero is faith phylogenetic diversity)
     -x string         two possible strings :"nexml" or "phyloxml" to select the xml output of the results
     -h 0 or 1         boolean to check if you want html output
     --QR              identify linneage present in the Query but not in the Reference (Need found observation to be tagged with "Query" prefix)
     --QRC             collapse branch with only query before analysis
     -e   
  
 run example using following call (for you Marica)
 ../esecutorePhyloHPandas.py -f test/Echinodermata.tree -s test/sampleTest  -g test/GroupTest -q 1 -r 2 -o xxx -e 0 -k 0
 
 The script is exposed behind a webservice. For details visit https://www.biodiversitycatalogue.org/rest_methods/143

Send e-mails to Saverio Vicario (the repo owner) for questions and comments