phyloH
======

The script estimate the phylogenetic entropy and the phylogenetic mutual information across sample given a reference tree and a table of observations

The script user the API https://github.com/albertyw/itol-api to load the reference 
tree and the estimate of entropy and mutual information on itol web site (http://itol.embl.de/) and produce an output.

HOW to use:
python esecutorePhyloH5.py [-option] [value]
option -f  -s, -g, -r are mandatory
to obtain the help call
python esecutorePhyloH5.py


 -f filename	   Use this file as the phylogeny file [phylo].
 -o filename       Use this file to record output
 -s filename       Use this file as the sample file [sample].
 -g filename       Use this file as the group file [group]
 -r INT		   Number of randomizations to use [999]
 -q float          q parameter in the hill series (q=1 index is beta is Chao phylogenetic entropy, q=2 Rao phylogenetic diversity,q  zero is faith phylogenetic diversity)
 -x string         two possible strings :"nexml" or "phyloxml" to select the xml output of the results
 -h 0 or 1         boolean to check if you want html output   
 
 The script is exposed behind a webservice. For details visit https://www.biodiversitycatalogue.org/rest_methods/143
