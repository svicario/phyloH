tar -xzf terr-ecoregions-TNC.tar.gz

python phyloH/esecutorePhyloHPandas.py -M TNC -G 1  -s SAMPLE_FILENAME -g "WWF_MHTNAM" -r 100 -o out -e 0 -k 0

# collect the output
tar -cvfz Fulloutput.tgz --ignore-failed-read --hard-dereference external/ out* PhyloCommunity.js sample sub* */decorated* SAMPLE_FILENAME*  

tar -cvfz --ignore-failed-read --hard-dereference Shapefile.tgz */decorated*

tar -cvfz --ignore-failed-read --hard-dereference HTMLReport.tgz external/ out* PhyloCommunity.js *.html   *.mibybranch *.TreeLabeled
