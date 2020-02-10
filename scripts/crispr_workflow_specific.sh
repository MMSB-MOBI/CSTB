#!/bin/bash
# unset HTTP_PROXY
# unset HTTPS_PROXY
# pam="NGG"
# CRISPR_TOOL_SCRIPT_PATH="../crispr/bin"
# URL_CRISPR="http://localhost:2345"
# sl="20"
# fileSet="../crispr/test/data/sg/output_c.txt"
# gi="Enterobacter sp. 638 GCF_000016325.1&Candidatus Blochmannia vafer str. BVAF GCF_000185985.2"
# gni=""
# rfg="../reference_genomes_pickle/"
# URL_TAXON="http://localhost:5984/taxon_db_size"
# URL_TREE="http://localhost:5984/taxon_tree_db"
# fileBlast="../crispr/test/data/sg/blast.xml"
# pid=70

set -e

queryFasta="query.fasta"
printf ">query\n$seq\n" > $queryFasta
queryIndex="query.index"
#Search and index sgRNA in gene
echo python -u $CRISPR_TOOLS_SCRIPT_PATH/index_sequence.py -f $squeryFasta -o $queryIndex > index_query.cmd
python -u $CRISPR_TOOLS_SCRIPT_PATH/index_sequence.py -f $queryFasta -o $queryIndex 2> index_query.err

slFlag=""
if [ "$sl" != "20" ]; then
    ((_sl = $sl + 3))
    slFlag="-d 23 -c ${_sl}"
fi

fileSet="set_index.txt"
echo setCompare $slFlag -i "$gi" -o "$gni" -l $rfg -e index -f $fileSet -s $queryIndex > setCompare.cmd
setCompare $slFlag -i "$gi" -o "$gni" -l $rfg -e index -f $fileSet -s $queryIndex 2> ./setCompare.err 1> ./setCompare.log

nb_hits=`sed -n "s/^# \([0-9]*\).*/\1/p" $fileSet`

blastOutput="blast_output.xml"
echo blastn -outfmt 5 -query $queryFasta -db $blastdb > blast.cmd
blastn -outfmt 5 -query $queryFasta -db $blastdb > $blastOutput

echo python -u $CRISPR_TOOLS_SCRIPT_PATH/post_processing.py --include "$gi" --exclude "$gni" --couch_endpoint "$COUCH_ENDPOINT" --taxon_db "$NAME_TAXON" --genome_db "$NAME_GENOME" --set_compare set_index.txt --length "$sl" --motif_broker_endpoint "$URL_CRISPR" --tag "$loc" --blast $blastOutput > post_processing.cmd

python -u $CRISPR_TOOLS_SCRIPT_PATH/post_processing.py --include "$gi" --exclude "$gni" --couch_endpoint "$COUCH_ENDPOINT" --taxon_db "$NAME_TAXON" --genome_db "$NAME_GENOME" --set_compare set_index.txt --length "$sl" --motif_broker_endpoint "$URL_CRISPR" --tag "$loc" --blast $blastOutput 2> post_processing.err

