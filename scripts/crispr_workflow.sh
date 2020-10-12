#!/bin/bash

# pam="NGG"
# CRISPR_TOOL_SCRIPT_PATH="../crispr/bin"
# URL_CRISPR="http://localhost:2346"
# sl="20"
# fileSet="../crispr/test/data/ag_simple/set_index.txt"
# gi="Buchnera aphidicola (Cinara tujafilina) GCF_000217635.1&Aliivibrio wodanis GCF_000953695.1"
# gni=""
# rfg="../reference_genomes_pickle/"
# NAME_TAXON="taxon_db"
# NAME_TREE="taxon_tree"
# URL_TREE_TAXON="http://localhost:2346/"

run_setCompare() {
    slFlag=""
    if [ "$sl" != "20" ]; then
        ((_sl = $sl + 3))
        slFlag="-c ${_sl}"
    fi
    fileSet="set_index.txt"
    echo setCompare $slFlag -i "$gi" -o "$gni" -l $rfg -e index -f $fileSet > setCompare.cmd
    setCompare $slFlag -i "$gi" -o "$gni" -l $rfg -e index -f $fileSet 2> ./setCompare.err 1> ./setCompare.log
}

run_post_processing(){
    if [[ -s $fileSet ]]; then
        loc=$(pwd | perl -ne '@tmp = split(/\//, $_); print "$tmp[$#tmp]";');
        echo python -u $CRISPR_TOOLS_SCRIPT_PATH/post_processing.py --include "$gi" --exclude "$gni" --couch_endpoint "$COUCH_ENDPOINT" --taxon_db "$NAME_TAXON" --genome_db "$NAME_GENOME" --set_compare set_index.txt --length "$sl" --motif_broker_endpoint "$MOTIF_BROKER_ENDPOINT" --tag "$loc" > post_processing.cmd
        python -u $CRISPR_TOOLS_SCRIPT_PATH/post_processing.py --include "$gi" --exclude "$gni" --couch_endpoint "$COUCH_ENDPOINT" --taxon_db "$NAME_TAXON" --genome_db "$NAME_GENOME" --set_compare set_index.txt --length "$sl" --motif_broker_endpoint "$MOTIF_BROKER_ENDPOINT" --tag "$loc" 2>> post_processing.err
    fi
}

error_json () {
    msg=$1
    echo "{\"error\": \"$1\"}"
}

empty_json() {
    msg=$1
    echo "{\"emptySearch\": \"$1\"}"
}

SECONDS=0
run_setCompare
if [[ -s setCompare.err ]]; then
   msg=$(cat setCompare.err)
	if [[ $msg == "intersect size is zero"* ]]; then
		empty_json "No hits found"
    else
        error_json "Error while setCompare - job $loc. Contact support with this job number : cecile.hilpert@ibcp.fr"
	fi
else 
    run_post_processing
fi

echo $SECONDS > running_time.txt



