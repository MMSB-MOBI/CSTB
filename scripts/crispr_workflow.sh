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

error_json () {
    echo "{\"emptySearch\": \"There is a problem, impossible to finish the program\"}" > fail.log
    cat fail.log
}

if [ "$pam" != "NGG" ]; then
    error_json

elif [ "$CRISPR_TOOL_SCRIPT_PATH" = "" ]; then
    error_json

elif  [ "$URL_CRISPR" = "" ]; then
    error_json

else
    slFlag=""
    #  shorter word
    if [ "$sl" != "20" ]; then
        ((_sl = $sl + 3))
        slFlag="-d 23 -c ${_sl}"
    fi
    printenv > env.log

    BASE_FOLDER=`pwd`

    pwd > pwd.log

    #echo "curl -X GET $URL_CRISPR/handshake" > handshake.cmd
    #curl -X GET $URL_CRISPR/handshake &> handshake.log

    echo $gi > f.gi
     fileSet="set_index.txt"
     setCompare $slFlag -i "$gi" -o "$gni" -l $rfg -e index -f $fileSet 2>> ./setCompare.err 1> ./setCompare.log

    echo python -u $CRISPR_TOOL_SCRIPT_PATH/post_processing.py -f $fileSet -sl $sl -pam "NGG" -gi "$gi" -gni "$gni" -r "$URL_CRISPR" -taxon_db "$NAME_TAXON" -tree_db "$NAME_TREE" -end_point "$URL_TREE_TAXON" -c 2000 --no-proxy > pp.cmd
    python -u $CRISPR_TOOL_SCRIPT_PATH/post_processing.py -f $fileSet -sl $sl -pam "NGG" -gi "$gi" -gni "$gni" -r "$URL_CRISPR" -taxon_db "$NAME_TAXON" -tree_db "$NAME_TREE" -end_point "$URL_TREE_TAXON" -c 2000 --no-proxy 2>> ./post_processing.err 1> ./post_processing.log


    if grep "Program terminated" ./post_processing.log > /dev/null;
    then
    perl -ne 'if ($_ =~ /Program terminated/){
        @error_split=split(/&/);
        $msg = $error_split[1];
        $msg =~ s/\n$//;
        print "{\"emptySearch\" :  \"$msg\" }";
    }' ./post_processing.log > ./fail.log;
    cat ./fail.log;
    else
        not_in=$(perl -ne 'BEGIN{$NR=0};$NR++; if($NR == 1){chomp;$_ =~ s/^[\s]*([\S].*[\S])[\s]*$/$1/;print $_; exit;}' ./post_processing.log);
        number_hits=$(perl -ne 'BEGIN{$NR=0};$NR++; if($NR == 3){chomp;$_ =~ s/^[\s]*([\S].*[\S])[\s]*$/$1/;print $_; exit;}' ./post_processing.log);
        tag=$(perl -ne 'BEGIN{$NR=0};$NR++; if($NR == 2){chomp;$_ =~ s/^[\s]*([\S].*[\S])[\s]*$/$1/;print $_; exit;}' ./post_processing.log);
        echo "$not_in" > ./stuff.log;
        echo "$number_hits" >> ./stuff.log;
        echo "$tag" >> ./stuff.log;
        loc=$(pwd | perl -ne '@tmp = split(/\//, $_); print "$tmp[$#tmp - 1]/$tmp[$#tmp]";');
    	#number_hits=lines[2].strip()
        echo "{\"out\" : {\"data_card\": $(cat ./results_by_org.json), \"data\" : $(cat ./results.json),  \"not_in\" : \""$not_in"\",\"gi\" : \""$gi"\",  \"number_hits\" : \""$number_hits"\", \"tag\" : \""$loc"\", \"size\" : $(cat ./size_org.json)}}"
    fi


    #cp -r $WORKDIR/* $BASE_FOLDER
    #cd $BASE_FOLDER
fi
