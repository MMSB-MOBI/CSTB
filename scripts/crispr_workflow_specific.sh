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


error_json () {
    echo "{\"emptySearch\": \"There is a problem, impossible to finish the program\"}" > fail.log
    cat fail.log
}

parse_logFile () {
    if grep "Program terminated" $1 > /dev/null;
    then
        perl -ne 'if ($_ =~ /Program terminated/){
            @error_split=split(/&/);
            $msg = $error_split[1];
            $msg =~ s/\n$//;
            print "{\"emptySearch\" :  \"$msg\" }";
        }' $1 > ./fail.log;
        cat ./fail.log;
        false
    else
         true
    fi
}

if [ "$pam" != "NGG" ]; then
    error_json

#elif [ "$sl" != "20" ]; then
#    error_json

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

    ## CREATE METAFILE ###
    queryFasta="query.fasta"
    printf ">query\n$seq\n" > $queryFasta

    metafileQuery="query"
    echo python -u $CRISPR_TOOL_SCRIPT_PATH/create_metafile.py -file $queryFasta -out $metafileQuery > sg.cmd
    python -u $CRISPR_TOOL_SCRIPT_PATH/create_metafile.py -file $queryFasta -out $metafileQuery 2>> ./metafile.err 1>./metafile.log
    # Check if sgRNA on the query gene
    parse_logFile ./metafile.log
    PRG_TERMINATED=$?
    if [ $PRG_TERMINATED = 0 ];then
        ### FILTER GENOMES ###
       # gi=$(python $CRISPR_TOOL_SCRIPT_PATH/filter_specie.py --ref $SPECIE_REF_JSON --query "$gi")
       # gni=$(python $CRISPR_TOOL_SCRIPT_PATH/filter_specie.py --ref $SPECIE_REF_JSON --query "$gni")
        echo $gi > f.gi

        ### INTERSECTION ###
        fileSet="set_index.txt"
        setCompare $slFlag -i "$gi" -o "$gni" -l $rfg -e index -f $fileSet -s $metafileQuery".index" 2>> ./setCompare.err 1> ./setCompare.log

        nb_hits=`sed -n "s/^# \([0-9]*\).*/\1/p" $fileSet`
    fi
    # Check if hits in intersection
    if [ $nb_hits = 0 ];then
        echo "{\"emptySearch\" :  \"No common hits between the sequence and included genomes \"}"  > ./fail.log
        cat ./fail.log
        PRG_TERMINATED=1
    fi

    ### BLAST N TO FIND HOMOLOGOUS GENES ###
    if [ $PRG_TERMINATED = 0 ];then
        fileBlast="blast_output.xml"
        echo blastn -outfmt 5 -query $queryFasta -db $blastdb > $fileBlast >> sg.cmd
        blastn -outfmt 5 -query $queryFasta -db $blastdb > $fileBlast

        ### PARSE BLAST OUTPUT ###
        parseBlast="parse_blast.p"
        echo python -u $CRISPR_TOOL_SCRIPT_PATH/parse_blast.py -blast $fileBlast -gi "$gi" -o $parseBlast -ip $pid >> sg.cmd
        python -u $CRISPR_TOOL_SCRIPT_PATH/parse_blast.py -blast $fileBlast -gi "$gi" -o $parseBlast -ip $pid 2>> ./parse_blast.err 1> ./parse_blast.log
        # check if hits
        parse_logFile ./parse_blast.log
        PRG_TERMINATED=$?
    fi

    ### CHECK IG SGRNA ARE ON HOMOLOGOUS GENES ###
    if [ $PRG_TERMINATED = 0 ];then
        # Post-processing with setCompare output and blast output
        echo python -u $CRISPR_TOOL_SCRIPT_PATH/specific_gene.py -f $fileSet -sl $sl -pam "NGG" -gi "$gi" -gni "$gni" -r "$URL_CRISPR" -taxon_db "$NAME_TAXON" -tree_db "$NAME_TREE" -end_point "$URL_TREE_TAXON" -c 2000 --no-proxy -blast $parseBlast >> sg.cmd
        python -u $CRISPR_TOOL_SCRIPT_PATH/specific_gene.py -f $fileSet -sl $sl -pam "NGG" -gi "$gi" -gni "$gni" -r "$URL_CRISPR"  -taxon_db "$NAME_TAXON" -tree_db "$NAME_TREE" -end_point "$URL_TREE_TAXON" -c 2000 --no-proxy -blast $parseBlast 2>> ./specific_gene.err 1> ./specific_gene.log
        # Check if exist sgrna on genes
        parse_logFile ./specific_gene.log
        PRG_TERMINATED=$?
    fi

    if [ $PRG_TERMINATED = 0 ];then
        not_in=$(perl -ne 'BEGIN{$NR=0};$NR++; if($NR == 1){chomp;$_ =~ s/^[\s]*([\S].*[\S])[\s]*$/$1/;print $_; exit;}' ./specific_gene.log);
        number_hits=$(perl -ne 'BEGIN{$NR=0};$NR++; if($NR == 3){chomp;$_ =~ s/^[\s]*([\S].*[\S])[\s]*$/$1/;print $_; exit;}' ./specific_gene.log);
        tag=$(perl -ne 'BEGIN{$NR=0};$NR++; if($NR == 2){chomp;$_ =~ s/^[\s]*([\S].*[\S])[\s]*$/$1/;print $_; exit;}' ./specific_gene.log);

        echo "$not_in" > ./stuff.log;
        echo "$number_hits" >> ./stuff.log;
        echo "$tag" >> ./stuff.log;
        loc=$(pwd | perl -ne '@tmp = split(/\//, $_); print "$tmp[$#tmp - 1]/$tmp[$#tmp]";');
        echo "{\"out\" : {\"data_card\": $(cat ./results_by_org.json), \"data\" : $(cat ./results.json),  \"not_in\" : \""$not_in"\",\"gi\" : \""$gi"\",  \"gene\" : $(cat ./genes.json), \"number_hits\" : \""$number_hits"\", \"tag\" : \""$loc"\", \"size\" : $(cat ./size_org.json)}}"

    fi
fi
