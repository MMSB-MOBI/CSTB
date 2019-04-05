#!/bin/bash

error_json () {
    echo "{\"emptySearch\": \"There is a problem, impossible to finish the program\"}" > fail.log
    cat fail.log
}

check_prg_no_terminated () {
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

    # Create Metafile
    queryFasta="query.fasta"
    echo ">query
$seq" > $queryFasta

    metafileQuery="query"
    echo python -u $CRISPR_TOOL_SCRIPT_PATH/create_metafile.py -file $queryFasta -out $metafileQuery > sg.cmd
    python -u $CRISPR_TOOL_SCRIPT_PATH/create_metafile.py -file $queryFasta -out $metafileQuery 2>> ./metafile.err 1>./metafile.log

    check_prg_no_terminated ./metafile.log
    PRG_TERMINATED=$?
    if [ $PRG_TERMINATED = 0 ];then
        # Filter genomes
        gi=$(python $CRISPR_TOOL_SCRIPT_PATH/filter_specie.py --ref $SPECIE_REF_JSON --query "$gi")
        gni=$(python $CRISPR_TOOL_SCRIPT_PATH/filter_specie.py --ref $SPECIE_REF_JSON --query "$gni")
        echo $gi > f.gi

        # Set Compare
        fileSet="set_index.txt"
        setCompare $slFlag -i "$gi" -o "$gni" -l $rfg -e index -f $fileSet -s $metafileQuery".index" 2>> ./setCompare.err 1> ./setCompare.log

        # Blast N to find homologous genes
        fileBlast="blast_output.xml"
        dbBlast="/data/databases/mobi/crispr_clean/big.fasta"
        echo blastn -outfmt 5 -query $queryFasta -db $dbBlast > $fileBlast >> sg.cmd
        blastn -outfmt 5 -query $queryFasta -db $dbBlast > $fileBlast

        # Parse the blast output
        parseBlast="parse_blast.p"
        echo python -u $CRISPR_TOOL_SCRIPT_PATH/parse_blast.py -blast $fileBlast -gi "$gi" -o $parseBlast -ip $pid >> sg.cmd
        python -u $CRISPR_TOOL_SCRIPT_PATH/parse_blast.py -blast $fileBlast -gi "$gi" -o $parseBlast -ip $pid 2>> ./parse_blast.err 1> ./parse_blast.log

        check_prg_no_terminated ./parse_blast.log
        PRG_TERMINATED=$?
    fi

    if [ $PRG_TERMINATED = 0 ];then
        # Post-processing with setCompare output and blast output
        echo python -u $CRISPR_TOOL_SCRIPT_PATH/specific_gene.py -f $fileSet -sl $sl -pam "NGG" -gi "$gi" -gni "$gni" -r "$URL_CRISPR"  -c 2000 --no-proxy -blast $parseBlast >> sg.cmd
        python -u $CRISPR_TOOL_SCRIPT_PATH/specific_gene.py -f $fileSet -sl $sl -pam "NGG" -gi "$gi" -gni "$gni" -r "$URL_CRISPR"  -c 2000 --no-proxy -blast $parseBlast 2>> ./specific_gene.err 1> ./specific_gene.log

        check_prg_no_terminated ./specific_gene.log
        PRG_TERMINATED=$?
    fi

    if [ $PRG_TERMINATED = 0 ];then
        not_in=$(perl -ne 'chomp;$_ =~ s/^[\s]*([\S].*[\S])[\s]*$/$1/;print $_; exit;' ./specific_gene.log);
        number_hits=$(perl -ne 'BEGIN{$NR=0};$NR++; if($NR == 3){chomp;$_ =~ s/^[\s]*([\S].*[\S])[\s]*$/$1/;print $_; exit;}' ./specific_gene.log);
        tag=$(perl -ne 'BEGIN{$NR=0};$NR++; if($NR == 2){chomp;$_ =~ s/^[\s]*([\S].*[\S])[\s]*$/$1/;print $_; exit;}' ./specific_gene.log);
        echo "$not_in" > ./stuff.log;
        echo "$number_hits" >> ./stuff.log;
        echo "$tag" >> ./stuff.log;
        loc=$(pwd | perl -ne '@tmp = split(/\//, $_); print "$tmp[$#tmp - 1]/$tmp[$#tmp]";');
        echo "{\"out\" : {\"data\" : $(cat ./results.json),  \"not_in\" : \""$not_in"\",  \"number_hits\" : \""$number_hits"\", \"tag\" : \""$loc"\"}}"
    fi
fi
