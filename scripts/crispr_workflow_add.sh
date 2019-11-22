#!/bin/bash

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


if [ "$URL_TAXON" = "" ]; then
    error_json
elif [ "$CRISPR_TOOL_SCRIPT_PATH" = "" ]; then
    error_json
elif  [ "$URL_CRISPR" = "" ]; then
    error_json
elif [ "$DB_TAXON_NAME" = "" ]; then
    error_json
elif [ "$FOLDER" = "" ]; then
    FOLDER=$(echo "$TAXID""_pickle")
    mkdir $FOLDER
else

    ### CREATE PICKLE AND INDEX METAFILE ###
    if [ $PLASMID = 0 ]; then
      argPlas="-plasmid"
    else
      argPlas=""
    fi
    echo python -u $CRISPR_TOOL_SCRIPT_PATH/create_metafile.py -file $FASTA_FILE -taxid $TAXID -gcf $GCF -rfg $rfg $argPlas > cmd.log
    python -u $CRISPR_TOOL_SCRIPT_PATH/create_metafile.py -file $FASTA_FILE -taxid "$TAXID" -gcf $GCF -rfg $rfg $argPlas 2> ./create_meta.err 1> ./create_meta.log
    echo cp $FASTA_FILE $rfg"/genome_fasta/"$TAXID"_"$GCF".fna" >> cmd.log
    cp "$FASTA_FILE" "$rfg/genome_fasta/$TAXID""_""$GCF.fna"

    # Check if any sgRNA has been found in this genome
    parse_logFile ./create_meta.log
    PRG_TERMINATED=$?

    ### ADD INTO DATABASE AND UPDATE JSON_TREE ###
    if [ $PRG_TERMINATED = 0 ]; then

        ## Create pickle file to insert into taxon_db ##
        echo python -u $CRISPR_TOOL_SCRIPT_PATH/add_scripts/create_file_taxondb.py single -gcf $GCF -taxid $TAXID -r $URL_TAXON -dbName $DB_TAXON_NAME >> cmd.log
        python -u $CRISPR_TOOL_SCRIPT_PATH/add_scripts/create_file_taxondb.py single -gcf $GCF -taxid "$TAXID" -r $URL_TAXON -dbName $DB_TAXON_NAME 2> ./taxondb_file.err 1> ./taxondb_file.log
        mv ./taxonDB_data $FOLDER/taxonDB_data

        parse_logFile ./taxondb_file.log
        PRG_TERMINATED=$?
        if [ $PRG_TERMINATED = 0 ]; then
          ## Update the Json_tree ##
          if [ $PLASMID = 0 ]; then
            argPlas="-name"
          else
            argPlas="-taxid"
          fi
          echo python -u $CRISPR_TOOL_SCRIPT_PATH/add_scripts/update_tree.py -url $URL_TREE -treeName $DB_TREE_NAME -taxonDB $URL_TAXON -taxonName $DB_TAXON_NAME $argPlas "$TAXID" --no-proxy>> cmd.log
          python -u $CRISPR_TOOL_SCRIPT_PATH/add_scripts/update_tree.py -url $URL_TREE -treeName $DB_TREE_NAME -taxonDB $URL_TAXON -taxonName $DB_TAXON_NAME $argPlas "$TAXID"  --no-proxy 2> ./update_tree.err 1> ./update_tree.log
          mv ./treeDB_data $FOLDER

          # Check if no problem to connect to the database (taxon_db and taxon_tree_db), to access to GCF
          parse_logFile ./update_tree.log
          PRG_TERMINATED=$?
          folder=`pwd`
          mail -s "[CRISPR] : Add new genome" sophie.lematre@ibcp.fr <<< "The genome with the taxid $TAXID  and the GCF  $GCF  is ready.        $folder"
        fi
    fi
fi




# ****************************
# *                          *
# * TO INSERT INTO DATABASES *
# *                          *
# ****************************
# NAME_FILE=`cat ./create_meta.log`
# ## Add the new genome into the CRISPR database
# echo python -u $CRISPR_TOOL_SCRIPT_PATH/add_scripts/couchBuild.py --url $URL_CRISPR --map $MAP_FILE --data $rfg"/genome_pickle/"$NAME_FILE".p" > cmd.log
# python -u $CRISPR_TOOL_SCRIPT_PATH/add_scripts/couchBuild.py --url $URL_CRISPR --map $MAP_FILE --data $rfg"/genome_pickle/"$NAME_FILE".p"
# if [ -f error_add_db.err ]; then
#     echo "{\"emptySearch\": \"There is a problem during the insertion into CRISPR database. Contact us \"}" > fail.log
#     cat fail.log
# else


# ## Add into taxon database ##
# echo python -u $CRISPR_TOOL_SCRIPT_PATH/add_scripts/couchBuild.py $DB_TAXON_NAME --url $URL_TAXON --data "$FOLDER/taxonDB_data/" >> cmd.log
# python -u $CRISPR_TOOL_SCRIPT_PATH/add_scripts/couchBuild.py $DB_TAXON_NAME --url $URL_TAXON --data "$FOLDER/taxonDB_data/"
#
# if [ $PRG_TERMINATED = 0 ]; then
#   echo python -u $CRISPR_TOOL_SCRIPT_PATH/add_scripts/couchBuild.py $DB_TREE_NAME --url $URL_TREE --data "$FOLDER/treeDB_data/" >> cmd.log
#   python -u $CRISPR_TOOL_SCRIPT_PATH/add_scripts/couchBuild.py $DB_TREE_NAME --url $URL_TREE --data "$FOLDER/treeDB_data/"
# fi
