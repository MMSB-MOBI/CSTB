#!/bin/bash

echo python -u $CRISPR_TOOL_SCRIPT_PATH/add_scripts/check_taxonomy.py -taxid $TAXID -gcf $GCF -url $URL_TAXON -dbName $DB_TAXON_NAME 2> ./check_taxon.err 1> ./check_taxon.log
python -u $CRISPR_TOOL_SCRIPT_PATH/add_scripts/check_taxonomy.py -taxid $TAXID -gcf $GCF -url $URL_TAXON -dbName $DB_TAXON_NAME 2> ./check_taxon.err 1> ./check_taxon.log

if grep "Be careful" ./check_taxon.log > /dev/null;
then
    perl -ne 'if ($_ =~ /Be careful/){
        @error_split=split(/&/);
        $msg = $error_split[1];
        $msg =~ s/\n$//;
        print "{\"becareful\" :  \"$msg\" }";
    }' ./check_taxon.log > ./careful.log;
    cat ./careful.log;
elif grep "Program terminated" ./check_taxon.log > /dev/null;
then
  perl -ne 'if ($_ =~ /Program terminated/){
      @error_split=split(/&/);
      $msg = $error_split[1];
      $msg =~ s/\n$//;
      print "{\"emptySearch\" :  \"$msg\" }";
  }' ./check_taxon.log > ./emptySearch.log;
  cat ./emptySearch.log;
else
  echo "{\"new\" : \" Check done : New taxon \"}"
fi
