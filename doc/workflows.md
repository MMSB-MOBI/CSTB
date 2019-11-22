# CRISPR workflow
![Worklofw script](https://github.com/sophielem/crispr/blob/dev_add_genome/doc/workflow_script.png)

First of all, use *setCompare.c* script to make a dichotomous search from index files.

## crispr workflow

```sh
usage: post_processing.py [-h] [-sl <int>] -pam <str> -gi <str> -gni <str> -r
                          <str> -taxon_db <str> -tree_db <str> -end_point
                          <str> -c <int> -f <str> [--no-proxy] [-nb_top <int>]

Post-processing results

optional arguments:
  -h, --help        show this help message and exit
  -sl <int>         The length of the sgrna, excluding pam
  -pam <str>        The pam motif
  -gi <str>         The organisms to search inclusion in.
  -gni <str>        The organisms to search exclusion from
  -r <str>          The end point
  -taxon_db <str>   The name of the taxon database
  -tree_db <str>    The name of the taxon database
  -end_point <str>  The end point of the taxon and tree database
  -c <int>          The length of the slice for the request
  -f <str>          The file index
  --no-proxy
  -nb_top <int>     The top hits to download
```

```sh
python bin/post_processing.py -f "test/data/ag_simple/set_index.txt" -sl 20 -pam "NGG" -gi "Buchnera aphidicola (Cinara tujafilina) GCF_000217635.1&Aliivibrio wodanis GCF_000953695.1"  -gni "" -r "http://localhost:2346/" -taxon_db "taxon_db" -tree_db "taxon_tree" -end_point "http://localhost:2346/" -c 2000 --no-proxy
```

## crispr workflow specific
Retrieve sgRNA on a specific gene given by user. Retrieve result from SetCompare script and
create files for displaying on webservice.

```sh
usage: specific_gene.py [-h] -gi <str> -gni <str> [-sl [<int>]] [-pam [<str>]]
                        -f <str> -blast <str> [-nb_top <int>] [--no-proxy] -r
                        <str> -c <int> -taxon_db <str> -tree_db <str>
                        -end_point <str>

Specific gene program.

optional arguments:
  -h, --help        show this help message and exit
  -gi <str>         The organisms to search inclusion in.
  -gni <str>        The organisms to search exclusion from
  -sl [<int>]       The length of the sgrna, excluding pam
  -pam [<str>]      The pam motif
  -f <str>          The file index
  -blast <str>      The blast output file
  -nb_top <int>     The top hits to download
  --no-proxy
  -r <str>          The end point
  -c <int>          The length of the slice for the request
  -taxon_db <str>   The name of the taxon database
  -tree_db <str>    The name of the taxon database
  -end_point <str>  The end point of the taxon and tree database
```

```sh
python bin/specific_gene.py -f "test/data/sg/output_c.txt" -sl 20 -pam "NGG" -gi "Enterobacter sp. 638 GCF_000016325.1&Candidatus Blochmannia vafer str. BVAF GCF_000185985.2" -gni "" -r "http://localhost:2346/" -taxon_db "taxon_db" -tree_db "taxon_db" -end_point "http://localhost:2346/" -c 2000 --no-proxy -blast "parse_blast.p"
```
### Parse blastn output
If no database for blast software, use this command to create it:
```sh
makeblastdb -in file.fasta -title "Database with all genomes" -dbtype nucl
```

Parse the output of a Blastn on the reference database. Create a BlastReport object which
contains only hits for organism which are in a given list, here the included genomes list, and
hits which have a percentage identity superior to a given percentage identity, by default 70.
The BlastReport contains a dictionary of organism with references which contains BlastHit object
This object contains the coordinates of the hit and its length.
```sh
usage: parse_blast.py [-h] -blast <str> [-ip [IP]] -gi <str> [-o [<str>]]

Parse the output of blastn

optional arguments:
  -h, --help    show this help message and exit
  -blast <str>  The path to the blastn output file
  -ip [IP]      identity percentage min for the research of homologous genes
                using blastn (default:70)
  -gi <str>     The organisms to search inclusion in.
  -o [<str>]    The output path for the pickle file
 ```

## Add genome
![Schema workflow add new genome](https://github.com/sophielem/crispr/blob/dev_add_genome/doc/add_new_genome.png)

### Check taxonomy
Check if a taxonomy ID given exists in the NCBI database. Then, check if this taxonomy ID is
present in a Taxonomy database and check if the GCF given for this taxonomy ID already exists in the Taxonomy database.
Print a message with the key word __becareful__ if the taxonomy is already in the database, __Program terminated__ if the ID does not exist in the NCBI database, or if the GCF is the current one. Otherwise, no message.

```sh
usage: check_taxonomy.py [-h] -taxid <str> -gcf <str> -r <str> -dbName <str>
                         -plasmid PLASMID

Check if taxonomy information given by user are correct

optional arguments:
  -h, --help        show this help message and exit
  -taxid <str>      The taxonomy ID
  -gcf <str>        The GCF associated to the taxonomy ID
  -r <str>          End point of the taxon database with slash
  -dbName <str>     Name of the taxon database
  -plasmid PLASMID  The genome is a plasmid
```

### Update tree
Load the MaxiTree object from the taxon_tree_db and insert a new member.
Then, create a pickle file to insert it into the database.

```sh
usage: update_tree.py [-h] -url <str> -treeName <str> -taxonDB <str>
                      -taxonName <str> [-taxid <str>] [-name <str>]
                      [--no-proxy]

Update the MaxiTree object from the database

optional arguments:
  -h, --help        show this help message and exit
  -url <str>        The endpoint of the tree database
  -treeName <str>   Name of the tree database where the MaxiTree object is
                    saved
  -taxonDB <str>    The endpoint to the taxon database
  -taxonName <str>  Name of the taxon database
  -taxid <str>      Taxonomy ID to update or to add
  -name <str>       Name of plasmid to add
  --no-proxy
```

## Date
July 25 2019
