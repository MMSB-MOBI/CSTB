This readme is obsolete

# CRISPR : search common sgrna sequences to target bacterial strains
This script is a part of a pipeline which allows to find common sgrna from a set of genomes of bacteria which cannot match with a set of other genomes.

With this, it is possible to select an heterogeneous population of bacteria and to kill a second population.
![Schema CRISPR service](https://github.com/sophielem/crispr/blob/dev_add_genome/doc/schema_crispr_service.png)

## Dependencies
To execute this script, you need few dependencies :
* docopt&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;== 0.6.2
* ete3&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;== 3.1.1
* biopython == 1.66

You also need blastn.

## Results files
Only genomes selected by user are conserved and results are displayed by two *json* file containing the first 100 hits organized by sgRNA or by organism and a *text* file containing the first 1,000 hits.
Another *json* file contains results sorted by organism name and not by sgRNAs.
##### Format of *json* files :
Organized by sequences
```json
{'sequence' : word, 'occurences' :
                                    {'org' : genome, 'all_ref' :
                                                                {'ref' : ref, 'coords' : [coordinates]
                                                               }
                                    }
}
```

Organized by organisms
```json
{org : {
          ref : {
                  sgRNA : [coordinates]
                }
        }
}
```

##### Format of the *text* file :
```
#ALL GENOMES
#Genomes included : *list of genomes*  ; Genomes excluded : *list of genomes*
#Parameters: pam: *PAM* ; sgrna size: *size*
        genome_in_1     genome_in_2     genomes_in_3...
word    ref : coordinates   ref : coordinates   ref : coordinates...
.
.
.
```

## Databases
[ReadMe to launch server](https://github.com/sophielem/crispr/blob/dev_add_genome/doc/README.md)

### Plain text
In a folder *crispr_clean*, index files are saved in *genome_index* folder on arwen-dev cluster and fasta files are saved in *genome_fasta* folder. Pickle files are saved in *genome_pickle* in the folder *crispr_clean* on arwen. These files are not needed for scripts, they are just a backup of the CRISPR database.
The database for BLAST is in *crispr_clean* on arwen-dev.
*genome_index* : name + space + GCF + .index
*genome_pickle* : name + space + GCF + .p
*genome_fasta* : GCF + undescrore + ASM + tar.gz

### CRISPR database
This database contains several volumes, 64. Each volume can be accessed from a regex, which is the length prefix 3 of the sgRNA. The entry is the sgRNA sequence and it contains a dictionary of organism and coordinates under the key **data**.

This database is full with pickle file and the *couchBuild.py* script.

#### Create Metafile (*pickle* and *index*)
```sh
usage: create_metafile.py [-h] -file <str> [-out <str>] [-taxid <str>]
                          [-gcf <str>] [-rfg [<str>]] [-plasmid]

Create pickle and index metafile

optional arguments:
  -h, --help    show this help message and exit
  -file <str>   The fasta file to parse
  -out <str>    Path for the output file without extension
  -taxid <str>  Taxonomy ID
  -gcf <str>    GCF ID
  -rfg [<str>]  The path to the database for index and pickle file
  -plasmid      If present, indicates the genome is a plasmid
```

##### Word detection
Detect sgRNAs by regex and dump results in a pickle file.

##### Indexing ATCG words

```sh
python wordIntegerIndexing.py data/example.pickle --out /any/folder/toto.index
```

will produce `/any/folder/toto.index`

!!! **Otherwise output file name is guessed from input and wrote locally**

```sh
python wordIntegerIndexing.py data/example.pickle
```

will produce `./example.index`

### Taxon database
It contains entries with taxonomy ID for name. This document contains the list of GCF under the key **GCF** and the current GCF under key **current**. The keys **date** and **user** give the date and the user who did the last update.

```JSON
{
  "_id": "100226",
  "_rev": "7-c85bed6b4690fa9f8c456490df9d3895",
  "GCF": [
    "GCF_000203835.1"
  ],
  "date": "2019-07-17 15:35",
  "user": "MMSB",
  "current": "GCF_000203835.1",
  "size": {
    "NC_003888.3": 8667507
  }
}
```

#### Create taxon file
```sh
usage: create_file_taxondb.py [-h] {scratch,single} ...

Convert the genomre_ref_taxid.json to a pickle file ready to be insert into
database

positional arguments:
  {scratch,single}  commands
    scratch         Create files to insert from the genomre_ref_taxid.json
                    file
    single          Create file to insert

optional arguments:
  -h, --help        show this help message and exit
```

Create a taxon file for a single genome. The url is needed to check if an entry exists for this taxon ID. If so, it takes the list of GCF, add the new GCF and make it the current one.
```sh
usage: create_file_taxondb.py single [-h] [-user [<str>]] -gcf <str> -taxid
                                     <str> -r <str> -dbName <str>
                                     [-out [<str>]] -rfg <str>

optional arguments:
  -h, --help     show this help message and exit
  -user [<str>]  The user
  -gcf <str>     GCF id
  -taxid <str>   Taxonomy id
  -r <str>       End point for taxon Database
  -dbName <str>  Name of the taxon Database
  -out [<str>]   The name of the outputfile
  -rfg <str>     Path of the fasta file
```

Convert the genome_ref_taxid.json to a pickle file ready to be insert into
database or list file with the format : taxonId \t GCF_id
```sh
usage: create_file_taxondb.py scratch [-h] -file <str> [-user [<str>]]
                                      [-out [<str>]] -rfg <str> [--no-proxy]
                                      [--json]

optional arguments:
  -h, --help     show this help message and exit
  -file <str>    The json file to convert or the list file in csv
  -user [<str>]  The user
  -out [<str>]   The name of the outputfile
  -rfg <str>     Path of the fasta file
  --no-proxy
  --json
```

The output will be :</br>
taxon_dt={"1234": {"GCF": ["GCF_111", "GCF_112", "GCF_113"], "date": "19-04-2019", "User": "Queen", "current": "GCF_111", "size":{"NC_123": 20003}},</br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"2345": {"GCF": ["GCF_211", "GCF_212", "GCF_213"], "date": "19-04-2019", "User": "Queen", "current": "GCF_211", "size":{"NC_132": 74563}}}

### Tree database
It contains only one entry, **maxi_tree**. This document contains a key **tree** with the tree in JSON format ({'text' : blabla, 'children': ['text': blabla....]}) and the date of the update under key **date**.  

```JSON
{
  "_id": "maxi_tree",
  "_rev": "22-8f8a7e827a00ea6605ce2c7b8a8d8aed",
  "tree": "{'text': 'root', 'children': [{'text': 'Gammaproteobacteria', 'children': [{'text': 'Enterobacterales',...",
  "date": "2019-05-09 12:44"
}
```

#### Maxi Tree object
Class MaxiTree object. This tree contains: a tree object from the package ete3 with the name
of the organism, the current GCF and its taxonomy Id if it's present into the CRISPR database. This taxon ID is
referenced in the name separate by a ':' and in a node's feature 'taxon'.
The name of a node is composed like this:
    name (without ' and /) GCF_ID : taxon_ID

Several constructor:
* from the database
* from the genome_ref_taxid.json file
* from a MaxiTree pickle file
* from the taxon database

Several methods are implemented to insert a new node, to write a json file of the tree
full (with taxonID) or not, to get the json string full or not...

The main function allows to create this Tree from the genome_ref_taxid.json file and to create a pickle MaxiTree file to insert into the taxon_tree_db.



## Workflows

[ReadMe for workflows](https://github.com/sophielem/crispr/blob/dev_add_genome/scripts/workflows.md)

## Date
June 24 2019
