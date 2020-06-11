This repository contains the back-end module for CSTB web service. 
It contains 2 bash workflows to process the 2 features of CSTB web service : All Genomes and Specific Gene. 
This workflows need a couchDB database with flat index files. See [CSTB_database_manager](https://github.com/MMSB-MOBI/CSTB_database_manager) for construction. 

## All Genomes workflow 
The purpose of All Genomes workflow is to search sgRNAs present in a pool of organisms and absent from an other one. 
All Genomes workflow inputs are a list of included organisms, a list of excluded organisms and parameters for sgRNA research (PAM motif, SGRNA length). 

1. Intersection of sgRNAs 
First, setCompare module is used to keep sgRNAs present in included organisms and absent from excluded organisms. To do that, the module use int representation of sgRNAs for each organism and compute efficiently int intersection. The stored indexes represents 23bp length sgRNAs and the module can deduce shorter sgRNAs from them. 

2. Post-processing steps 
* Retrieve organism metadata  
The first step is to interrogate couchDB database to make the link between organism identifiant used as inputs (it's a couchDB uuid) and common names of organisms to display to user. 

* Parse setCompare results  
Convert setCompare results which are a list of int indexes into sgRNAs sequences 

* Search sgRNAs occurences  
With sgRNAs sequences, interrogate couchDB database to retrieve the coordinates of it inside each genome. 

* Format results   
Format worklow results into json to be read by CSTB client.

* Serialize results   
Create report of results as tsv flat file to be downloaded by user. 

## Specific Gene workflow 
Specific Gene workflow purpose is to search sgRNAs present in a gene and its homologous from a pool of bacteria and absent from an other pool of bacteria. 
Its inputs are the gene sequence, a list of included organisms, a list of excluded organisms, parameters for sgRNA (length, PAM motif) and parameters for homologous research (identity percentage). 

1. Search and encode sgRNA of the gene. 
Search sgRNAs in input gene and encode them with 2-bits int representation. 

2. Intersection of sgRNAs 
Intersect sgRNAs found in gene and all sgRNAs present in included organisms and absent from excluded organisms with setCompare. 

3. Identify homologous genes 
Use blast to identify homologous genes in included organisms. We works directly from nucleotides sequences to keep sgRNAs, so only very strong homology can be identified, which is not the general definition of homologous genes, more based on identity and similarity of protein sequence. We can call our homologous genes "near identical genes". 

4. Post-processing steps 
* Retrieve organism metadata  
Same as All Genomes. 

* Parse setCompare results   
Same as All Genomes. 

* Search sgRNAs occurences.   
Same as All Genomes but we also kept the information of presence or absence of sgRNA in homologous gene. 

* Format and serialize results.   

#### [Detailed documentation](https://mmsb-mobi.github.io/CSTB/)

## Dependencies 
* [setCompare](https://github.com/glaunay/crispr-set) : Cython library to intersect sgRNAs by using their int 2-bits encoded representation.
* [CSTB-core](https://github.com/MMSB-MOBI/CSTB_core) : Python library to decode and encode sgRNAS.
* [pyCouch](https://github.com/MMSB-MOBI/pyCouch) : Python library to interrogate couchDB database.
* [motif-broker](https://github.com/glaunay/motif-broker-2) : NodeJS service to efficiently retrieve sgRNAs from couch database. 
* [requests](https://pypi.org/project/requests/) : Python library to send http requests. Here, used to interrogate motif-broker service. 

## Web service parts 
* Back-end code (here) 
* [Server and interrogation web page](https://github.com/MMSB-MOBI/CSTB_server)
* [Result web page](https://github.com/MMSB-MOBI/result_page_crispr)
* [Database manager](https://github.com/MMSB-MOBI/CSTB_database_manager)
