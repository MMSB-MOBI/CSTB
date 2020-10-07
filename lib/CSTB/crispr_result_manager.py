import logging
import pycouch.error
import re
from collections import OrderedDict
from CSTB.engine.crispr_hit import Hit
import CSTB.utils.error as error
import requests
import time
from typing import List, Dict
from CSTB.engine.crispr_blast import BlastReport
import motif_broker_request.request as mb_request


SESSION = requests.Session()
SESSION.trust_env = False

'''TO DO
- make module that interrogate CSTB database (function get_taxon_name and get_genomes_size)
'''

class CrisprResultManager():
    """Manager for crispr results
    
    Attributes:
        wrapper (pycouch.wrapper.Wrapper): pycouch wrapper to interrogate couchDB
        motif_broker_endpoint (str): motif-broker endpoint to get sgRNA occurences from couchDB with mapping
        taxondb (str): taxon database name
        genomedb (str): genome database name
        tag (str): tag for result storage
        nb_total_hits (int): Number of hits found by setCompare. Set with parse_set_compare() call.
        nb_treated_hits (int): Number of setCompare hits treated by this manager. Set with parse_set_compare() call.
        hits_collection (List[Hit]) : List of Hit objects
        include_taxon (Dict) : Dict with included genomes uuid as keys and corresponding taxon name as value
        exclude_taxon (Dict) : Dict with included genomes uuid as keys and corresponding taxon name as value

    """

    def __init__(self, pycouch_wrapper, taxondb, genomedb, motif_broker_endpoint, tag, include_taxon = {}, exclude_taxon = {}, nb_total_hits = None, nb_treated_hits = None, hits_collection = [], homolog_genes = []):
        self.wrapper = pycouch_wrapper
        self.motif_broker_endpoint = motif_broker_endpoint
        self.taxondb = taxondb
        self.genomedb = genomedb
        self.tag = tag
        self.nb_total_hits = nb_total_hits
        self.nb_treated_hits = nb_treated_hits
        self.hits_collection = hits_collection
        self.include_taxon = include_taxon
        self.exclude_taxon = exclude_taxon
        self.homolog_genes = homolog_genes
        mb_request.configure(self.motif_broker_endpoint)

    #Move this in some consumer
    def get_taxon_name(self, list_uuid):
        """Get taxon names for a list of genome uuid. Will interrogate couchDB genome db and taxon db.
        
        Args:
            list_uuid (List[str]): list of genomes uuid
        
        Raises:
            error.CouchNotFound: Raise when couch document is not found
        """
        correspondance_genome_taxon = {}
        logging.debug(f"Get taxon name for {list_uuid}")
        for g_uuid in list_uuid:
            logging.debug(g_uuid)
            resp = self.wrapper.couchGetDoc(self.genomedb, g_uuid)
            if not resp:
                raise error.CouchNotFound(f"{self.wrapper.end_point}/{self.genomedb}/{g_uuid} not found")
            taxon_uuid = resp["taxon"]

            resp = self.wrapper.couchGetDoc(self.taxondb, taxon_uuid)
            if not resp:
                raise error.CouchNotFound(f"{self.wrapper.end_point}/{self.taxondb}/{taxon_uuid} not found")
            correspondance_genome_taxon[g_uuid] = resp["name"]

        return correspondance_genome_taxon

    #Move this in some consumer
    def get_genomes_size(self, list_uuid):
        """ Get genome size for a list of genomes uuid. Interrogate couchDB genome database.
        
        Args:
            list_uuid (List[str]): list of genomes uuid
        
        Raises:
            error.CouchNotFound: [description]
        
        Returns:
            Dict: Dict with genome_uuid as key and size dict as value. Size dict contains subsequence name as key and list of coordinates as value. 
        """
        sizes = {}
        for g_uuid in list_uuid:
            resp = self.wrapper.couchGetDoc(self.genomedb, g_uuid)
            if not resp:
                raise error.CouchNotFound(f"{self.wrapper.end_point}/{self.genomedb}/{g_uuid} not found")
            size = resp["size"]
            sizes[g_uuid] = size
        return sizes

    def get_fasta_metadata(self, list_uuid):
        """Get metadata for fasta subsequences (sizes and headers) of genomes given in list_uuid

        Args:
            list_uuid (str[]): [description]

        Raises:
            error.CouchNotFound: Raise if couch document doesn't exist
            error.FastaMetadataError: Raise if something is wrong with given list of genomes or couch document.

        Returns:
            Dict: { genome_uuid : { fasta_subsequence_ref : { "size" : size_value, "header" : header_value } } }. Dictionnary with header and size associated with each fasta subsequence of each genome from given list.
        """
        metadata = {}
        for g_uuid in list_uuid:
            resp = self.wrapper.couchGetDoc(self.genomedb, g_uuid)
            if not resp:
                raise error.CouchNotFound(f"{self.wrapper.end_point}/{self.genomedb}/{g_uuid} not found")
            
            if not resp.get("size"):
                raise error.FastaMetadataError(f"{self.wrapper.end_point}/{self.genomedb}/{g_uuid} has not 'size' attribute")
            
            if not resp.get("headers"):
                raise error.FastaMetadataError(f"{self.wrapper.end_point}/{self.genomedb}/{g_uuid} has not 'headers' attribute")
            if g_uuid in metadata:
                raise error.FastaMetadataError(f"{g_uuid} appears several times in list_uuid")
            metadata[g_uuid] = {}
            for fasta_seq in resp["size"]:
                if fasta_seq in metadata[g_uuid]:
                    raise error.FastaMetadataError(f"{fasta_seq} is duplicate in {self.wrapper.end_point}/{self.genomedb}/{g_uuid}")
                metadata[g_uuid][fasta_seq] = {}
                metadata[g_uuid][fasta_seq]["size"] = resp["size"][fasta_seq]
                metadata[g_uuid][fasta_seq]["header"] = resp["headers"][fasta_seq]
            
        return metadata
        
            
    def set_taxon_names(self, include, exclude):
        #type: (List[str], List[str]) -> None
        """
        Get taxon names for include and exclude genomes uuid. Call get_taxon_name() and set include_taxon and exclude_taxon attributes. 

        :param include: List of include genome uuid
        :param exclude: List of exclude genome uuid
        """

        if include:
            self.include_taxon = self.get_taxon_name(include)
        if exclude:
            self.exclude_taxon = self.get_taxon_name(exclude)

    def search_occurences(self, genomes_include):
        """Search sgRNA occurences in couchDB with motif broker and complete Hit objects.
        
        :param genomes_include: list of genome uuid
        :param len_slice: Length of packet to interrogate couchDB, defaults to 2000

        :raises error.PingError: Raise when motif-broker can't be reach
        """
        try:
            SESSION.get(self.motif_broker_endpoint + "/handshake")
        except:
            raise error.PingError(f"Can't handshake motif-broker at {self.motif_broker_endpoint}")
        
        all_seqs = [seq for hit in self.hits_collection for seq in hit.to_request_sequences]

        results = mb_request.get(all_seqs)
        
        for hit in self.hits_collection:
            hit.store_occurences([results[seq] for seq in hit.to_request_sequences], genomes_include)
            
           

    def parse_set_compare(self, setCompare_file, word_length, to_keep=None):
        """ Parse results from setCompare. Will call a different parsing function for word with size 20 and for shorter words. Separate functions are required because setCompare results format is not the same with with 20-length words and shorter words. Initialize hits_collection and nb_treated_hits
        
        :param setCompare_file: path to setCompare result file 
        :type setCompare_file: str
        :param word_length: Word length
        :type word_length: int
        :param to_keep: Number of hits to keep
        :type to_keep: int
        """
        
        
        self.hits_collection=[]
        logging.info(f"Parse set compare for word length {word_length}")
        if word_length == 20:
            self._parse_set_compare_20(setCompare_file, to_keep, word_length)
        else:
            self._parse_set_compare_other(setCompare_file, to_keep, word_length)
        logging.info(f"Nb hits collection {len(self.hits_collection)}")

        self.nb_treated_hits = len(self.hits_collection)
            
    def _parse_set_compare_20(self, setCompare_file, to_keep, word_length):
        """Parse setCompare for word of size 20. 
        
        :param setCompare_file: path to setCompare result file 
        :type setCompare_file: str
        :param to_keep: Number of hits to keep
        :type to_keep: int
        :param word_length: sgRNA length
        :type word_length: int
        """

        with open(setCompare_file, "r") as filin:
            text = filin.readlines()

        self.nb_total_hits = re.search("[0-9]+", text[-2]).group()
        if int(self.nb_total_hits) == 0:
            return

        index_dic = OrderedDict()
        i = 0
        for rank_occ in text[-1].strip().split(","):
            if to_keep and i == to_keep: break
            self.hits_collection.append(Hit(rank_occ.split(":")[0],rank_occ.split(":")[1], word_length + 3))
            #index_dic[int(rank_occ.split(":")[0])] = rank_occ.split(":")[1]
            i += 1

    def _parse_set_compare_other(self, setCompare_file, to_keep, word_length):
        logging.debug("parse set compare other")
        with open(setCompare_file, "r") as filin:
            for line in filin:
                logging.debug(line)
                regex_nb_hits = re.search("^# ([0-9]+)", line)
                if regex_nb_hits:
                    self.nb_total_hits = int(regex_nb_hits.group(1))

                    if self.nb_total_hits == 0:
                        return
                    break

            index_dic = OrderedDict()
            i = 0
            for rank_occ in filin:
                if (to_keep and i == to_keep) or rank_occ == "\n": break
                rank_splitted = rank_occ.split(":")
                
                index_sgrna = rank_splitted[0]
                weight = rank_splitted[1].split("[")[0]
                index_longer_sgrna = [int(index.replace("'","")) for index in rank_occ.split("[")[1].rstrip("]\n").split(",")]
                self.hits_collection.append(Hit(index_sgrna, weight, word_length + 3, longer_index = index_longer_sgrna))
                i += 1
            logging.debug(f"First hit\n{self.hits_collection[0]}")

    def generate_json_data(self):
        """Generate json data for client (to display into table) from Hit collection. 
        
        :return: json data, 
        :rtype: List[Dict]
        """

        results = []
        sorted_hits = sorted(self.hits_collection, key=lambda hit: hit.number_occurences)
        #logging.debug([hit.number_occurences for hit in self.hits_collection])
        for hit in sorted_hits:
            results.append({"sequence" : hit.sequence, "occurences" : hit.list_occ(self.include_taxon)})
        return results

    def generate_json_data_card(self):
        """Generate json data for client to display genomic card from Hit collection. 
        
        :return: json data
        :rtype: Dict { organism : { fasta_header : { sgRNA_sequence : List[str] } } }
        """

        def insert_in_data_card(genome_uuid, subseq, sequence, coords):
            genome_name = self.include_taxon[genome_uuid]
            if genome_name not in data_card:
                data_card[genome_name] = {}
            
            if subseq not in data_card[genome_name]:
                data_card[genome_name][subseq] = {}
            
            data_card[genome_name][subseq][hit.sequence] = hit.occurences[genome_uuid][subseq]

        data_card = {}
        for hit in self.hits_collection:
            for genome in hit.occurences:
                for subseq in hit.occurences[genome]:
                    insert_in_data_card(genome, subseq, hit.sequence, hit.occurences[genome][subseq])
        
        return data_card

    def generate_json_size(self):
        """[summary]
        
        :return: json data
        :rtype: Dict { organism : { fasta_header : length } }
        """
        json_size = {}
        sizes = self.get_genomes_size(list(self.include_taxon.keys()))
        for g_uuid in sizes:
            json_size[self.include_taxon[g_uuid]] = sizes[g_uuid]  
        
        return json_size

    def generate_fasta_metadata(self):
        """Format json for fasta metadata

        Returns:
            Metadata[]: List of Metadata. Metadata is dictionnary with organism_name, fasta reference, size and header. { "org" : organism_name, "fasta_ref" : fasta_reference, "size" : size value, "header" : header value }
        """
        json = []
        metadata = self.get_fasta_metadata(list(self.include_taxon.keys()))
        for g_uuid in metadata:
            g_name = self.include_taxon[g_uuid]
            for fasta_ref in metadata[g_uuid]:
                data = {"org" : g_name, \
                    "fasta_ref" : fasta_ref, \
                    "size" : metadata[g_uuid][fasta_ref]["size"], \
                    "header" : metadata[g_uuid][fasta_ref]["header"]}
                json.append(data)
        return json


        
    def format_results(self, blast = False):
        """Create the final json for client
        
        :return: final json
        :rtype: Dict
        """
        final_json = {}
        final_json["gi"] = "&".join(list(self.include_taxon.values()))
        final_json["not_in"] = ",".join(list(self.exclude_taxon.values()))
        final_json["number_hits"] = self.nb_total_hits
        final_json["number_treated_hits"] = self.nb_treated_hits
        final_json["data"] = self.generate_json_data()
        final_json["data_card"] = self.generate_json_data_card()
        final_json["tag"] = self.tag
        final_json["size"] = self.generate_json_size() #To delete
        final_json["fasta_metadata"] = self.generate_fasta_metadata()

        if blast:
            final_json["gene"] = self.generateGeneData()

        return final_json

    def parseBlast(self, blast_xml, identity, included_genomes):
        """Take blast xml file and parse it. Just keep homolog genes with more than defined identity and just keep coordinates that maps with homolog genes.
        
        :param blast_xml: path to blast xml result
        :type blast_xml: str
        :param identity: minimum identity percentage
        :type identity: int
        :param included_genomes: list of included genomes uuid
        :type included_genomes: List[str]
        :raises error.NoBlastHit: Raise if there is no blast hit at all
        :raises error.NoHomolog: Raise if at least 1 included genome has not an homolog gene.
        """
        blast_report = BlastReport(blast_xml, identity, included_genomes)
        if not blast_report.is_hit():
            raise error.NoBlastHit("No blast hit")

        logging.debug(f"Blast homolog genes {blast_report.homolog_genes}")

        # Check if all organisms have homolog genes
        if set(included_genomes) != blast_report.organisms:
            raise error.NoHomolog(set(included_genomes).difference(blast_report.organisms))

        self.homolog_genes = blast_report.filterGenes(included_genomes)
        
        for hit in self.hits_collection:
            hit.storeOnGeneOccurences(self.homolog_genes)

    def filterOnGeneOccurences(self):
        """Filter results to just keep hits with occurences on Gene. Return a new CrisprResultsManager object with a new hits collection.
        
        :return: new CrisprResultManager with same attributes as the current one, except for hit collection. 
        :rtype: CrisprResultManager
        """
        new_hit_collection = [hit for hit in self.hits_collection if hit.on_gene_occurences]
        new_results = CrisprResultManager(self.wrapper, self.taxondb, self.genomedb, self.motif_broker_endpoint, self.tag, self.include_taxon, self.exclude_taxon, self.nb_total_hits, self.nb_treated_hits, new_hit_collection, self.homolog_genes)
        return new_results

    def generateGeneData(self):
        """Generate json data for homolog genes when blast is used
        
        :return: json data
        :rtype: Dict -> {organism_name (str): {fasta_header (str) : [{'start' : gene.start(int), 'end' : gene.end(int)}] }}
        """
        json = {}
        for gene in self.homolog_genes:
            org_name = self.include_taxon[gene.org_uuid]
            if org_name not in json:
                json[org_name] = {}
            if gene.fasta_header not in json[org_name]:
                json[org_name][gene.fasta_header] = []
            
            json[org_name][gene.fasta_header].append({"start": gene.start, "end": gene.end})
        
        return json

    def serializeResults(self, file, gene = False):
        o = open(file,"w")
        o.write("#SgRNA sequence\tOrganism\tFasta sequence reference\tCoordinates")
        if gene:
            o.write("\tOn at least 1 homologous gene")
        o.write("\n")
        for hit in self.hits_collection:
            for org in hit.occurences:
                org_name = self.include_taxon[org]
                for fasta_seq in hit.occurences[org]: 
                    for coord in hit.occurences[org][fasta_seq]["coords"]:
                        o.write(f"{hit.sequence}\t{org_name}\t{fasta_seq}\t{coord['coord']}")
                        if gene :
                            if coord["is_on_gene"]:
                                o.write("\tTrue")
                            else:
                                o.write("\tFalse")
                        o.write("\n")
        o.close() 


        
    




