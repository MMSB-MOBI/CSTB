import logging
import pycouch.error
import re
from collections import OrderedDict
from CSTB.engine.crispr_hit import Hit
import CSTB.utils.error as error
import requests
import time
from typing import List, Dict, TypedDict


SESSION = requests.Session()
SESSION.trust_env = False

'''TO DO
- make module that interrogate CSTB database (function get_taxon_name and get_genomes_size)
- make module with wordIntegerIndexing, also use in crispr-manager
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

    def __init__(self, pycouch_wrapper, taxondb, genomedb, motif_broker_endpoint, tag):
        self.wrapper = pycouch_wrapper
        self.motif_broker_endpoint = motif_broker_endpoint
        self.taxondb = taxondb
        self.genomedb = genomedb
        self.tag = tag
        self.nb_total_hits = None
        self.nb_treated_hits = None
        self.hits_collection = []
        self.include_taxon = {}
        self.exclude_taxon = {}

    #Move this in some consumer
    def get_taxon_name(self, list_uuid):
        """Get taxon names for a list of genome uuid. Will interrogate couchDB genome db and taxon db.
        
        Args:
            list_uuid (List[str]): list of genomes uuid
        
        Raises:
            error.CouchNotFound: Raise when couch document is not found
        """
        correspondance_genome_taxon = {}
        logging.debug(f"Try to get taxon name for {list_uuid}")
        for g_uuid in list_uuid:
            logging.debug(g_uuid)
            resp = self.wrapper.couchGetDoc(self.genomedb, g_uuid)
            if not resp:
                raise error.CouchNotFound(f"{self.wrapper.endpoint}/{self.genomedb}/{g_uuid} not found")
            taxon_uuid = resp["taxon"]

            resp = self.wrapper.couchGetDoc(self.taxondb, taxon_uuid)
            if not resp:
                raise error.CouchNotFound(f"{self.wrapper.endpoint}/{self.taxondb}/{taxon_uuid} not found")
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
                raise error.CouchNotFound(f"{self.wrapper.endpoint}/{self.genomedb}/{g_uuid} not found")
            size = resp["size"]
            sizes[g_uuid] = size
        return sizes
            
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

    def search_occurences(self, genomes_include, len_slice = 2000):
        #type: (List[str], int) -> None
        """Search sgRNA occurences in couchDB with motif broker and complete Hit objects.
        
        :param genomes_include: list of genome uuid
        :param len_slice: Length of packet to interrogate couchDB, defaults to 2000
        :raises error.PingError: Raise when motif-broker can't be reach
        """
        try:
            SESSION.get(self.motif_broker_endpoint + "/handshake")
        except:
            raise error.PingError(f"Can't handshake motif-broker at {self.motif_broker_endpoint}")

        results = {}
        for i in range(0, len(self.hits_collection), len_slice):
            joker = 0
            request_sliced = {"keys" : [hit.sequence for hit in self.hits_collection[i : i + len_slice]]}
            while True:
                try:
                    results.update(SESSION.post(self.motif_broker_endpoint + "/bulk_request",json=request_sliced).json()["request"])
                except:
                    joker += 1
                    if joker > 3:
                        raise Exception(f"Can't interrogate motif-broker at {self.motif_broker_endpoint} after 3 tries")
                    time.sleep(5)
                    continue
                break
        
        for hit in self.hits_collection:
            hit.store_occurences(results[hit.sequence], genomes_include)
           

    def parse_set_compare(self, setCompare_file, word_length, to_keep):
        """ Parse results from setCompare. Will call a different parsing function for word with size 20 and for shorter words. Separate functions are required because setCompare results format is not the same with with 20-length words and shorter words. Initialize hits_collection and nb_treated_hits
        
        :param setCompare_file: path to setCompare result file 
        :type setCompare_file: str
        :param word_length: Word length
        :type word_length: int
        :param to_keep: Number of hits to keep
        :type to_keep: int
        """
        
        self.nb_treated_hits = to_keep
        self.hits_collection=[]
        logging.info(f"Parse set compare for word length {word_length}")
        if word_length == 20:
            self._parse_set_compare_20(setCompare_file, to_keep)
        else:
            logging.debug("Handle this")
            
    def _parse_set_compare_20(self, setCompare_file, to_keep):
        """Parse setCompare for word of size 20. 
        
        :param setCompare_file: path to setCompare result file 
        :type setCompare_file: str
        :param to_keep: Number of hits to keep
        :type to_keep: int
        """

        with open(setCompare_file, "r") as filin:
            text = filin.readlines()

        self.nb_total_hits = re.search("[0-9]+", text[-2]).group()
        if int(self.nb_total_hits) == 0:
            return

        index_dic = OrderedDict()
        i = 0
        for rank_occ in text[-1].strip().split(","):
            if i == to_keep: break
            self.hits_collection.append(Hit(rank_occ.split(":")[0],rank_occ.split(":")[1]))
            #index_dic[int(rank_occ.split(":")[0])] = rank_occ.split(":")[1]
            i += 1

    def _parse_set_compare_other(self, setCompare_file, to_keep):
        with open(setCompare_file, "r") as filin:
            for line in filin:
                regex_nb_hits = re.search("^# ([0-9]+)", line)
                if regex_nb_hits:
                    nb_hits = int(regex_nb_hits.group(1))
                    if nb_hits == 0:
                        return
                    break

            index_dic = OrderedDict()
            i = 0
            for rank_occ in filin:
                if i == to_keep or rank_occ == "\n": break
                rank_splitted = rank_occ.split(":")
                rankw20_occ = rank_splitted[1].split("[")
                index_dic[int(rank_splitted[0])] = [rankw20_occ[0], rankw20_occ[1][:-2].split(",")]
                i += 1
        return index_dic, nb_hits

    def generate_json_data(self):
        """Generate json data for client (to display into table) from Hit collection. 
        
        :return: json data, 
        :rtype: List[Dict]
        """

        results = []
        sorted_hits = sorted(self.hits_collection, key=lambda hit: hit.number_occurences)
        #logging.debug([hit.number_occurences for hit in self.hits_collection])
        logging.debug(max([hit.number_occurences for hit in self.hits_collection]))
        for hit in sorted_hits:
            results.append({"sequence" : hit.sequence, "occurences" : hit.list_occ(self.include_taxon)})
        logging.debug(results[-1])
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
        
    def format_results(self):
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
        final_json["size"] = self.generate_json_size()

        return final_json
       
            





