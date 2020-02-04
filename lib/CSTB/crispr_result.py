import logging
import pycouch.error
import re
from collections import OrderedDict
import CSTB.display_result as dspl
import CSTB.error as error
import requests
import time

SESSION = requests.Session()
SESSION.trust_env = False

class CrisprResult:
    def __init__(self, include, exclude):
        self.genomes_include = include
        self.genomes_exclude = exclude

    def serialize_genomes_info(self):
        pass


class CrisprResultManagerOld:
    def __init__(self, pycouch_wrapper, taxondb, genomedb, motif_broker_endpoint, tag):
        self.wrapper = pycouch_wrapper
        self.motif_broker_endpoint = motif_broker_endpoint
        self.taxondb = taxondb
        self.genomedb = genomedb
        self.tag = tag
        self.nb_total_hits = None
        self.hits_collection = []
        self.include_taxon = {}
        self.exclude_taxon = {}

    #Move this in some consumer
    def get_taxon_name(self, list_uuid):
        taxon_names = {}
        for g_uuid in list_uuid:
            resp = self.wrapper.couchGetDoc(self.genomedb, g_uuid)
            if not resp:
                raise error.CouchNotFound(f"{self.wrapper.endpoint}/{self.genomedb}/{g_uuid} not found")
            taxon_uuid = resp["taxon"]

            resp = self.wrapper.couchGetDoc(self.taxondb, taxon_uuid)
            if not resp:
                raise error.CouchNotFound(f"{self.wrapper.endpoint}/{self.taxondb}/{taxon_uuid} not found")
            taxon_names[g_uuid] = resp["name"]

        return taxon_names

    #Move this in some consumer
    def get_genomes_size(self, list_uuid):
        sizes = {}
        for g_uuid in list_uuid:
            resp = self.wrapper.couchGetDoc(self.genomedb, g_uuid)
            if not resp:
                raise error.CouchNotFound(f"{self.wrapper.endpoint}/{self.genomedb}/{g_uuid} not found")
            size = resp["size"]
            sizes[g_uuid] = size
        return sizes
            
    def set_taxon_names(self, include, exclude):
        self.include_taxon = self.get_taxon_name(include)
        self.exclude_taxon = self.get_taxon_name(exclude)

    def search_occurences(self, genomes_include, len_slice = 2000):
        try:
            SESSION.get(self.motif_broker_endpoint + "/handshake")
        except:
            error.error_exit("Can't handshake motif-broker")

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
                        error.error_exit("Can't access to database after 3 tries")
                    time.sleep(5)
                    continue
                break
        
        for hit in self.hits_collection:
            hit.store_occurences(results[hit.sequence], genomes_include)
           
    def compute(self, include, exclude, setCompare_f, word_length, to_keep = 5000):
        include_name = self.get_taxon_name(include)
        exclude_name = self.get_taxon_name(exclude)

        self.stored_taxon_name = {**include_name, **exclude_name}

        try:
            self.parse_set_compare(setCompare_f, word_length, to_keep)
        except : 
            error.error_exit("Error while parsing set compare")

        for hit in self.hits_collection:
            logging.debug(f"Hit object:\n{hit}")
            break
        
        logging.info(f"{self.nb_total_hits} hits found. {to_keep} first are kept.")

        try:
            self.search_occurences(include)
        except:
            error.error_exit("Error while searching occurences in genomes.")
        
        for hit in self.hits_collection:
            logging.debug(f"Hit object:\n{hit}")
            break

        results_json = self.format_results()

        json = {"gi" : "&".join(list(include_name.values())), "not_in" : ",".join(list(exclude_name.values())), "data": results_json}

        return json

    def parse_set_compare(self, setCompare_file, word_length, to_keep):
        self.hits_collection=[]
        logging.info(f"Parse set compare for word length {word_length}")
        if word_length == 20:
            self._parse_set_compare_20(setCompare_file, to_keep)
        else:
            logging.debug("Handle this")
            
    def _parse_set_compare_20(self, setCompare_file, to_keep):
        with open(setCompare_file, "r") as filin:
            text = filin.readlines()

        self.nb_total_hits = re.search("[0-9]+", text[-2]).group()
        if int(self.nb_total_hits) == 0:
            error.error_exit("No hits")

        index_dic = OrderedDict()
        i = 0
        for rank_occ in text[-1].strip().split(","):
            if i == to_keep: break
            self.hits_collection.append(dspl.Hit(rank_occ.split(":")[0],rank_occ.split(":")[1]))
            #index_dic[int(rank_occ.split(":")[0])] = rank_occ.split(":")[1]
            i += 1

    def _parse_set_compare_other(self, setCompare_file, to_keep):
        with open(setCompare_file, "r") as filin:
            for line in filin:
                regex_nb_hits = re.search("^# ([0-9]+)", line)
                if regex_nb_hits:
                    nb_hits = int(regex_nb_hits.group(1))
                    if nb_hits == 0:
                        error.error_exit("No hits")
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
        results = []
        for hit in self.hits_collection:
            results.append({"sequence" : hit.sequence, "occurences" : hit.list_occ(self.include_taxon)})
        return results

    def generate_json_data_card(self):

        def insert_in_data_card(genome, subseq, sequence, coords):
            if genome not in data_card:
                data_card[genome] = {}
            
            if subseq not in data_card[genome]:
                data_card[genome][subseq] = {}
            
            data_card[genome][subseq][hit.sequence] = hit.occurences[genome][subseq]

        data_card = {}
        for hit in self.hits_collection:
            for genome in hit.occurences:
                for subseq in hit.occurences[genome]:
                    insert_in_data_card(genome, subseq, hit.sequence, hit.occurences[genome][subseq])
        
        return data_card

    def generate_json_size(self):
        json_size = {}
        sizes = self.get_genomes_size(list(self.include_taxon.keys()))
        for g_uuid in sizes:
            json_size[self.include_taxon[g_uuid]] = sizes[g_uuid]  
        
        return json_size
        
    def format_results(self):
        final_json = {}
        final_json["gi"] = "&".join(list(self.include_taxon.values()))
        final_json["not_in"] = ",".join(list(self.exclude_taxon.values()))
        final_json["number_hits"] = self.nb_total_hits
        final_json["data"] = self.generate_json_data()
        final_json["data_card"] = self.generate_json_data_card()
        final_json["tag"] = self.tag
        final_json["size"] = self.generate_json_size()

        return final_json
       
            





