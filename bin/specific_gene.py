#!/usr/bin/env python3
"""
Retrieve sgRNA on a specific gene given by user. Retrieve result from SetCompare script and
create files for display in webservice.
"""

import sys
import argparse
import pickle
import re
import json
from collections import OrderedDict
import wordIntegerIndexing as decoding
import display_result as dspl
import post_processing as pp
from parse_blast import BlastReport, BlastHit


class ResumeSeq(object):
    """docstring for ResumeSeq."""
    def __init__(self):
        self.proportion = 0
        self.dic_org = {}

    def __repr__(self):
        return "Proportion : {}\tTotal occurences : {}".format(self.proportion, self.total_occ())

    def __getitem__(self, k):
        if k in self.dic_org:
            return self.dic_org[k]
        return None

    def cal_proportion(self, nb_total):
        """
        Calcul proportion of organism present in the dictionary
        """
        self.proportion = len(self.dic_org) / int(nb_total)
        return self.proportion

    def total_occ(self, org_name=None):
        """
        Calcul the total number of occurences for the sequence or for a given organism
        for a sequence
        """
        occ = 0
        if org_name:
            for coord_seq in self.dic_org[org_name]:
                occ += coord_seq.total_occ()
        else:
            for name in self.dic_org:
                occ += self.total_occ(name)
        return occ

    def update(self, dict_org):
        """
        Update the dictionary
        """
        if dict_org:
            self.dic_org.update(dict_org)

    def org_names(self):
        """
        Return a list of organism name which have this sequence
        """
        return list(self.dic_org.keys())

    def write(self, genomes_in):
        """
        Return a string of coordinates for each reference for each organism
        given an order by the list
        """
        to_write = ""
        for gen_in in genomes_in:
            if gen_in in self.org_names():
                for coord_seq in self[gen_in]:
                    to_write += coord_seq.ref + ':' + ','.join(coord_seq.list_coord) + ';'
            else:
                to_write += "None ;"
        to_write += str(self.total_occ())
        return to_write.strip(";")

    def list_ref(self, org_name):
        """
        Return a list of dictionaries containing coordinates of the sequence for
        a reference organism
        """
        list_ref = [{"ref": coord_seq.ref, "coords": coord_seq.list_coord} for coord_seq in self.dic_org[org_name]]
        return list_ref

    def list_occ(self):
        """
        Return a list of dicionaries containing a list_ref for all references for
        each organism
        """
        list_occurences = [{'org': genome, 'all_ref': self.list_ref(genome)} for genome in self.org_names()]
        return list_occurences



class CoordSeq(object):
    """docstring for CoordSeq."""
    def __init__(self, ref):
        self.list_coord = []
        self.ref = ref

    def __repr__(self):
        return "Total occurences : {}\tFor {} coordinates".format(self.total_occ(), len(self))

    def __len__(self):
        return len(self.list_coord)

    def update(self, coord):
        """
        Update the coordinate list
        """
        if coord not in self.list_coord:
            self.list_coord.append(coord)
            return True
        return False

    def total_occ(self):
        """
        The number of coordinates for the reference
        """
        return len(self.list_coord)


### METHODS FOR THE ARGUMENTS ###
def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    parser = argparse.ArgumentParser(__file__, description="Specific gene program.")
    parser.add_argument("-gi", metavar="<str>",
                        help="The organisms to search inclusion in.",
                        required=True)
    parser.add_argument("-gni", metavar="<str>",
                        help="The organisms to search exclusion from",
                        required=True)
    parser.add_argument("-sl", metavar="<int>",
                        help="The length of the sgrna, excluding pam",
                        nargs="?", default=20)
    parser.add_argument("-pam", metavar="<str>",
                        help="The pam motif",
                        nargs="?", default="NGG")
    parser.add_argument("-f", metavar="<str>",
                        help="The file index",
                        required=True)
    parser.add_argument("-blast", metavar="<str>",
                        help="The blast output file",
                        required=True)
    parser.add_argument("-nb_top", metavar="<int>",
                        help="The top hits to download",
                        default=1000)
    parser.add_argument('--no-proxy', action='store_true')
    parser.add_argument("-r", metavar="<str>",
                        help="The end point",
                        required=True)
    parser.add_argument("-c", metavar="<int>",
                        help="The length of the slice for the request",
                        required=True)
    parser.add_argument("-taxon_db", metavar="<str>",
                        help="The name of the taxon database",
                        required=True)
    parser.add_argument("-tree_db", metavar="<str>",
                        help="The name of the taxon database",
                        required=True)
    parser.add_argument("-end_point", metavar="<str>",
                        help="The end point of the taxon and tree database",
                        required=True)

    return parser.parse_args()


def on_gene(sgrna, gene):
    """
    Check if a given sgrna is on a given gene
    """
    try:
        sgrna_start = int(re.search("[+-]\(([0-9]*),", sgrna).group(1))
        sgrna_end = int(re.search(",([0-9]*)", sgrna).group(1))
        return gene.start <= sgrna_start and gene.end >= sgrna_end
    except:
        return []


def coord_on_gene(list_coord, list_gene, ref):
    """
    For a reference organism check if all sgrna are on genes
    Return a CoordSeq object if at least one sgrna is on a gene
    """
    coord_seq2 = CoordSeq(ref)
    for coord in list_coord:
        nb_on_gene = list(filter(lambda x: on_gene(coord, x), list_gene))
        if nb_on_gene:
            coord_seq2.update(coord)
    if coord_seq2.list_coord:
        return coord_seq2
    return None


def check_on_gene(blastoutput, dic_index, nb_gi):
    """
    Return a dictionary of ResumeSeq object for each sequence which is at least
    on one gene in one organism.
    """
    resume_seq = {}
    # For each sgrna
    for seq in dic_index:
        # Create a ResumeSeq to save info
        resume_seq[seq] = ResumeSeq()
        # For each organism find where an homologous gene has been found
        for org in blastoutput.org_names():
            tmp = {org: []}
            for ref in blastoutput.ref_names(org):
                if org in dic_index[seq].genomes_dict and ref in dic_index[seq].genomes_dict[org]:
                    # Find the list of coordinates on gene
                    coord_seq = coord_on_gene(dic_index[seq].genomes_dict[org][ref],
                                              blastoutput.gene_coords(org, ref), ref)
                else:
                    coord_seq = None
                if coord_seq:
                    tmp[org].append(coord_seq)
            if tmp[org] != []:
                resume_seq[seq].update(tmp)
        # If no coordinate in the ResumeSeq object, then delete the sgrna from the dictionary
        if not resume_seq[seq].dic_org:
            del resume_seq[seq]
        else:
            resume_seq[seq].cal_proportion(nb_gi)
    # If resume_seq is not empty, then at least one sgrna is on one gene
    if resume_seq:
        return resume_seq
    # None sgrna on homologous gene
    print("Program terminated&No hits on genes")
    sys.exit()


if __name__ == "__main__":
    PARAM = args_gestion()
    GENOMES_IN = PARAM.gi.split("&")
    GENOMES_NOTIN = PARAM.gni.split("&")
    DIC_HITS = pp.create_dic_hits(PARAM, GENOMES_IN)

    # Keep sgrna which are on gene
    BLASTOUTPUT = pickle.load(open(PARAM.blast, "rb"))
    RESUME_SEQ = check_on_gene(BLASTOUTPUT, DIC_HITS, len(GENOMES_IN))
    # Sort sgrna by the proportion of organism containing this sgrna on a gene
    LIST_ORDERED = sorted(RESUME_SEQ, key=lambda hit: RESUME_SEQ[hit].proportion, reverse=True)

    # Display the result for the navigator
    dspl.display_hits(RESUME_SEQ, GENOMES_IN, GENOMES_NOTIN,
                      PARAM.pam, int(PARAM.sl), ".", int(PARAM.nb_top), False, LIST_ORDERED)

    pp.get_size(PARAM.end_point, PARAM.tree_db, PARAM.taxon_db, GENOMES_IN)

    print(','.join(GENOMES_NOTIN))
    print("TASK_KEY")
    print(len(RESUME_SEQ))
    json.dump(BLASTOUTPUT.json_str(), open("genes.json", "w"))
