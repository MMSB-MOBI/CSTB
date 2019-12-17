#!/usr/bin/env python3
"""
Parse the output of a Blastn on the reference database. Create a BlastReport object which
contains only hits for organism which are in a given list, here the included genomes list, and
hits which have a percentage identity superior to a given percentage identity, by default 70.
The BlastReport contains a dictionary of organism with references which contains BlastHit object
This object contains the coordinates of the hit and its length.
"""

import argparse
import pickle
import re
import xml.etree.ElementTree as ET
import re

REGEX_ID = [re.compile("^GCF_[0-9]{9}.[0-9]{1}$"), re.compile("CP[0-9]{6}.[0-9]{1}")]

class BlastHit(object):
    """docstring for BlastHit."""
    def __init__(self, hsp):
        hsp_from = int(hsp.find("Hsp_hit-from").text) - 1
        hsp_to = int(hsp.find("Hsp_hit-to").text) - 1
        self.start = min(hsp_from, hsp_to)
        self.end = max(hsp_from, hsp_to)
        self.len = int(hsp.find("Hsp_align-len").text)

    def __repr__(self):
        return "\nStart : {}\nEnd : {}\nLength aln : {}".format(self.start, self.end, self.len)


class BlastReport(object):
    """docstring for BlastHit."""
    def __init__(self, f_name, id_min, genomes_in):
        tree = ET.parse(f_name)
        self.root = tree.getroot()
        self.homolog_gene = self._parse_(id_min, genomes_in) if self.is_hit() else {}

    def __getitem__(self, k):
        if k in self.homolog_gene:
            return self.homolog_gene[k]
        return None

    def _parse_(self, id_min, genomes_in):
        data = {}
        len_query = int(self.root.find("BlastOutput_query-len").text)
        for hit in self.root.iter(tag="Hit"):
            org_name = hit.find("Hit_def").text
            ref = org_name.split(" ")[0]
            subdata = {ref : []}
            # Find the name given in the database and if it is in the genome_in list
            real_name = "".join(list(filter(lambda x: check_in_gi(x, org_name), genomes_in)))
            if real_name:
                for hsp in hit.iter(tag="Hsp"):
                    # Check if the identity percentage is superior to the identity given by user
                    if (int(hsp.find("Hsp_identity").text)/ len_query) * 100 > id_min:
                        # Create a list of BlastHit object for the reference
                        subdata[ref].append(BlastHit(hsp))
                # Check if the organism has already been seen because an organism can have
                # several references so need to update not to assign
                if real_name not in data:
                    data[real_name] = {}
                data[real_name].update(subdata)
        return data

    def is_hit(self):
        """
        Return true if there are hits else return False
        """
        msg = self.root.find("./BlastOutput_iterations/Iteration/Iteration_message")
        if not hasattr(msg, "text"):
            return True
        msg = msg.text
        print(msg)
        return False

    def org_names(self):
        """
        Return the list of organism name which have homologous genes
        """
        return list(self.homolog_gene.keys())

    def ref_names(self, org_name):
        """
        Return reference which have homologous gene for a given organism
        which have homologous gene
        """
        if org_name in self.org_names():
            return list(self[org_name].keys())
        return None

    def gene_coords(self, org_name, ref_name):
        """
        Return a list of object BlastHit for a reference and organism given
        """
        if ref_name in self.ref_names(org_name):
            return list(self[org_name][ref_name])
        return None

    def json_str(self):
        """
        Return a string in JSON format
        """
        tmp = {}
        for org in self.homolog_gene:
            tmp[org] = {}
            for ref in self.homolog_gene[org]:
                tmp[org][ref] = []
                for blastObj in self.homolog_gene[org][ref]:
                    tmp[org][ref].append({"start": str(blastObj.start), "end":str(blastObj.end)})
        return tmp


def check_in_gi(gi_name, blast_name):
    """
    Return a match or none if the name of the blast hit is in the list
    of included genomes
    """
    # Remove the GCF code of the end if exists
    putative_gcf = gi_name.split(" ")[-1]
    check_id = [regex.match(putative_gcf) for regex in REGEX_ID if regex.match(putative_gcf)]
    if check_id:
        gi_name = " ".join(gi_name.split(" ")[:-1]).strip()
    return re.search(gi_name + "[ ,$]", blast_name)


def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    parser = argparse.ArgumentParser(description="Parse the output of blastn")
    parser.add_argument("-blast", metavar="<str>",
                        help="The path to the blastn output file",
                        required=True)
    parser.add_argument('-ip', type=int,
                        help='identity percentage min for the research of\
                              homologous genes using blastn (default:70)',
                        nargs="?", default=70)
    parser.add_argument("-gi", metavar="<str>",
                        help="The organisms to search inclusion in.",
                        required=True)
    parser.add_argument("-o", metavar="<str>",
                        help="The output path for the pickle file",
                        nargs="?", default="output_blast.p")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    PARAM = args_gestion()
    GENOMES_IN = PARAM.gi.split("&")
    OUTPUT_BLAST = BlastReport(PARAM.blast, PARAM.ip, GENOMES_IN)
    if OUTPUT_BLAST.is_hit():
        pickle.dump(OUTPUT_BLAST, open(PARAM.o, "wb"), protocol=3)
    else:
        print("Program terminated&No homologous gene found")
