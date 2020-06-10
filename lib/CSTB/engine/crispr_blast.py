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
import logging


REGEX_ID = [re.compile("^GCF_[0-9]{9}.[0-9]{1}$"), re.compile("CP[0-9]{6}.[0-9]{1}")]

#class ResumeSeq():


class BlastHit(object):
    """docstring for BlastHit."""
    def __init__(self, org_uuid, fasta_header, hsp):
        self.org_uuid = org_uuid
        self.fasta_header = fasta_header
        hsp_from = int(hsp.find("Hsp_hit-from").text) - 1
        hsp_to = int(hsp.find("Hsp_hit-to").text) - 1
        self.start = min(hsp_from, hsp_to)
        self.end = max(hsp_from, hsp_to)
        self.len = int(hsp.find("Hsp_align-len").text)

    #def __repr__(self):
    #    return "\nStart : {}\nEnd : {}\nLength aln : {}".format(self.start, self.end, self.len)


class BlastReport(object):
    """docstring for BlastHit."""
    def __init__(self, f_name, id_min, genomes_in):
        tree = ET.parse(f_name)
        self.root = tree.getroot()
        self.homolog_genes = self._parse_(id_min, genomes_in) if self.is_hit() else {}

    def __getitem__(self, k):
        return [gene for gene in self.homolog_genes if gene.org_uuid == k]

    def _parse_(self, id_min, genomes_in):
        list_hits = []
        len_query = int(self.root.find("BlastOutput_query-len").text)
        for hit in self.root.iter(tag="Hit"):
            org_uuid = hit.find("Hit_def").text.split("|")[0]
            fasta_header = hit.find("Hit_def").text.split("|")[1].split(" ")[0] #Maybe have a way to identify this properly ?
            if org_uuid in genomes_in:
                #Check hsp and store if %id is superior to id_min
                for hsp in hit.iter(tag="Hsp"):
                    if (int(hsp.find("Hsp_identity").text)/ len_query) * 100 > id_min:
                        # Create a list of BlastHit object for the reference
                        list_hits.append(BlastHit(org_uuid, fasta_header, hsp))
        
        return list_hits

    def is_hit(self):
        """
        Return true if there are hits else return False
        """
        msg = self.root.find("./BlastOutput_iterations/Iteration/Iteration_message")
        if not hasattr(msg, "text"):
            return True
        msg = msg.text
        return False
    
    @property
    def organisms(self):
        """
        Return the list of organism name which have homologous genes
        """
        return set([gene.org_uuid for gene in self.homolog_genes])

    def filterGenes(self, orgs, fasta_header = None):
        """ Filter genes from blast report. Just keep genes corresponding to a list of organisms and to fasta header if given.
        
        :param orgs: List of organisms for which we want to return gene
        :type orgs: List[str]
        :param fasta_header: List of fasta header for which we want to return gene, defaults to None
        :type fasta_header: str, optional
        :return: list of blast hits corresponding to kept genes
        :rtype: List[BlastHit]
        """
        filtered_genes = [gene for gene in self.homolog_genes if gene.org_uuid in orgs]
        filtered_genes = sorted(filtered_genes, key=lambda gene:gene.start) #genes in apparition on genome order.
        if fasta_header:
            return [gene for gene in filtered_genes if gene.fasta_header == fasta_header]
        
        return filtered_genes