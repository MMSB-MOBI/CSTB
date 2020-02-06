#!/usr/bin/env python3
"""
Display results:
Contain object Hit, which contains the sgRNA, its score and the dictionary containing
genomes and references where it is present and its coordinates

Write json and text files with results from the intersection of genomes

****** Format of the json file ******
{'sequence' : word, 'occurences' :
                                    {'org' : genome, 'all_ref' :
                                                                {'ref' : ref, 'coords' : [coordinates]
                                                               }
                                    }
}

****** Format of the text file ******
#ALL GENOMES
#Genomes included : *list of genomes*  ; Genomes excluded : *list of genomes*
#Parameters: pam: *PAM* ; sgrna size: *size*
        genome_in_1     genome_in_2     genomes_in_3...
word    ref : coordinates   ref : coordinates   ref : coordinates...
.
.
.

"""

import sys
import json
import itertools
import CSTB.engine.wordIntegerIndexing as decoding
import logging
import CSTB.utils.error as error



class Hit():
    """Hit object
    
    :ivar index: setCompare index
    :vartype index: int
    :ivar weight: setCompare weigth, correspond to number of sgRNA occurences in all genomes
    :vartype weight: int
    :ivar sequence: sgRNA sequence in nucleotides
    :vartype sequence: str
    :ivar occurences: Occurences of sgRNA in genomes
    :vartype occurences: Dict { genome : { fasta_header : List[str] } }
   
    """
    def __init__(self, index, weight):
        """ Initialize an Hit object
        
        :param index: setCompare index
        :type index: int
        :param weight: setCompare weight, correspond to number of sgRNA occurences in all genomes
        :type weight: int
        """
        self.index = int(index)
        self.weight = weight
        self.sequence = self.decode()
        self.occurences = {} #This dictionnary will be {"organism": {"subsequence" : coords[]}}

    def _list_ref(self, org_name):
        """Format occurences for an organism 
        """
        list_ref = [{"ref": ref, "coords": self.occurences[org_name][ref]} for ref in self.occurences[org_name]]
        return list_ref

    def list_occ(self, taxon_name):
        """Format occurences for an organism
        """
        list_occurences = [{'org': taxon_name[genome], 'all_ref': self._list_ref(genome)} for genome in self.occurences]
        return list_occurences

    @property
    def number_occurences(self):
        """Give the total number of hit occurences
        
        :return: Total number of hit occurences
        :rtype: int
        """
        nb_occ = 0
        for genome in self.occurences:
            for subseq in self.occurences[genome]:
                nb_occ += len(self.occurences[genome][subseq])
        return nb_occ

    def __str__(self):
        return f"index:{self.index}\nweight:{self.weight}\nsequence:{self.sequence}\noccurences:{self.occurences}"

    def store_occurences(self, couch_doc, genomes_in):
        """Set occurences attributes from couch document
        
        :param couch_doc: couch response
        :type couch_doc: Dict
        :param genomes_in: list of genomes uuid
        :type genomes_in: List[str]
        :raises error.ConsistencyError: Raise if at least 1 given genome is not in couch document. 
        """
        db_genomes = set(couch_doc.keys())
        logging.debug(f"store_occurences\ncouch_doc: {couch_doc}\ngenomes_in: {genomes_in}")
        if not set(genomes_in).issubset(db_genomes):
            raise error.ConsistencyError(f"Consistency error. Genomes included {genomes_in} are not in couch database for {self.sequence}")
        
        self.occurences = {genome:couch_doc[genome] for genome in genomes_in}

    def decode(self):
        """Decode setCompare index into nucleotide sequence
        
        :return: sgRNA sequence
        :rtype: str
        """
        return decoding.decode(self.index, ["A", "T", "C", "G"], 23)


def write_to_file(genomes_in, genomes_not_in, dic_hits, pam, non_pam_motif_length, workdir, nb_top, hit_obj, list_ordered):
    """
    Write results in a file.
    The file is a tabulated file, with first column=sgrna sequence,
    then one column for each included organisms with list of positions
    of the sequence in each.

    OBSOLETE, maybe we need this if we keep the option to download results.
    """
    eprint('\n\n-- Write results to file --')
    rep_rslt_file = workdir + '/results_allgenome.txt'
    output = open(rep_rslt_file, 'w')
    gen_i = ','.join(genomes_in)
    gen_ni = ','.join(genomes_not_in)
    if gen_ni == '':
        gen_ni = 'None'
    output.write('#ALL GENOMES\n#Genomes included :' + gen_i +
                 ' ; Genomes excluded :' + gen_ni + '\n'+'#Parameters: pam:' +
                 pam + ' ; sgrna size:' + str(non_pam_motif_length) + '\n')
    output.write('sgrna sequence')
    for genome_i in genomes_in:
        output.write('\t' + genome_i)
    output.write('\n')
    i = 0
    for sgrna in list_ordered:
        if i == nb_top : break
        i += 1
        hit = dic_hits[sgrna]
        output.write(sgrna+'\t')
        to_write = hit.write(genomes_in)
        output.write(to_write + '\n')
    output.close()