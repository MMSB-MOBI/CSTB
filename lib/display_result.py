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


class Hit():
    """
    Object Hit containing the CRISPR sequence, indexes where it is found
    and its associated score
    """
    def __init__(self, sequence, score):
        self.sequence = sequence
        self.genomes_dict = {}
        self.score = score

    def set_genomes_dict(self, dic_seq):
        self.genomes_dict = dic_seq

    def write(self, genomes_in):
        to_write = ""
        for gi in genomes_in:
            for ref in self.genomes_dict[gi]:
                to_write += ref + ':' + ','.join(self.genomes_dict[gi][ref]) + ';'
        return to_write.strip(';')

    def list_ref(self, org_name):
        list_ref = [{"ref": ref, "coords": self.genomes_dict[org_name][ref]} for ref in self.genomes_dict[org_name]]
        return list_ref

    def list_occ(self):
        list_occurences = [{'org': genome, 'all_ref': self.list_ref(genome)} for genome in self.genomes_dict]
        return list_occurences


def eprint(*args, **kwargs):
    """
    For printing only on the terminal
    """
    print(*args, file=sys.stderr, **kwargs)


def write_to_file(genomes_in, genomes_not_in, dic_hits, pam, non_pam_motif_length, workdir, nb_top, hit_obj, list_ordered):
    """
    Write results in a file.
    The file is a tabulated file, with first column=sgrna sequence,
    then one column for each included organisms with list of positions
    of the sequence in each.
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


def create_list_ref(dic_ref):
    """
    Create a list of references
    """
    list_ref = []
    for ref in dic_ref:
        dic = {"ref": ref, "coords": dic_ref[ref]}
        list_ref.append(dic)
    return list_ref


def create_list_occurences(dic_occurences):
    """
    Create a list of occurences
    """
    list_occurences = []
    for genome in dic_occurences:
        dic_genome = {'org': genome, 'all_ref': create_list_ref(dic_occurences[genome])}
        list_occurences.append(dic_genome)
    return list_occurences


def output_interface(dic_hits, workdir, nb_top):
    """
    Reformat the results to create a json file.
    It will be parsed in javascript to display it in interface.
    """
    # eprint('\n\n-- Construct results for graphical interface --')
    json_result_file = workdir + '/results.json'
    list_dic = [{'sequence': sgrna, 'occurences':dic_hits[sgrna].list_occ()} for i, sgrna in itertools.takewhile(lambda x: x[0] < nb_top, enumerate(dic_hits)) ]
    list_dic.reverse()
    with open(json_result_file, 'w') as filout:
        json.dump(list_dic, filout, indent=4)


def display_hits(dic_hits, genomes_in, genomes_not_in, pam, non_pam_motif_length, workdir, nb_top, hit_obj, list_ordered):
    """
    write output for interface
    """
    # Put results in local file for access via the interface.
    write_to_file(genomes_in, genomes_not_in, dic_hits, pam,
                  non_pam_motif_length, workdir, nb_top, hit_obj, list_ordered)

    # Output formatting for printing to interface
    output_interface(dic_hits, workdir, 100)

    dic_sorted_by_org = {genome : {} for genome in genomes_in}
    if hit_obj:
        for seq in dic_hits:
            for genome in genomes_in:
                for ref in dic_hits[seq].genomes_dict[genome]:
                    if ref not in dic_sorted_by_org[genome]:
                        dic_sorted_by_org[genome][ref] = {}
                    dic_sorted_by_org[genome][ref][seq] = dic_hits[seq].genomes_dict[genome][ref]
    else:
        for seq in dic_hits:
            for genome in genomes_in:
                for coord_seq in dic_hits[seq].dic_org[genome]:
                    ref = coord_seq.ref
                    if ref not in dic_sorted_by_org[genome]:
                        dic_sorted_by_org[genome][ref] = {}
                    dic_sorted_by_org[genome][ref][seq] = coord_seq.list_coord
    json.dump(dic_sorted_by_org, open("results_by_org.json", "w"))
