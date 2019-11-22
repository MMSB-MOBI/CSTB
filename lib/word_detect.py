#!/usr/bin/env python3
"""
With regex, find sgRNA
"""

import os
import re
import argparse
import pickle
from Bio import SeqIO


def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    command = sys.argv[1] if len(sys.argv) > 1 else None
    # Argparsing
    parser = argparse.ArgumentParser(description="Find sgRNA sequences")
    parser.add_argument("-file", metavar="<str>",
                        help="The fasta file to parse",
                        required=True)
    parser.add_argument("-org", metavar="<str>",
                        help="Name of the organism",
                        required=True)
    parser.add_argument("-out", metaver="<str>",
                        help="The output file",
                        required=True)
    parser.add_argument("-pam", metaver="<str>",
                        help="The PAM motif, default : NGG",
                        nargs='?', const="NGG")
    parser.add_argument("-sl", metaver="int",
                        help="The length of the motif without PAM, default : 20",
                        nargs='?', const=20)
    args = parser.parse_args()
    return args


def complement_seq(sequence):
    """
    Function for turning a 5'-3' nucleotidic sequence into
    its 5'-3' reverse complement.
    """
    rev_comp = []
    complement = {"A" : "T", "C" : "G", "T" : "A", "G" : "C"}
    rev_comp = [complement[sequence[idx]] if sequence[idx] in complement.keys() else "N" for idx in range(len(sequence) - 1, -1, -1)]
    return "".join(rev_comp)


def build_expression(seq):
    """
    Build the regular expression by replacing the letter by iupac code or
    the letter itself
    """
    result = ''
    iupac_code = {'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
                  'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
                  'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'}
    for letter in seq:
        if letter in iupac_code:
            result = result + iupac_code[letter]
        else:
            result = result + letter
    return result


def find_indices_sgrna(seq, pam):
    """
    Uses Regular expression matching of the pam motif to the reference genome
    to get the start positions (0-based) of each match
    """
    reg_exp = build_expression(pam)
    indices = [m.start() for m in re.finditer('(?=' + reg_exp + ')', seq, re.I)]
    return indices


def find_sgrna_seq(seq_list, len_seq, reverse, str_reverse, seq_dict, genome_seq, organism, ref):
    """
    Uses start index and the length of the sequence to get the sequence and
    fill the dictionary containing all sgRNA sequences with coordinates
    """
    for indice in seq_list:
        end = indice + len_seq
        seq = genome_seq[indice:end] if reverse else genome_seq[indice:end].reverse_complement()
        seq = str(seq)
        if seq not in seq_dict:
            seq_dict[seq] = {organism: {}}
        if ref not in seq_dict[seq][organism]:
            seq_dict[seq][organism][ref] = []
        seq_dict[seq][organism][ref].append(str_reverse + str(indice+1) + ',' +
                                            str(end) + ')')
    return seq_dict


def construct_in(fasta_file, pickle_file, organism, pam="NGG", non_pam_motif_length=20):
    """
    Construct the sequences for first organism,
    with python regular expression research
    """
    sgrna = "N" * non_pam_motif_length + pam
    seq_dict = {}

    for genome_seqrecord in SeqIO.parse(fasta_file, "fasta"):
        genome_seq = genome_seqrecord.seq
        ref = genome_seqrecord.id
        seq_list_forward = find_indices_sgrna(str(genome_seq),
                                              complement_seq(sgrna))
        seq_list_reverse = find_indices_sgrna(str(genome_seq), sgrna)

        seq_dict = find_sgrna_seq(seq_list_forward, len(pam) + non_pam_motif_length,
                                  False, "+(", seq_dict, genome_seq, organism, ref)
        seq_dict = find_sgrna_seq(seq_list_reverse, len(pam) + non_pam_motif_length,
                                  True, "-(", seq_dict, genome_seq, organism, ref)

    pickle.dump(seq_dict, open(pickle_file, "wb"), protocol=3)
    return seq_dict


if __name__ == '__main__':
    PARAM = args_gestion()
    construct_in(PARAM.file, PARAM.out, PARAM.org,PARAM.pam, PARAM.sl)

    # filin = "../../test/reference_genomes/fasta/GCF_001022195.1_ASM102219v1/GCF_001022195.1_ASM102219v1_genomic.fna"
    # filout = "../../test/reference_genomes/test"
    # PAM = "NGG"
    # SL = 20
    # construct_in(filin, filout, pam="NGG", non_pam_motif_length=20)
