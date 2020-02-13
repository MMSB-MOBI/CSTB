'''
create_metafile.py at the origin
'''
from CSTB_core.engine.word_detect import sgRNAfastaSearch
from CSTB_core.engine.wordIntegerIndexing import indexAndOccurence, writeIndexes
import os.path
import argparse
import logging
logging.basicConfig(filename = "index_sequence.log", level = logging.DEBUG, format='%(levelname)s\t%(message)s')
from CSTB.utils.error import empty_exit


def args_gestion():
    parser = argparse.ArgumentParser(description="Search sgrna in nucleotide sequence and index it for setCompare")
    parser.add_argument("-f", "--fasta", metavar = "<file>", help = "Fasta file with sequence", required = True)
    parser.add_argument("-o", "--output", metavar = "<file>", help = "Results file with indexes", required = True)
    return parser.parse_args()

if __name__ == "__main__":
    logging.info("== index_sequence.py")
    ARGS = args_gestion()
    if not (os.path.isfile(ARGS.fasta)):
        raise FileNotFoundError(f"{ARGS.fasta} doesn't exist")
    sgRNA_seqs = sgRNAfastaSearch(ARGS.fasta, 'gene')
    if not sgRNA_seqs:
        empty_exit("No sgRNA found in gene")
    indexData = indexAndOccurence(sgRNA_seqs)
    writeIndexes(indexData, ARGS.output)
    
    

    
