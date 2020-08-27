'''
create_metafile.py at the origin
'''
from CSTB_core.engine.word_detect import sgRNAfastaSearch
from CSTB_core.engine.wordIntegerIndexing import indexAndMayOccurence, setEncoding
from CSTB_core.utils.io import sgRNAIndexWriter
import os.path
import argparse
import logging
logging.basicConfig(filename = "index_sequence.log", level = logging.DEBUG, format='%(levelname)s\t%(message)s')
import CSTB.utils.error as error


def args_gestion():
    parser = argparse.ArgumentParser(description="Search sgrna in nucleotide sequence and index it for setCompare")
    parser.add_argument("-f", "--fasta", metavar = "<file>", help = "Fasta file with sequence", required = True)
    parser.add_argument("-c", "--codec", metavar = "<str>", help = "Encoding mode (twobits or pow2)", default = "twobits")
    parser.add_argument("-o", "--output", metavar = "<file>", help = "Results file with indexes", required = True)

    args = parser.parse_args()

    if args.codec not in ["twobits", "pow2"]:
        raise error.ArgumentError("-c/--codec has to be twobits or pow2")

    return args

if __name__ == "__main__":
    logging.info("== index_sequence.py")
    ARGS = args_gestion()
    if not (os.path.isfile(ARGS.fasta)):
        raise FileNotFoundError(f"{ARGS.fasta} doesn't exist")
    sgRNA_seqs = sgRNAfastaSearch(ARGS.fasta, 'gene')
    if not sgRNA_seqs:
        error.empty_exit("No sgRNA found in gene")
    logging.info(sgRNA_seqs)
    setEncoding(ARGS.codec)
    indexData, word_length = indexAndMayOccurence(sgRNA_seqs)
    sgRNAIndexWriter(indexData, ARGS.output, word_length, ARGS.codec)
    #writeIndexes(indexData, ARGS.output)
    
    

    
