'''
create_metafile.py at the origin
'''
from CSTB.engine.word_detect import sgRNAfastaSearch
from CSTB.engine.wordIntegerIndexing import indexAndOccurence, writeIndexes
import os.path
import argparse

def args_gestion():
    parser = argparse.ArgumentParser(description="Search sgrna in nucleotide sequence and index it for setCompare")
    parser.add_argument("-f", "--fasta", metavar = "<file>", help = "Fasta file with sequence", required = True)
    parser.add_argument("-o", "--output", metavar = "<file>", help = "Results file with indexes", required = True)
    return parser.parse_args()

if __name__ == "__main__":
    ARGS = args_gestion()
    if not (os.path.isfile(ARGS.fasta)):
        raise FileNotFoundError(f"{ARGS.fasta} doesn't exist")
    sgRNA_seqs = sgRNAfastaSearch(ARGS.fasta, 'gene')
    indexData = indexAndOccurence(sgRNA_seqs)
    writeIndexes(indexData, ARGS.output)
    
    

    
