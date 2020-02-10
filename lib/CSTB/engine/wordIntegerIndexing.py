"""
    Translate a dictionary of constant-length CRISPR motifs into an ordered set of integer, using base4 encoding.
    Usage:
        wordIntegerIndexing.py <pickledDictionary> [--out=<outFile> --occ]

    Options:
        -h --help                               Show this screen
        -o <outFile>, --out <outFile>           Name of the output file [default : ./pickledDictionary.index]

"""

import os
import pickle
import math
from docopt import docopt

def occWeight(data):
    k,datum = data
    n = 0
    for o in datum:
        for _o in datum[o]:
            n += len(datum[o][_o])
    #print(k,n)
    return n

# same as index pickle, coding and order wise, 
# we just add a second field to each wordCode line, the occurence number
def indexAndOccurencePickle(file_path, target_file):

    p_data = pickle.load(open(file_path, "rb"))
    word_list = list(p_data.keys())
    wordLen = len(word_list[0])
    data = sorted( [ (weightWord(w, ["A", "T", "C", "G"], wordLen), occWeight((w, p_data[w])) ) for w in word_list], key=lambda x: x[0])

    with open(target_file, "w") as filout:
        filout.write(str(len(data)) + "\n")
        for datum in data:
            filout.write( ' '.join([str(d) for d in datum]) + "\n")

    return len(data)

# same as index pickle, coding and order wise, 
# we just add a second field to each wordCode line, the occurence number
def indexAndOccurence(data):
    word_list = list(data.keys())
    wordLen = len(word_list[0])
    data = sorted( [ (weightWord(w, ["A", "T", "C", "G"], wordLen), occWeight((w, data[w])) ) for w in word_list], key=lambda x: x[0])
    return data
    
def indexPickle(file_path, target_file):
    """
    Take a pickle file, code it and write it in a file
    """
    p_data = pickle.load(open(file_path, "rb"))
    word_list = list(p_data.keys())
    data = sorted([weightWord(w, ["A", "T", "C", "G"], len(word_list[0])) for w in word_list])

    with open(target_file, "w") as filout:
        filout.write(str(len(data)) + "\n")
        for coding_int in data:
            filout.write(str(coding_int) + "\n")

    return len(data)


def weightWord(word, alphabet, length=None):
    """
    Code a word by base len(alphabet) and return this int
    """
    rank = 0
    if length:
        if length != len(word):
            raise ValueError("Irregular word length " + str(len(word)) +
                             " (expected " + str(length) + ")")
    for i, letter in enumerate(word[::-1]):
        wei = alphabet.index(letter)
        base = len(alphabet)
        rank += wei * pow(base, i)
    return rank

def project(value, lenFrom, lenTo, alphabet="ATCG"):
    base = len(alphabet)
    _value = value
    offset = 0
    for i in range(lenFrom, lenTo-1, -1):
        w = math.trunc(value / pow(base, i))
        offset += w * pow(base, i)
        value = value % pow(base, i)
    
    return _value - offset

def decode(rank, alphabet, length=20):
    """
    Decode the rank (int) to a word according to an alphabet given
    """
    word = ""
    base = len(alphabet)
    for i in range(length - 1, -1, -1):
        index = math.trunc(rank / pow(base, i))
        word += alphabet[index]
        rank = rank % pow(base, i)
    return word

def writeIndexes(indexData, output):
    with open(output, "w") as o:
        o.write(f"{len(indexData)}\n")
        for datum in indexData:
            o.write(f"{' '.join([str(d) for d in datum])}\n")

if __name__ == "__main__":
    ARGUMENTS = docopt(__doc__, version='wordIntegerIndexing 1.0')

    TARGET_FILE = ('.'.join(os.path.basename(ARGUMENTS['<pickledDictionary>']).split('.')[0:-1])
                   + '.index')
    if ARGUMENTS['--out']:
        TARGET_FILE = ARGUMENTS['--out']
    indexFn = indexPickle
    if ARGUMENTS['--occ']:
        indexFn = indexAndOccurencePickle
    TOTAL = indexFn(ARGUMENTS['<pickledDictionary>'], TARGET_FILE)
    print("Successfully indexed", TOTAL, "words\nfrom:",
          ARGUMENTS['<pickledDictionary>'], "\ninto:", TARGET_FILE)
