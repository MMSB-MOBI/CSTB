"""

"""

import sys
import re
import urllib
from Bio import SeqIO
from Bio import Entrez
from ete3 import NCBITaxa
import display_result as dspl
DEBUG = True
# CHECK IF THE TAXON ID IS ALREADY PRESENT AND COPY THE FASTA FILE
def valid_taxid(taxid):
    """
    Check if the taxon id given by the user is in the NCBI taxonomy
    database
    """
    ncbi = NCBITaxa()
    try:
        ncbi.get_lineage(taxid)
        return taxid
    except Exception as err:
        print("Program terminated&The taxon id given ({})\
               is not in the NCBI taxonomy database !".format(taxid))
        sys.exit()


def name_organism(gb_data, accession, name):
    """
    Return the name of the organism with the GCF at the end
    """
    name = name.replace(accession, "")
    name = name.split(",")[0]
    name = name.replace("/", "_")
    name = name.replace("'", "")
    acc = re.search("^\.[0-9] ", name).group()
    name = name.replace(acc, "")
    name = name.strip()
    return name + " " + get_gcf_id(gb_data)


def get_accession_number(name):
    """
    Get the assession number of the head of fasta file
    """
    try:
        accession = re.search("([A-Z]+[\_]?[A-Za-z0-9]+)\.[0-9]*", name).group(1)

        if DEBUG: print("Accession  : " + accession)
        return accession
    except Exception as e:
        print("Program terminated&No accession number found in head of fasta file")
        sys.exit()


def get_taxon_id(gb_data):
    """
    Get taxonomy ID from a genbank data, check if this id is in the NCBI
    database and return it
    """
    try:
        taxon = gb_data.features[0].qualifiers["db_xref"]
        if DEBUG: print("Longueur taxon list  : " + str(len(taxon)))
        for tax in taxon:
            if re.search("taxon", tax):
                if DEBUG: print(tax)
                taxid = re.search("[0-9]+", tax).group()
                valid_taxid(taxid)
                if DEBUG: print("Taxon ID  :  " + taxid)
                return taxid
    except Exception as e:
        pass
    print("Program terminated&No taxonomy ID found")
    sys.exit()


def get_gcf_id(gb_data):
    """
    Get the GCF id, which is the ID for the assembly with annotations from a
    genbank data
    """
    try:
        gcf = re.search("GCF\_[0-9]+\.[0-9]+", gb_data.dbxrefs[0]).group()
        if DEBUG: print("GCF  :  " + gcf)
        return gcf
    except Exception:
        dspl.eprint("No GCF id found")
        return "None"


def get_asm_id(gcf):
    """
    Get the ASM id, which is the name of the Assembly from a GCF id
    """
    try:
        handle = Entrez.esearch(db="assembly", term=gcf)
        for line in handle:
            if re.search("<Id>", line):
                id_report = re.search("[0-9]+", line).group()
        handle = Entrez.esummary(db="assembly", id=id_report, report="full")
        for line in handle:
            if re.search("<AssemblyName>", line):
                asm = re.search("<AssemblyName>(.*)<\/AssemblyName>", line).group(1)
                if DEBUG: print("ASM    :  " + asm)
                return asm
    except Exception:
        pass
    dspl.eprint("No ASM id found")
    return "None"


def get_gcf_taxid(filename):
    """
    Get the accession Number, the GCF id, the ASM id and the taxonomy id from a
    fasta file
    """
    seq_record = SeqIO.parse(filename, "fasta")
    name = next(seq_record).description
    accession = get_accession_number(name)
    # Get the genbank data from the accession number
    Entrez.email = "example@gmail.com"
    try:
        Entrez.einfo()
        res = Entrez.efetch(db="nuccore", id=accession, rettype="gb", seq_start=1, seq_stop=1)
        gb_data = SeqIO.read(res, "genbank")
    except urllib.error.HTTPError:
        print("Program terminated&The accession number not present in NCBI")
        sys.exit()
    except urllib.error.URLError:
        print("Program terminated&No connection with the NCBI")
        sys.exit()
    except:
        print("Program terminated&Problem to parse the genbank file, Contact the support")
        sys.exit()

    # Get all neccessary informations
    gcf = get_gcf_id(gb_data)
    taxid = get_taxon_id(gb_data)
    name = name_organism(gb_data, accession, name)
    asm = get_asm_id(gcf) if gcf != "None" else name
    return gcf, asm, taxid, name
