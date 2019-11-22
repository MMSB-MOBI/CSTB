"""
Create pickle file to inset into Taxonomy database
It can create from a GCF and Taxonomy ID given, from the json file genome_ref_taxid.json file or
from a list file in csv format

OUTUPUT:
taxon_dt={"1234": {"GCF": ["GCF_111", "GCF_112", "GCF_113"], "date": "19-04-2019", "User": "Queen"},
          "2345": {"GCF": ["GCF_211", "GCF_212", "GCF_213"], "date": "19-04-2019", "User": "Queen"}}
"""

import os
import sys
import datetime
import json
import pickle
import argparse
import csv
import requests
import display_result as dspl
from Bio import SeqIO


def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    command = sys.argv[1] if len(sys.argv) > 1 else None
    parser = argparse.ArgumentParser(description="Convert the genomre_ref_taxid.json to a pickle\
                                                  file ready to be insert into database")
    subparsers = parser.add_subparsers(help="commands")

    # Parser from SCRATCH
    scratch_parser = subparsers.add_parser("scratch", help="Create files to insert from\
                                                            the genomre_ref_taxid.json file")
    scratch_parser.add_argument("-file", metavar="<str>", required=True,
                                help="The json file to convert or the list file in csv")
    scratch_parser.add_argument("-user", metavar="<str>", nargs="?", help="The user",
                                default="MMSB")
    scratch_parser.add_argument("-out", metavar="<str>", nargs="?",
                                help="The name of the outputfile", default="taxon_dt.p")
    scratch_parser.add_argument("-rfg", metavar="<str>", required=True,
                                help="Path of the fasta file")
    scratch_parser.add_argument('--no-proxy', action='store_true')
    scratch_parser.add_argument('--json', action='store_true')

    # Parser for a SINGLE genome
    single_parser = subparsers.add_parser("single", help="Create file to insert")
    single_parser.add_argument("-user", metavar="<str>", help="The user", nargs="?", default="MMSB")
    single_parser.add_argument("-gcf", metavar="<str>", required=True, help="GCF id")
    single_parser.add_argument("-taxid", metavar="<str>", required=True, help="Taxonomy id")
    single_parser.add_argument("-r", metavar="<str>", required=True,
                               help="End point for taxon Database")
    single_parser.add_argument("-dbName", metavar="<str>", required=True,
                               help="Name of the taxon Database")
    single_parser.add_argument("-out", metavar="<str>", nargs="?",
                               help="The name of the outputfile", default="taxon_dt.p")
    single_parser.add_argument("-rfg", metavar="<str>", required=True,
                                help="Path of the fasta file")

    return parser.parse_args(), command


def set_env(json_file=None):
    """
    Create a folder taxonDB_data and check if the json_file given exists
    """
    try:
        os.mkdir("./taxonDB_data/")
    except OSError:
        print("Be careful : The directory taxonDB_data exists")

    if json_file and not os.path.isfile(json_file):
        sys.exit("The file does not exist, check your path")


def size_fasta(fasta_file):
    tmp = {}
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            tmp[record.id] = len(record.seq)
    except:
        with open(".taxon.log", "a") as filout:
            filout.write(fasta_file + "\n")
        return 1
    return tmp


def init_taxondt(gcfs, user, taxon_id, fasta_path, gcf_given):
    """
    Initialize the dictionary of taxon and return it
    """
    tmp_dic = {}
    tmp_dic["size"] = size_fasta(fasta_path + taxon_id + "_" + gcf_given + ".fna")
    if(tmp_dic["size"] == 1): return 1
    tmp_dic["GCF"] = gcfs
    tmp_dic["date"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    tmp_dic["user"] = user
    tmp_dic["current"] = gcfs[0]
    return tmp_dic


def load_json_file(file_name, user, fasta_path) :
    taxon_dt = {}
    duplicate = []
    dic_ref = json.load(open(file_name, "r"))
    for org in dic_ref:
        taxid = dic_ref[org][1]
        gcf = "_".join(dic_ref[org][0].split("_")[:2])
        list_gcf = [gcf] if taxid not in taxon_dt.keys() else taxon_dt[taxid]["GCF"] + [gcf]
        if len(list_gcf) > 1: duplicate.append(taxid)
        tmp_taxon_dt = init_taxondt(list_gcf, user, taxid, fasta_path, gcf)
        if (tmp_taxon_dt != 1):
            taxon_dt[taxid] = tmp_taxon_dt

    if duplicate:
        with open("duplicate_taxon.log", "w") as filout:
            [filout.write("{}\n".format(tax)) for tax in duplicate]

    return taxon_dt


def load_list_file(file_name, user, fasta_path):
    taxon_dt = {}
    duplicate_taxon = []
    try:
        filin = open(file_name, 'r')
        content = csv.reader(filin, delimiter="\t")
        for line in content:
            taxid = line[0]
            gcf = line[1]
            list_gcf = [gcf] if taxid not in taxon_dt.keys() else taxon_dt[taxid]["GCF"] + [gcf]
            if len(list_gcf) > 1: duplicate_taxon.append(taxid)
            tmp_taxon_dt = init_taxondt(list_gcf, user, taxid, fasta_path, gcf)
            if (tmp_taxon_dt != 1): taxon_dt[taxid] = tmp_taxon_dt

        if duplicate_taxon:
            with open("duplicate_taxon.log", "w") as filout:
                [filout.write("{}\n".format(tax)) for tax in duplicate_taxon]
    except:
        sys.exit("Problem with file.\n{}\n* FORMAT *\n{}\ntaxId{}GCF_id".format("*"*10, "*"*10, r"\t"))
    filin.close()
    return taxon_dt


def check_connexion(end_point):
    req_func = requests.Session()
    req_func.trust_env = False
    try:
        res = req_func.get(end_point + "handshake").json()
        dspl.eprint("HANDSHAKE PACKET (tree_db) : {}".format(res))
        return req_func
    except Exception as e:
        dspl.eprint("Could not perform handshake, exiting")
        print("Program terminated&No handshake with taxon database")
        sys.exit(1)


if __name__ == '__main__':
    PARAM, COMMAND = args_gestion()

    # Insert from SCRATCH
    if COMMAND == "scratch":
        set_env(PARAM.file)
        TAXON_DT = load_json_file(PARAM.file, PARAM.user, PARAM.rfg) if PARAM.json else load_list_file(PARAM.file, PARAM.user, PARAM.rfg)

    # Insert for a SINGLE genome
    else:
        set_env()
        TAXON_DT = {}
        req_func = check_connexion(PARAM.r)
        DOC = req_func.post(PARAM.r + PARAM.dbName, json={"keys": [PARAM.taxid]}).json()["request"]
        LIST_GCF = [PARAM.gcf] + DOC[PARAM.taxid]["GCF"] if DOC else [PARAM.gcf]
        tmp_taxon_dt = init_taxondt(LIST_GCF, PARAM.user, PARAM.taxid, PARAM.rfg, PARAM.gcf)
        if(tmp_taxon_dt != 1): TAXON_DT[PARAM.taxid] = tmp_taxon_dt

    pickle.dump(TAXON_DT, open("./taxonDB_data/" + PARAM.out, "wb"), protocol=3)
