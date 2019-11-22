#!/usr/bin/env python3
"""
Post-processing after the SetCompare script. Retrieve coordinates of sgRNA with requests into
database and display it for the webservice.
"""

import sys
import os
import time
import argparse
import re
import maxi_tree as MT
import pycouch.wrapper as couchdb
import json
import operator
from collections import OrderedDict
import requests
import wordIntegerIndexing as decoding
import display_result as dspl


def valid_file(parser, filename):
    """
    Check if the file exists
    """
    # Check if the file exists
    if not os.path.isfile(filename):
        parser.error("The file {} does not exist !".format(filename))
    # If all is good, return the filename
    return filename


def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    # Argparsing
    parser = argparse.ArgumentParser(description="Post-processing results")
    # parser.add_argument("-file", metavar="FILE", type=lambda x: valid_file(parser, x),
    #                     help="The path to the fasta file")
    parser.add_argument("-sl", metavar="<int>",
                        help="The length of the sgrna, excluding pam", default=20)
    parser.add_argument("-pam", metavar="<str>",
                        help="The pam motif",
                        required=True, default="NGG")
    parser.add_argument("-gi", metavar="<str>",
                        help="The organisms to search inclusion in.",
                        required=True)
    parser.add_argument("-gni", metavar="<str>",
                        help="The organisms to search exclusion from",
                        required=True)
    parser.add_argument("-r", metavar="<str>",
                        help="The end point",
                        required=True)
    parser.add_argument("-taxon_db", metavar="<str>",
                        help="The name of the taxon database",
                        required=True)
    parser.add_argument("-tree_db", metavar="<str>",
                        help="The name of the taxon database",
                        required=True)
    parser.add_argument("-end_point", metavar="<str>",
                        help="The end point of the taxon and tree database",
                        required=True)
    parser.add_argument("-c", metavar="<int>",
                        help="The length of the slice for the request",
                        required=True)
    parser.add_argument("-f", metavar="<str>",
                        help="The file index",
                        required=True)
    parser.add_argument('--no-proxy', action='store_true')
    parser.add_argument("-nb_top", metavar="<int>",
                        help="The top hits to download",
                        default=1000)
    return parser.parse_args()


def couchdb_search(sgrna_list, end_point, len_slice, no_poxy_bool):
    """
    Check if it can connect to the dabase
    Requests a list of sgrna and return the result
    """
    req_func = requests
    if no_poxy_bool:
        req_func = requests.Session()
        req_func.trust_env = False

    joker = 0
    results = {"request": {}}
    try:
        res = req_func.get(end_point + "/handshake")
        data = res.json()
        dspl.eprint("HANDSHAKE PACKET:\n" + str(data) + "\n")
    except Exception as e:
        dspl.eprint("Could not perform handshake, exiting")
        print("Program terminated&No handshake")
        sys.exit(1)

    while True:
        try:
            for i in range(0, len(sgrna_list), len_slice):
                dspl.eprint(i)
                request_sliced = {"keys" :sgrna_list[i : i + len_slice]}
                results["request"].update(req_func.post(end_point + "/bulk_request",
                                                        json=request_sliced).json()["request"])
        except Exception as e:
            dspl.eprint("Something wrong append '{}', retrying time {}".format(end_point, joker))
            dspl.eprint("Error LOG is ", str(e))
            joker += 1
            if joker > 3:
                print("Progam terminated&after 50 tries to access to the database")
                sys.exit(1)
            time.sleep(5)
            continue
        break

    # dspl.eprint("motif-broker ans:\n", str(results))
    return results


def treat_db_search_20(dic_hits, genomes_in, end_point, len_slice, no_poxy_bool):
    """
    Treat directly sequences with 20 length and return an update of a dictionary containing
    Hit objects
    """
    results = couchdb_search(list(dic_hits.keys()), end_point, len_slice, no_poxy_bool)
    for sgrna in results["request"]:
        dic_seq = {}
        for org_name in results["request"][sgrna]:
            nobackslash_org_name = org_name.replace('/', '_')
            if nobackslash_org_name in genomes_in:
                dic_seq[nobackslash_org_name] = results["request"][sgrna][org_name]
        dic_hits[sgrna].set_genomes_dict(dic_seq)

    return dic_hits


def treat_db_search_other(dic_hits, dic_index, genomes_in, end_point, len_slice, no_poxy_bool, len_seq):
    """
    Find the result for word_20 associated to a word_15. Then merge results and update the
    Hit object
    """
    sgrna_list = [word for w15 in dic_index for word in dic_index[w15]]
    results = couchdb_search(sgrna_list, end_point, len_slice, no_poxy_bool)
    # Loop on word_15 to find their word_20 associated
    for sgrna in dic_hits:
        dic_seq = {}
        for word_20 in dic_index[sgrna]:
            for org_name in results["request"][word_20]:
                # To remove when the database is clean
                nobackslash_org_name = org_name.replace('/', '_')
                if nobackslash_org_name in genomes_in:
                    results["request"][word_20][org_name] = update_coord(results["request"][word_20][org_name], len_seq)
                    if nobackslash_org_name in dic_seq:
                        dic_seq[nobackslash_org_name] = merge_dic(dic_seq[nobackslash_org_name], results["request"][word_20][org_name])
                    else:
                        dic_seq[nobackslash_org_name] = results["request"][word_20][org_name]
        dic_hits[sgrna].set_genomes_dict(dic_seq)
    return dic_hits


def update_coord(dic_results, len_seq):
    """
    Update coordinates with the length of sgRNA
    """
    offset = 23 - len_seq
    return {ref: [replace_coord("[+-]\(([0-9]*),", operator.add, coord, offset) if coord[0] == "+" else replace_coord(",([0-9]*)", operator.sub, coord, offset) for coord in dic_results[ref]] for ref in dic_results}


def replace_coord(regex, op_func, coord, offset):
    """
    Change the coordinate string
    """
    sgrna_start = int(re.search(regex, coord).group(1))
    return coord.replace(str(sgrna_start), str(op_func(sgrna_start, offset)))


def merge_dic(dic_seq_ref, dic_results):
    """
    Merge dictionaries. Because several word_20 are associated to word_15
    so the same reference can be write
    """
    for ref in dic_results:
        if ref in dic_seq_ref.keys():
            dic_seq_ref[ref] += dic_results[ref]
        else:
            dic_seq_ref[ref] = dic_results[ref]
    return dic_seq_ref


def check_find_database(dic_hits):
    """
    Check if all sgrna were found in database
    """
    for sgrna in dic_hits:
        if not dic_hits[sgrna].genomes_dict:
            dspl.eprint("Not find in database : {}".format(sgrna))


def parse_setcompare_other(output_c, nb_top):
    """
    Parse the output of setCompare C script and return a dictionary with the coding word
    and its occurence if the length is 20 or a tuple containing its occurence and
    a list of word_20 associated
    """
    with open(output_c, "r") as filin:
        for line in filin:
            regex_nb_hits = re.search("^# ([0-9]+)", line)
            if regex_nb_hits:
                nb_hits = int(regex_nb_hits.group(1))
                if nb_hits == 0:
                    print("Program terminated&No hits remain")
                    sys.exit()
                break

        index_dic = OrderedDict()
        i = 0
        for rank_occ in filin:
            if i == nb_top or rank_occ == "\n": break
            rank_splitted = rank_occ.split(":")
            rankw20_occ = rank_splitted[1].split("[")
            index_dic[int(rank_splitted[0])] = [rankw20_occ[0], rankw20_occ[1][:-2].split(",")]
            i += 1
    return index_dic, nb_hits


def parse_setcompare_out(output_c, nb_top):
    """
    Parse the output of setCompare C script and return a dictionary with the coding word
    and its occurence if the length is 20 or a tuple containing its occurence and
    a list of word_20 associated
    """
    with open(output_c, "r") as filin:
        text = filin.readlines()

    nb_hits = re.search("[0-9]+", text[-2]).group()
    if int(nb_hits) == 0:
        print("Program terminated&No hits remain")
        sys.exit()

    index_dic = OrderedDict()
    i = 0
    for rank_occ in text[-1].strip().split(","):
        if i == nb_top: break
        index_dic[int(rank_occ.split(":")[0])] = rank_occ.split(":")[1]
        i += 1
    return index_dic, nb_hits


def create_dic_hits(param, genomes_in):
    """
    Treat the search into database and return a dictionary of hits
    """
    dic_index, nb_hits = parse_setcompare_out(param.f, int(param.nb_top)) if int(param.sl) == 20 else parse_setcompare_other(param.f, int(param.nb_top))
    dspl.eprint("NB_HITS  ==>  {}     Length DIC_INDEX  ==>  {}".format(nb_hits, len(dic_index)))

    len_seq = int(param.sl) + int(len(param.pam))
    dic_hits = OrderedDict()
    dic_new = {}
    # Decoding of each index into sequence
    for rank in dic_index:
        sequence = decoding.decode(rank, ["A", "T", "C", "G"], len_seq)
        occ = dic_index[rank] if int(param.sl) == 20 else dic_index[rank][0]
        if int(param.sl) != 20:
            dic_new[sequence] = [decoding.decode(int(w20), ["A", "T", "C", "G"], 23) for w20 in dic_index[rank][1]]
        dic_hits[sequence] = dspl.Hit(sequence, occ)

    # Search coordinates for each sgrna in each organism
    if int(param.sl) == 20:
        dic_hits = treat_db_search_20(dic_hits, genomes_in, param.r, int(param.c),
                                      param.no_proxy)
    else:
        dic_hits = treat_db_search_other(dic_hits, dic_new, genomes_in, param.r,
                                         int(param.c), param.no_proxy, len_seq)
    return dic_hits


def get_size(end_point, tree_db_name, taxon_db, genomes_in):
    tree = MT.MaxiTree.from_database(end_point, tree_db_name)
    json_tree = tree.get_json(True)
    taxon_orgname = {}

    req_func = requests.Session()
    req_func.trust_env = False
    try:
        res = req_func.get(end_point + "handshake").json()
        dspl.eprint("HANDSHAKE PACKET (taxon_db) : {}".format(res))
    except Exception as e:
        dspl.eprint("Could not perform handshake, exiting")
        print("Program terminated&No handshake with taxon database")
        sys.exit(1)

    try:
        for org in genomes_in:
            taxon_orgname[re.search(org.replace("(", "\(").replace(")", "\)") + " : ([^']*)", json_tree).group(1)] = org

        size_dic = {}
        joker = 0
        while joker < 3 :
            try:
                res = req_func.post(end_point + taxon_db ,json={"keys" : list(taxon_orgname.keys())}).json()["request"]
                for taxon in res :
                    size_dic[taxon_orgname[taxon]] = res[taxon]["size"]
                break
            except:
                dspl.eprint("something wrong append '{}', {} times".format(list(taxon_orgname.keys()), joker))
                joker += 1

    except:
        # print("Program terminated&Error with genomes size")
        json.dump("", open("size_org.json", "w"))
        # sys.exit()

    json.dump(size_dic, open("size_org.json", "w"))


if __name__ == '__main__':
    PARAM = args_gestion()
    GENOMES_IN = PARAM.gi.split("&")
    GENOMES_NOTIN = PARAM.gni.split("&")
    DIC_HITS = create_dic_hits(PARAM, GENOMES_IN)

    # Display the result for the navigator
    dspl.display_hits(DIC_HITS, GENOMES_IN, GENOMES_NOTIN,
                      PARAM.pam, int(PARAM.sl), ".", int(PARAM.nb_top),
                      True, list(DIC_HITS.keys()))

    get_size(PARAM.end_point, PARAM.tree_db, PARAM.taxon_db, GENOMES_IN)

    print(','.join(GENOMES_NOTIN))
    print("TASK_KEY")
    print(len(DIC_HITS))
