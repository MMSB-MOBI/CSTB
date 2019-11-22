#!/usr/bin/env python3
"""
Create a folder with fasta files and necessary informations to add them into the database
It can do it from scratch with the json tree file and the name of the node or with the
MaxiTree from the database and the name of the node
"""

import os
import sys
import re
import json
import argparse
import shutil
import requests
import maxi_tree as mt


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
    command = sys.argv[1] if len(sys.argv) > 1 else None
    parser = argparse.ArgumentParser(description="Retrieve all children from a given node and\
                                                 create configure file to add them to the database")
    subparsers = parser.add_subparsers(help="commands")

    # Parser from json tree file and genome_ref_taxid file
    parser_scratch = subparsers.add_parser("scratch", help="From the genome_ref file and the\
                                                           json tree file")
    parser_scratch.add_argument("-rfg", metavar="<str>",
                                help="The path to the reference genome folder",
                                required=True)
    parser_scratch.add_argument("-tree", metavar="<str>", type=lambda x: valid_file(parser_scratch, x),
                                help="The path to the json tree",
                                required=True)
    parser_scratch.add_argument("-node", metavar="<str>",
                                help="The node to search after",
                                required=True)
    parser_scratch.add_argument("-dir", metavar="<str>",
                                help="The workdir where configure files will be created",
                                nargs="?", default="./")

    # Parser from databases
    parser_db = subparsers.add_parser("db", help="From the genome_ref file and the\
                                                           json tree file")
    parser_db.add_argument("-rfg", metavar="<str>",
                           help="The path to the reference genome folder",
                           required=True)
    parser_db.add_argument("-tree", metavar="<str>",
                           help="End point to the Tree Database",
                           required=True)
    parser_db.add_argument("-treeName", metavar="<str>",
                           help="Name of the Tree Database",
                           required=True)
    parser_db.add_argument("-taxon", metavar="<str>",
                           help="End point to the Taxon Database",
                           required=True)
    parser_db.add_argument("-taxonName", metavar="<str>",
                           help="Name of the Taxon Database",
                           required=True)
    parser_db.add_argument("-node", metavar="<str>",
                           help="The node to search after",
                           required=True)
    parser_db.add_argument("-dir", metavar="<str>",
                           help="The workdir where configure files will be created",
                           nargs="?", default="./")

    return parser.parse_args(), command


def create_folder(path_folder):
    """
    Create folder to create configure files
    """
    if os.path.exists(path_folder):
        shutil.rmtree(path_folder)
    os.makedirs(path_folder)


def search_subtree(tree, value):
    """
    Return a subtree from a given node
    """
    # Node name is the searched node
    if tree["text"] == value: return tree
    # No children, so return None
    if "children" not in list(tree.keys()): return None
    # Traverse the tree via children
    for subref in tree["children"]:
        res = search_subtree(subref, value)
        if res != None:
            return res


def search_leaves(tree, list_child):
    """
    Return a list of all children from a given tree
    """
    # The node has not children, so it is a leave else traverse children
    if "children" not in list(tree.keys()):
        list_child.append(tree["text"])
        return list_child
    # Traverse the tree via children
    for subref in tree["children"]:
        list_child = search_leaves(subref, list_child)
    return list_child


def configfile(list_leaves, rfg, path_folder):
    """
    Copy fasta file and create config file for each genome to add
    """
    path_fasta = rfg + "/genome_fasta/"

    try:
        with open(rfg + "/genome_ref_taxid.json", "r") as json_data:
            dic_ref = json.load(json_data)
    except:
        sys.exit("No genome_ref_taxid file in the path {}, can't copy fasta file".format(rfg))


    for leaf in list_leaves:
        try:
            ref = dic_ref[leaf][0]
            gcf = re.search("GCF_[0-9]+\.[0-9]+", ref)[0]
            taxid = dic_ref[leaf][1]
            # Copy the fasta file into the tmp folder
            name_fasta = taxid + "_" + gcf + ".fna"
            shutil.copy(path_fasta + name_fasta, path_folder)
            # Write it into config_file
            with open(path_folder + taxid + "_" + gcf, "w") as config_file:
                config_file.write("GCF\tTaxon ID\n")
                config_file.write("{}\t{}".format(gcf, taxid))
        except KeyError:
            print("The organism {} is not in the genome_ref_taxid file, can't find its reference".format(leaf))
        except FileNotFoundError:
            print("No fasta file found for {} at {}".format(leaf, path_fasta + name_fasta)
        except:
            print("Problem with the organism {}".format(leaf))


def configfile_from_db(list_leaves, rfg, end_point, db_taxon, path_folder):
    """
    Write a config file with the taxonID and the GCF for each organism to add them in the database
    """
    req_func = requests.Session()
    req_func.trust_env = False
    try:
        res = req_func.get(end_point + "handshake").json()
        dspl.eprint("HANDSHAKE PACKET (tree_db) : {}".format(res))
        return True
    except Exception as e:
        dspl.eprint("Could not perform handshake, exiting")
        print("Program terminated&No handshake with taxon database")
        sys.exit(1)

    error = []
    path_fasta = rfg + "/genome_fasta/"
    for leaf in list_leaves:
        taxid = leaf.split(":")[-1]
        doc = req_func.post(end_point + db_taxon, json={"keys": [str(taxid)]}).json()["request"][str(taxid)]
        if not doc:
            error.append(taxid)
        else:
            gcf = doc["current"]
            with open(path_folder + taxid + "_" + gcf, "w") as config_file:
                config_file.write("GCF\tTaxonID\n")
                config_file.write("{}\t{}".format(gcf, taxid))
            shutil.copy(path_fasta + taxid + "_" + gcf + ".fna", path_folder)
    if error:
        with open("config_file.log", "w") as filout:
            [filout.write(err + "\n") for err in error]


def load_tree(url_tree, db_tree):
    """
    Load the tree from the database and return it under a dictionary
    """
    MAXITREE = mt.MaxiTree.from_database(url_tree, db_tree)
    # Get the tree with taxonID to find the GCF in taxon database
    json_tree = MAXITREE.get_json(True)
    return json.loads(json_tree)


if __name__ == '__main__':
    PARAM, COMMAND= args_gestion()
    PATH_FOLDER = PARAM.dir + "/genomes_add/"
    create_folder(PATH_FOLDER)
    TREE = json.load(open(PARAM.tree, "r")) if COMMAND == "scratch" else load_tree(PARAM.tree, PARAM.treeName)
    SUBTREE = search_subtree(TREE, PARAM.node)
    if not SUBTREE:
        sys.exit("No subtree found, check the name of the node")
    LIST_LEAVES = search_leaves(SUBTREE, [])
    with open("genomes_from_node.log", "w") as filout:
        for leaf in LIST_LEAVES:
            filout.write(leaf + "\n")
    configfile(LIST_LEAVES, PARAM.rfg, PATH_FOLDER) if COMMAND == "scratch" else configfile_from_db(LIST_LEAVES, PARAM.rfg, PARAM.taxon, PARAM.taxonName, PATH_FOLDER)
