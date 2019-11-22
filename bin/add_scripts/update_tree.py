#!/usr/bin/env python3
"""
Load the MaxiTree object from the taxon_tree_db and insert a new member.
Then, create a pickle file to insert it into the database
"""

import os
import argparse
import requests
import maxi_tree as mt


def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    parser = argparse.ArgumentParser(description="Update the MaxiTree object from the database")
    parser.add_argument("-url", metavar="<str>", required=True,
                        help="The endpoint of the tree database")
    parser.add_argument("-treeName", metavar="<str>", required=True,
                        help="Name of the tree database where the MxiTree object is saved")
    parser.add_argument("-taxonDB", metavar="<str>", required=True,
                        help="The endpoint to the taxon database")
    parser.add_argument("-taxonName", metavar="<str>", required=True,
                        help="Name of the taxon database")
    parser.add_argument("-taxid", metavar="<str>",
                        help="Taxonomy ID to update or to add")
    parser.add_argument("-name", metavar="<str>",
                        help="Name of plasmid to add")
    parser.add_argument('--no-proxy', action='store_true')
    return parser.parse_args()


if __name__ == '__main__':
    PARAM = args_gestion()
    print(PARAM)
    if PARAM.no_proxy:
        req_func = requests.Session()
        req_func.trust_env = False

    MAXITREE = mt.MaxiTree.from_database(PARAM.url, PARAM.treeName)
    MAXITREE.insert(PARAM.taxid, PARAM.taxonDB, PARAM.taxonName) if PARAM.taxid else MAXITREE.insert_plasmid(PARAM.name, PARAM.taxonDB, PARAM.taxonName)

    try:
        os.mkdir("./treeDB_data/")
    except OSError:
        print("Be careful : The directory treeDB_data exists")

    MAXITREE.dump("./treeDB_data/maxi_tree.p")
