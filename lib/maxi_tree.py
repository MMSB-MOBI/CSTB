"""
Class MaxiTree object. This tree contains: a tree object from the package ete3 with the name
of the organism, the current GCF and its taxonomy Id if it's present into the CRISPR database. This taxon ID is
referenced in the name separate by a ':' and in a node's feature 'taxon'.
The name of a node is composed like this:
    name (without ' and /) GCF_ID : taxon_ID

Several constructor:
    1. from the database
    2. from the genomre_ref_taxid.json file
    3. from a MaxiTree pickle file
    4. from the taxonomy database

Several methods are implemented to insert a new node, to write a json file of the tree
full (with taxonID) or not, to get the json string full or not...

The main function allow to create this Tree from the genomre_ref_taxid.json file and to create a
pickle MaxiTree file to insert into the taxon_tree_db
"""

import os
import sys
import re
import json
import pickle
import argparse
import datetime
import pickle
import requests
import pycouch.wrapper as couchdb
import display_result as dspl
from ete3 import Tree,NCBITaxa

SESSION = requests.session()
SESSION.trust_env = False

class MaxiTree(object):
    """docstring for MaxiTree."""
    def __init__(self, dt_tree, all_spc):
        self.tree = dt_tree
        self.all_spc = all_spc

    @classmethod
    # NO VERIFICATION : PROBABLY OBSOLETE
    def from_database(cls, end_point, db_name):
        req_func = check_connexion(end_point)
        try:
            # Try to get the MaxiTree object by the key maxi_tree
            maxi_tree = req_func.post(end_point + db_name, json={"keys": ["maxi_tree"]}).json()["request"]["maxi_tree"]
            maxi_tree = json.loads(maxi_tree["tree"].replace("'", '"'))
        except Exception as e:
            print("Program terminated&Problem with the maxi_tree in the database {}".format (e))
            sys.exit()
        # Convert the json string into a Tree object
        tree_topo = Tree()
        cls.convert_json_to_tree(cls, tree_topo, maxi_tree)
        # List all taxon ID from leaves
        list_taxid = [int(node.taxon) for node in tree_topo.iter_descendants() if hasattr(node, "taxon")]
        return cls(tree_topo, list_taxid)

    @staticmethod
    # NO VERIFICATION : PROBABLY OBSOLETE
    def convert_json_to_tree(cls, node, dic):
        # Give the name to the root
        if node.is_root() and node.name == '':
            node.name = dic["text"]
            new_node = node
        else:
        # Give the name to the node
            new_node = node.add_child(name=dic["text"])
        try:
        # Try to retrieve the taxon ID if it exists
            taxid = re.search(" : ([0-9]+)", dic["text"]).group(1)
            new_node.add_feature("taxon", taxid)
        except:
            pass
        # The current node is a leaf, so it is the end of the process for this subtree
        if not "children" in dic.keys(): return None
        # For each child, traverse it
        for i in dic["children"]:
            cls.convert_json_to_tree(cls, new_node, i)

    @classmethod
    # NO VERIFICATION : PROBABLY OBSOLETE
    def from_maxitree(cls, maxitree_file):
        # Check if the pickle file exists
        if not os.path.isfile(maxitree_file): sys.exit("*** File does not exist ***")
        try:
            tree_topo = pickle.load(open(maxitree_file, "rb"))
        except:
            sys.exit("Problem to open the pickle file {}".format(maxitree_file))
        # Generate the list of taxonID from leaves
        list_taxid = [int(node.taxon) for node in tree_topo.iter_descendants() if hasattr(node, "taxon")]
        return cls(tree_topo, list_taxid)

    @classmethod
    # NO VERIFICATION : PROBABLY OBSOLETE
    def from_gen_ref(cls, gen_ref_file, end_point, db_name):
        # Check if the genome_ref_taxid.json file exists
        if not os.path.isfile(gen_ref_file): sys.exit("*** File does not exist ***")
        try:
            dic_ref = json.load(open(gen_ref_file))
            all_taxid = [int(dic_ref[org][1]) for org in dic_ref]
        except:
            sys.exit("Problem with the format of the file {}".format(gen_ref_file))
        # Construct the Tree object
        # Add the taxonID of plasmid
        tree_topo = cls.construct_tree(cls, all_taxid + [36549], end_point, db_name)
        # Generate the list of TaxonID from leaves
        list_taxid = [int(node.taxon) for node in tree_topo.iter_descendants() if hasattr(node, "taxon")]
        return cls(tree_topo, list_taxid)

    @staticmethod
    # VERIFICATED : OK
    def construct_tree(cls, list_taxid, end_point, db_name):
        # Check if it can connect to the database
        ncbi = NCBITaxa()
        # Create the topology Tree with taxon ID
        tree_topo = ncbi.get_topology(list_taxid)
        # Add the feature taxon for leaves
        for node in tree_topo.iter_descendants():
            if int(node.name) in list_taxid:
                if node.children:
                    node.add_child(name=node.name)
                else:
                    node.add_feature("taxon", node.name)

            node.name = rename_node(node.name, ncbi)
            if hasattr(node, "taxon"):
                # Add the GCF and taxon ID to the name
                taxon_doc = json.loads(SESSION.get(end_point + db_name + "/" + node.taxon).text)
                gcf = taxon_doc["current"]
                taxondb_name = taxon_doc["current_name"]
                node.name = "{} {} : {}".format(taxondb_name, gcf, node.taxon)
        # Name for the root
        tree_topo.name = ncbi.get_taxid_translator([int(tree_topo.name)])[int(tree_topo.name)]

        return tree_topo

    @classmethod
    # VERIFICATED : OK
    def from_taxon_database(cls, end_point, db_name):
        couchdb.setServerUrl(end_point + db_name)
        if not couchdb.couchPing():
            print("Program terminated&Can't connect to the Taxon database")
            sys.exit()
        # Retrieve all taxon ID and plasmid name
        all_entries = couchdb.couchGetRequest("_all_docs")
        list_taxon = [row['key'] for row in all_entries['rows']]
        # Only keep taxon id and cast into int
        list_taxond_id = [int(i)  for i in list_taxon if re.match("^[0-9]*$", i)]
        
        # Only keep plasmid name
        list_plasmids = [i  for i in list_taxon if not re.match("^[0-9]*$", i)]
        # Construct the Tree object and Add the taxonID of plasmid
        tree_topo = cls.construct_tree(cls, list_taxond_id + ["36549"], end_point, db_name)
        # Insert plasmid under node plasmid
        for plasmid in list_plasmids:
            tree_topo = cls.insert_plasmid(cls, plasmid, end_point, db_name, tree_topo)
        # Generate the list of TaxonID from leaves
        list_taxid = [int(node.taxon) for node in tree_topo.iter_descendants() if hasattr(node, "taxon")]
        return cls(tree_topo, list_taxid)

    ### Write AND Representation of tree ###
    def __repr__(self):
        return "Number of leaves : {} from the root \"{}\"".format(len(self.all_spc), self.tree.name)

    def write_json_full(self):
        json_tree = self.get_json(True)
        with open("jsontree_full.json", "w") as filout:
            filout.write(json_tree)

    def write_json(self):
        json_tree = self.get_json(False)
        with open("jsontree.json", "w") as filout:
            filout.write(json_tree)

    def get_json(self, full):
        return str(self.__get_json__(self.tree, full))

    def __get_json__(self, node, full):
        if full:
            json = {"text": node.name}
        else:
            # Remove the taxon ID if it is present
            json = {"text": ":".join(node.name.split(":")[:-1])} if hasattr(node, "taxon") else {"text": node.name}
        # If node has children, create a list of it and traverse it, else do not create attribute
        # children and return the json string
        if node.children:
            json["children"] = []
            for ch in node.children:
                json["children"].append(self.__get_json__(ch, full))
        return json

    def dump(self, filout):
        # Get the Tree in json string with taxonID
        json_tree = self.get_json(True)
        # Create a dictionary with the MaxiTree json and the date
        dic_tree = {}
        dic_tree["maxi_tree"] = {}
        dic_tree["maxi_tree"]["tree"] = json_tree
        dic_tree["maxi_tree"]["date"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
        pickle.dump(dic_tree, open(filout, "wb"), protocol=3)

    ### Modify Tree ###
    def is_member(self, taxid):
        return taxid in self.all_spc

    def insert(self, taxid, end_point, db_name):
        # Check if can connect to the database
        req_func = check_connexion(end_point)

        # If the taxonID is already present, it is just a change of the GCF
        if self.is_member(taxid):
            # search the node in the original tree
            node = self.tree.search_nodes(taxon=str(taxid))[0]
            # Retrieve the GCF and compare it to the GCF given by the user
            gcf_old = re.search("GCF_[0-9]+\.?[0-9]*", node.name)[0]
            gcf_new = req_func.post(end_point + db_name, json={"keys": [str(taxid)]}).json()["request"][str(taxid)]["current"]
            if gcf_old == gcf_new:
                print("Program terminated&Error, the new GCF is the same than the old GCF")
                sys.exit()
            # Replace GCF
            node.name = node.name.replace(gcf_old, gcf_new)
            return True

        # Construct the topology tree with taxonID
        ncbi = NCBITaxa()
        tree_topo = ncbi.get_topology(self.all_spc + [taxid])
        # Rename node
        for node in tree_topo.iter_descendants():
        # Retrieve name from the original tree, faster than create it, just replace the GCF id
        # by the one given by the user which is already is the taxon_db
            if int(node.name) in self.all_spc:
                if node.children:
                    node.add_child(name=node.name)
                    node.name = rename_node(node.name, ncbi)
                else:
                    node.add_feature("taxon", node.name)
                    node.name = self.tree.search_nodes(taxon=node.name)[0].name
            elif node.name == str(taxid):
                if node.children:
                    node.add_child(name=node.name)
                    node.name = rename_node(node.name, ncbi)
                else:
                     try:
                         node.add_feature("taxon", node.name)
                         node.name = rename_node(node.name, ncbi)
                         gcf = req_func.post(end_point + db_name, json={"keys": [node.taxon]}).json()["request"][node.taxon]["current"]
                         node.name = "{} {} : {}".format(node.name, gcf, node.taxon)
                     except:
                         print("Program terminated&Problem with the taxon {} in the Taxon database \
    (WiP : insert a new member in the Tree)".format(taxid))
                         sys.exit()
            else:
                node.name = rename_node(node.name, ncbi)

        # Search for the node plasmid and replace it in the new tree by detaching the old plasmid
        # node without children
        plasmid_origin = self.tree.search_nodes(taxon="36549")[0]
        plasmid_new = tree_topo.search_nodes(taxon="36549")[0]
        plasmid_new.detach()
        try:
            tree_topo.add_child(plasmid_origin)
        except:
            pass

        # Name for the root
        tree_topo.name = ncbi.get_taxid_translator([int(tree_topo.name)])[int(tree_topo.name)]
        self.tree = tree_topo
        self.all_spc = [int(node.taxon) for node in tree_topo.iter_descendants() if hasattr(node, "taxon")]
        return True

    def insert_plasmid(self, name, end_point, db_name, tree):
        # Check if can connect to the database
        #req_func = check_connexion(end_point)
        #gcf = req_func.post(end_point + db_name, json={"keys": [name]}).json()["request"][name]["current"]
        taxon_doc = json.loads(SESSION.get(end_point + db_name + "/" + name).text)
        gcf = taxon_doc["current"]
        name = "{} {}".format(name, gcf)
        plasmid_node = tree.search_nodes(name="plasmids")[0]
        plasmid_node.add_child(name=name)
        new_plasmid = plasmid_node.search_nodes(name=name)[0]
        new_plasmid.add_feature("plasmid", name)
        return tree


def rename_node(taxid, ncbi):
    """
    Rename node by searching the organism name from its taxonomy ID and
    replacing single quote by space and slash by underscore
    """
    node_name = ncbi.get_taxid_translator([int(taxid)])[int(taxid)]
    node_name = node_name.replace("'", '')
    node_name = node_name.replace("/", '_')
    return node_name


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


def args_gestion():
    """
    Takes and treat arguments given
    """
    parser = argparse.ArgumentParser(description="Construct the MaxiTree and save it in the database")
    parser.add_argument("-file", metavar="<str>", required=False, help="The genome_ref_taxid.json file")
    parser.add_argument("-url", metavar="<str>", required=True, help="End_point for taxon_database")
    parser.add_argument("-dbName", metavar="<str>", required=True, help="Name for the database")
    return parser.parse_args()

if __name__ == '__main__':
    PARAM = args_gestion()

    TREE = MaxiTree.from_gen_ref(PARAM.file, PARAM.url, PARAM.dbName) if PARAM.file else MaxiTree.from_taxon_database(PARAM.url, PARAM.dbName)
    JSON_TREE = TREE.get_json(True)
    DIC_TREE = {}
    DIC_TREE["maxi_tree"] = {}
    DIC_TREE["maxi_tree"]["tree"] = JSON_TREE
    DIC_TREE["maxi_tree"]["date"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    pickle.dump(DIC_TREE, open("jsontree_full.p", "wb"), protocol=3)
