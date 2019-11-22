import sys
from ete3 import Tree,NCBITaxa
import json

"""
Last change : 5pm 11 dec. 2016

Run thois script with the commande below:

	python3 tax2json.py <tax_file>

The <tax_file> must be a text file witch has content like this:

	Buchnera aphidicola str. APS (Acyrthosiphon pisum)
	Borrelia burgdorferi B31
	Treponema denticola ATCC 35405
	<NCBI taxonomy name>
	...

The output file will be named 'json_tree.json' as json format or
flat file.

* Part of script adapted from https://gist.github.com/jhcepas/9205262
"""

def get_json(node):
	node.name = node.name.replace("'", '')
	json = {"text": node.name}
	if node.children:
		json["children"] = []
		for ch in node.children:
			json["children"].append(get_json(ch))
	return json

if __name__ == '__main__':

	Dict=json.load(open(sys.argv[1] + '/genome_ref_taxid.json','r'))	##This Dictionary contains as keys the names of all the bacterial genomes in the database.

	ncbi = NCBITaxa()
	allspe = []

	list_ID=[]

	invert_dic={}

	for sp in Dict:
		ID = int(Dict[sp][1])
		allspe.append(ID)
		invert_dic[ID]=sp

	# Build topologic tree
	treeTopo = ncbi.get_topology(allspe)


	# Convert to Newick tree as string format
	treeNwk = treeTopo.write(format_root_node=True, format=8)
	#print(treeNwk)

	# build tree object from Newick format string
	newTree = Tree(treeNwk, format=1)	# in = str
	# get nodes and leaves name
	for i in newTree.iter_descendants():
		if int(i.name) in invert_dic:
			name=invert_dic[int(i.name)]
			i.name=name
			#i.name=i.name.split(' ')[-2:]
		else:
			i.name = ncbi.get_taxid_translator([int(i.name)])[int(i.name)]
	# get root's name
	newTree.name = ncbi.get_taxid_translator([int(newTree.name)])[int(newTree.name)]


	# Convert to json format
	json_tree = str(get_json(newTree)).replace("'", '"')

	# output to json file (or .txt format)
	output = open(sys.argv[2], 'w')
	output.write(json_tree)
	output.close()
