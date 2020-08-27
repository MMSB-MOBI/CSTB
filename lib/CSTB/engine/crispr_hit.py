#!/usr/bin/env python3

import CSTB_core.engine.wordIntegerIndexing as decoding
import logging
import CSTB.utils.error as error
import operator
import re
import copy

class Occurence():
    def __init__(self, sgRNA, genome, fasta_header, coords):
        self.sgRNA = sgRNA
        self.genome = genome
        self.fasta_header = fasta_header
        self.coords = coords

    def keepOnGene(self, genes):
        """Keep coordinates from occurence that are on corresponding homolog gene. Return None if occurence is not on gene.
        
        :param genes: Blast hits for included genomes
        :type genes: List[BlastHit]
        :raises error.SeveralGenes: raise if we have more than one hit for one genome and one fasta_header. Probably needs to be handled.
        :return: list of new Occurence with just on gene coordinates if the current occurence is on gene
        :rtype: List[Occurence] or None
        """

        homolog_genes = [ gene for gene in genes if gene.org_uuid == self.genome ]
        if not homolog_genes:
            return None
        
        list_new_occurences = []
        for homolog_gene in homolog_genes:
            new_coords = []
            for coord in self.coords: 
                if isOnGene(coord, homolog_gene):
                    new_coords.append(coord)
            if new_coords:
                list_new_occurences.append(Occurence(self.sgRNA, self.genome, self.fasta_header, new_coords))
        return list_new_occurences

    def storeOnGene(self, genes):
        homolog_genes = [ gene for gene in genes if gene.org_uuid == self.genome ]
        for coord_obj in self.coords:
            nb_homolog = 0
            for homolog_gene in homolog_genes:
                nb_homolog += 1
                gene_id = f"homolog_gene_{nb_homolog}"
                if isOnGene(coord_obj["coord"], homolog_gene):
                    coord_obj["is_on_gene"].append(gene_id)

    def updateCoords(self, bp_to_del:int):
        def replace_coord(regex, op_func, coord, offset):
            sgrna_start = int(re.search(regex, coord).group(1))
            new_coord = coord.replace(str(sgrna_start), str(op_func(sgrna_start, offset)))
            return new_coord

        coords_save = copy.deepcopy(self.coords)
        self.coords = []
        for coord_obj in coords_save: 
            coord = coord_obj["coord"]
            new_coord = replace_coord("[+-]\(([0-9]*),", operator.add, coord, bp_to_del) if coord[0] == "+" else replace_coord(",([0-9]*)", operator.sub, coord, bp_to_del)
            self.coords.append({"coord":new_coord, "is_on_gene":[]})


class Hit():
    """Hit object
    
    :ivar index: setCompare index
    :vartype index: int
    :ivar weight: setCompare weigth, correspond to number of sgRNA occurences in all genomes
    :vartype weight: int
    :ivar sequence: sgRNA sequence in nucleotides
    :vartype sequence: str
    :ivar occurences: Occurences of sgRNA in genomes
    :vartype occurences: Dict { genome : { fasta_header : List[str] } }
    :ivar index_longer: setCompare index for 20-length words
    :vartype index_longer: List[int]
    :ivar codec: codec used to encode indexes
    :vartype codec: twobits|pow2
   
    """
    
    def __init__(self, index, weight, len_sgrna, longer_index = [], codec = "twobits"):
        """ Initialize an Hit object
        
        :param index: setCompare index
        :type index: int
        :param weight: setCompare weight, correspond to number of sgRNA occurences in all genomes
        :type weight: int
        """
        self.index = int(index)
        self.weight = weight
        self.sequence = self.decode(len_sgrna, codec)
        self.list_occurences = []
        #self.on_gene_occurences = []
        self.len_sgrna = len_sgrna
        self.longer_index = longer_index
        self.to_request_sequences = self.decode_longer(codec) if self.longer_index else [self.sequence]
        self.on_gene = []
        self.not_on_gene = []

    @property
    def occurences(self): #This dictionnary will be {"organism": {"subsequence" : coords[]}}
        """A way to keep old comportment for the creation of output dictionnaries.
        """
        self._format_occurences()
        return self.formated_occurences

    def _format_occurences(self):
        occurences = {}
        for occ in self.list_occurences:
            if occ.genome not in occurences: 
                occurences[occ.genome] = {}
            
            dic = {'coords' : occ.coords}
            occurences[occ.genome][occ.fasta_header] = dic
            
        self.formated_occurences = occurences

    def _list_ref(self, org_name):
        """Format occurences for an organism 
        """
        list_ref = []
        for ref in self.occurences[org_name]:
            dic = {"ref" : ref, "coords" : self.occurences[org_name][ref]["coords"]}
            list_ref.append(dic)
        return list_ref

    def list_occ(self, taxon_name):
        """Format occurences for an organism
        """
        list_occurences = [{'org': taxon_name[genome], 'all_ref': self._list_ref(genome)} for genome in self.occurences]
        return list_occurences

    @property
    def number_occurences(self):
        """Give the total number of hit occurences
        
        :return: Total number of hit occurences
        :rtype: int
        """
        nb_occ = 0
        for genome in self.occurences:
            for subseq in self.occurences[genome]:
                nb_occ += len(self.occurences[genome][subseq]["coords"])
        return nb_occ

    def __str__(self):
        return f"index:{self.index}\nweight:{self.weight}\nsequence:{self.sequence}\noccurences:{self.occurences}\nindex_longer:{self.longer_index}\nto_request_sequences:{self.to_request_sequences}"

    def store_occurences(self, list_couch_doc, genomes_in):
        """Set occurences attributes from couch document
        
        :param list_couch_doc: couch responses
        :type couch_doc: List[Dict]
        :param genomes_in: list of genomes uuid
        :type genomes_in: List[str]
        :raises error.ConsistencyError: Raise if at least 1 given genome is not in couch document. 
        """
        logging.debug(f"Store couch doc\n {list_couch_doc}")

        # If we have more than 1 couch document, merge into one
        if len(list_couch_doc) > 1:
            couch_doc = self.merge_occurences(list_couch_doc)
        else:
            couch_doc = list_couch_doc[0]
        
       
        logging.debug(f"After merge :\n{couch_doc}")

        db_genomes = set(couch_doc.keys())
        if not set(genomes_in).issubset(db_genomes):
            raise error.ConsistencyError(f"Consistency error. Genomes included {genomes_in} are not in couch database for {self.sequence}")

        for genome in genomes_in:  
            for fasta_header in couch_doc[genome]:
                formated_coords = []
                for coord in couch_doc[genome][fasta_header]:
                    formated_coords.append({"coord":coord, "is_on_gene":[]})
                self.list_occurences.append(Occurence(self, genome, fasta_header, formated_coords))

        if self.longer_index:
            self._update_coords()

    def _update_coords(self):
        offset = 23 - self.len_sgrna
        for occ in self.list_occurences:
            occ.updateCoords(offset)

    def merge_occurences(self, list_couchdoc):
        """Merge a list of couch document into one
        
        :param list_couchdoc: List of couch document corresponding to sgRNA entry (with just organism key and their values)
        :type list_couchdoc: List[ {org(str) : {fasta_header(str) : coords(List[str]) } } ]
        :return: Merged couch document
        :rtype: Dict {org(str) : {fasta_header(str) : coords(List[str]) } }
        """
        logging.debug("Merge occurences")
        merged_couchdoc = {}
        for couchdoc in list_couchdoc:
            for genome in couchdoc:
                if genome not in merged_couchdoc:
                    merged_couchdoc[genome] = {}
                for fasta_header in couchdoc[genome]:
                    if fasta_header not in merged_couchdoc[genome]:
                        merged_couchdoc[genome][fasta_header] = []
                    merged_couchdoc[genome][fasta_header] += couchdoc[genome][fasta_header]
        return merged_couchdoc

    def decode(self, len_seq, codec):
        """Decode setCompare index into nucleotide sequence
        
        :return: sgRNA sequence
        :rtype: str
        """
        return decoding.decode(self.index, len_seq, codec)

    def decode_longer(self, codec, len_seq = 23):
        """Decode words from self.longer_index
        
        :param len_seq: Word length, defaults to 23
        :type len_seq: int, optional
        :return: list of sgRNA sequence
        :rtype: List[str]
        """
        if not self.longer_index:
            return []
        return [decoding.decode(index, len_seq, codec) for index in self.longer_index]


    def storeOnGeneOccurences(self, genes):
        """Fill on gene occurences. From current occurences, just keep the occurences that are on gene with just corresponding coordinates.
        
        :param genes: List of blast hits corresponding to homolog genes of included genomes.
        :type genes: List[BlastHit]
        """
        for occ in self.list_occurences:
            occ.storeOnGene(genes)

def write_to_file(genomes_in, genomes_not_in, dic_hits, pam, non_pam_motif_length, workdir, nb_top, hit_obj, list_ordered):
    """
    Write results in a file.
    The file is a tabulated file, with first column=sgrna sequence,
    then one column for each included organisms with list of positions
    of the sequence in each.

    OBSOLETE, maybe we need this if we keep the option to download results.
    """
    eprint('\n\n-- Write results to file --')
    rep_rslt_file = workdir + '/results_allgenome.txt'
    output = open(rep_rslt_file, 'w')
    gen_i = ','.join(genomes_in)
    gen_ni = ','.join(genomes_not_in)
    if gen_ni == '':
        gen_ni = 'None'
    output.write('#ALL GENOMES\n#Genomes included :' + gen_i +
                 ' ; Genomes excluded :' + gen_ni + '\n'+'#Parameters: pam:' +
                 pam + ' ; sgrna size:' + str(non_pam_motif_length) + '\n')
    output.write('sgrna sequence')
    for genome_i in genomes_in:
        output.write('\t' + genome_i)
    output.write('\n')
    i = 0
    for sgrna in list_ordered:
        if i == nb_top : break
        i += 1
        hit = dic_hits[sgrna]
        output.write(sgrna+'\t')
        to_write = hit.write(genomes_in)
        output.write(to_write + '\n')
    output.close()

def isOnGene(coord, gene):
    start = int(re.search("[+-]\(([0-9]*),", coord).group(1))
    end = int(re.search(",([0-9]*)", coord).group(1))
    return gene.start <= start and gene.end >= end