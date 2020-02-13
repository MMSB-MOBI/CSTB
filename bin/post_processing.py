import logging
logging.basicConfig(filename = "post_processing.log", level = logging.DEBUG, format='%(levelname)s\t%(message)s')
import argparse
from CSTB.crispr_result_manager import CrisprResultManager
import pycouch.wrapper as couch_wrapper
import json
import traceback
import CSTB.utils.error as error

'''TO DO
- Add verbosity parameter
- Change error document and dump a json
'''

def args_gestion():
    """
    Take and treat arguments that user gives in command line
    """
    # Argparsing
    parser = argparse.ArgumentParser(description="Post-processing results")
    parser.add_argument("--exclude", metavar="<str>",
                        help="The organisms to search exclusion from",
                        required=True)
    parser.add_argument("--include", metavar="<str>",
                        help="The organisms to search inclusion in.",
                        required=True)
    parser.add_argument("--couch_endpoint", metavar="<str>",
                        help="The end point of the taxon and tree database",
                        required=True)
    parser.add_argument("--taxon_db", metavar="<str>",
                        help="The name of the taxon database",
                        required=True)
    parser.add_argument("--genome_db", metavar="<str>",
                        help="The name of the genome database",
                        required=True)
    parser.add_argument("--set_compare", metavar="<path>", help="setCompare result file", required = True)
    parser.add_argument("--length", metavar="<int>", help = "sgRNA length (exclude pam)", required = True, type=int)
    parser.add_argument("--motif_broker_endpoint", metavar="<url>", help = "Motif broker endpoint", required = True)
    parser.add_argument("--tag", metavar = "<str>", help = "tag for outputs", required = True)
    parser.add_argument("--blast", metavar = "<file>", help = "Blast results (xml format) if specific gene")
    return parser.parse_args()

def error_exit(message): 
    json_dic = {"error" : message}
    print(json.dumps(json_dic)) #Need to be json dumped
    traceback.print_exc()
    exit()

def empty_exit(message):
    json_dic = {"emptySearch" : message}
    print(json.dumps(json_dic))
    exit()

def main():
    logging.info("== post_processing.py")
    global PARAM
    PARAM = args_gestion()

    #Get genomes uuid from args
    include = PARAM.include.split("&")
    exclude = PARAM.exclude.split("&") if PARAM.exclude else []

    logging.debug(PARAM.include)
    logging.debug(PARAM.exclude)

    logging.info(f"Include genomes : {include} ({len(include)})\nExclude genomes : {exclude} ({len(exclude)})\n")

    #Initialize pycouch wrapper
    logging.info("= Initialize pycouch wrapper")
    wrapper = couch_wrapper.Wrapper(PARAM.couch_endpoint)
    if not wrapper.couchPing():
        error_exit("Can't ping couch database")

    results = CrisprResultManager(wrapper, PARAM.taxon_db, PARAM.genome_db, PARAM.motif_broker_endpoint, PARAM.tag)

    logging.info("= Interrogate couchDB to retrieve taxon name")
    try:
        results.set_taxon_names(include, exclude)
    except:
        error_exit("Error while set taxon names")

    logging.info("= Parse setCompare")
    try:
        results.parse_set_compare(PARAM.set_compare, PARAM.length, 1000)
    except:
        error_exit("Error while parse setCompare")

    if not results.hits_collection:
        empty_exit("No hits")
    
    logging.info("= Search sgrna occurences in couchDB")
    try:
        results.search_occurences(include)
    except:
        error_exit("Error while search sgRNA occurences in couchDB")

    if PARAM.blast:
        logging.info("= Parse Blast")
        try: 
            results.parseBlast(PARAM.blast, 70, include)
        except error.NoBlastHit:
            empty_exit("No blast hit for your gene.")
        except error.NoHomolog :
            empty_exit(f"Some organisms don't have homolog gene.")
        except:
            error_exit("Error while parse blast")
        
        try:
            filtered_results = results.filterOnGeneOccurences()
        except:
            error_exit("Error while filter blast results")
    
    else:
        filtered_results = results

    logging.info(f"{len(filtered_results.hits_collection)}")

    logging.info("= Format results")
    try:
        blast = False
        if PARAM.blast:
            blast = True
        json_results = filtered_results.format_results(blast)
    except:
        error_exit("Error while format results")
    
    print(json.dumps(json_results))


if __name__ == '__main__':
    main()
    