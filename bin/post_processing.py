import logging
logging.basicConfig(filename = "post_processing.log", level = logging.DEBUG, format='%(levelname)s\t%(message)s')
import argparse
from CSTB.crispr_result import CrisprResultManagerOld
import pycouch.wrapper as couch_wrapper
import json
import traceback

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
    return parser.parse_args()

def error_exit(message): #WHERE TO PUT THIS ?
    print({"emptySearch" : message})
    traceback.print_exc()
    exit()

if __name__ == '__main__':
    logging.info("== post_processing.py")
    global PARAM
    PARAM = args_gestion()

    #Get genomes uuid from args
    include = PARAM.include.split("&")
    exclude = PARAM.exclude.split("&")

    #Initialize pycouch wrapper
    logging.info("= Initialize pycouch wrapper")
    wrapper = couch_wrapper.Wrapper(PARAM.couch_endpoint)
    if not wrapper.couchPing():
        print({"emptySearch": "Can't ping couch database"})

    results = CrisprResultManagerOld(wrapper, PARAM.taxon_db, PARAM.genome_db, PARAM.motif_broker_endpoint, PARAM.tag)

    logging.info("= Interrogate couchDB to retrieve taxon name")
    try:
        results.set_taxon_names(include, exclude)
    except:
        error_exit("Error while set taxon names")

    logging.info("= Parse setCompare")
    try:
        results.parse_set_compare(PARAM.set_compare, PARAM.length, 5000)
    except:
        error_exit("Error while parse setCompare")
    
    logging.info("= Search sgrna occurences in couchDB")
    try:
        results.search_occurences(include)
    except:
        error_exit("Error while search occurences in couchDB")

    logging.info("= Format results")
    try:
        json_results = results.format_results()
    except:
        error_exit("Error while format results")
    
    print(json.dumps(json_results))

