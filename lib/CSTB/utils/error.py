import json
import traceback

class CouchNotFound(Exception):
    pass

class ConsistencyError(Exception):
    pass

class PingError(Exception):
    pass

class NoHomolog(Exception):
    pass

class NoBlastHit(Exception):
    pass

class FastaMetadataError(Exception):
    """Raise if fasta metadata has wrong format
    """
    pass

class ArgumentError(Exception):
    """Raise if error in command line argument
    """
    pass

def empty_exit(message, job_number):
    """Print json with emptySearch key and exit
    
    :param message: Message to display
    :type message: str
    """
    json_dic = {"emptySearch" : message, "tag":job_number}
    print(json.dumps(json_dic))
    exit()

def error_exit(message, job_number): 
    """Print json with error key, traceback error and exit
    
    :param message: Message to display
    :type message: str
    """
    json_dic = {"error" : message + f"\n Contact support with the job number {job_number} : cecile.hilpert@ibcp.fr", "tag" : job_number}
    print(json.dumps(json_dic)) #Need to be json dumped
    traceback.print_exc()
    exit()