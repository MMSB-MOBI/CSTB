import json

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

def empty_exit(message):
    json_dic = {"emptySearch" : message}
    print(json.dumps(json_dic))
    exit()

def error_exit(message): 
    json_dic = {"error" : message}
    print(json.dumps(json_dic)) #Need to be json dumped
    traceback.print_exc()
    exit()