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

class SeveralGenes(Exception):
    pass