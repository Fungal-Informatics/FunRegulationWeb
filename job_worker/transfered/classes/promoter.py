"""
    Class defining a Promoter
"""
class Promoter:
    def __init__(self, locus_tag, strand, source, start, stop):
        self.locus_tag = locus_tag
        self.strand = strand
        self.source = source
        self.start = start
        self.stop = stop
