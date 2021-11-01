"""
    Class defining a Gene
"""
class Gene:
    def __init__(self, organism, locus_tag, symbol_gene, description, is_tf):
        self.organism = organism
        self.locus_tag = locus_tag
        self.symbol_gene = symbol_gene
        self.description = description
        self.is_tf = is_tf
        