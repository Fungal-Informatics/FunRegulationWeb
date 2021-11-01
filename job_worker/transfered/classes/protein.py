"""
    Class defining a Promoter
"""
class Protein:
    def __init__(self, locus_tag, id, interpro, pfam, go, gene3d, reactome, panther, uniprot, kegg_enzyme, cazy, uniparc):
        self.locus_tag = locus_tag
        self.id = id
        self.interpro = interpro
        self.pfam = pfam
        self.go = go
        self.gene3d = gene3d
        self.reactome = reactome
        self.panther = panther
        self.uniprot = uniprot
        self.kegg_enzyme = kegg_enzyme
        self.cazy = cazy
        self.uniparc = uniparc
