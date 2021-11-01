"""
    Class defining an Organism
"""
class Organism:
    def __init__(self, order, genus, species, strain, taxon, assembly_name, assembly_accession, is_model, cis_bp):
        self.order = order
        self.genus = genus
        self.species = species
        self.strain = strain
        self.taxon = taxon
        self.assembly_name = assembly_name
        self.assembly_accession = assembly_accession
        self.is_model = is_model
        self.cis_bp = cis_bp
