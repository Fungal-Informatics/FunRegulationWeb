class Organism:
    def __init__(self, assembly_accession, order, genus, species, strain, is_model, cis_bp):
        self.assembly_accession = assembly_accession
        self.order = order
        self.genus = genus
        self.species = species
        self.strain = strain
        self.is_model = is_model
        self.cis_bp = cis_bp

"""
    Class defining a Gene
"""
class Gene:
    def __init__(self, organism_accession, locus_tag, symbol_gene, is_tf):
        self.organism_accession = organism_accession
        self.locus_tag = locus_tag
        self.symbol_gene = symbol_gene
        self.is_tf = is_tf
        
"""
    Class defining a Promoter
"""
class Promoter:
    def __init__(self, organism_accession, locus_tag, strand, source, start, stop, promoter_seq):
        self.organism_accession = organism_accession
        self.locus_tag = locus_tag
        self.strand = strand
        self.source = source
        self.start = start
        self.stop = stop
        self.promoter_seq = promoter_seq

"""
    Class defining a Protein
"""
class Protein:
    def __init__(self, organism_accession, locus_tag, id, product, interpro, pfam, go, gene3d, reactome, panther, uniprot, ec_number, cazy):
        self.organism_accession = organism_accession
        self.locus_tag = locus_tag
        self.id = id
        self.product = product
        self.interpro = interpro
        self.pfam = pfam
        self.go = go
        self.gene3d = gene3d
        self.reactome = reactome
        self.panther = panther
        self.uniprot = uniprot
        self.ec_number = ec_number
        self.cazy = cazy

"""
    Class defining an Orthology
"""
class Orthology:
    def __init__(self, model_assembly_accession, model_locus_tag, model_protein, target_assembly_accession, target_locus_tag, target_protein):
        self.model_assembly_accession = model_assembly_accession
        self.model_locus_tag = model_locus_tag
        self.model_protein = model_protein
        self.target_assembly_accession = target_assembly_accession
        self.target_locus_tag = target_locus_tag
        self.target_protein = target_protein
        
        
"""
    Class defining a PWM
"""
class Pwm:
    def __init__(self, id, organism_accession, locus_tag, motif_id, status, tf_family, motif_type, msource_author, msource, pubmedid):
        self.id = id
        self.organism_accession = organism_accession
        self.locus_tag = locus_tag
        self.motif_id = motif_id
        self.status = status
        self.tf_family = tf_family
        self.motif_type = motif_type
        self.msource_author = msource_author
        self.msource = msource
        self.pubmedid = pubmedid
        
"""
    Class defining a Model Regulatory Interaction
"""
class ModelRegulatory:
    def __init__(self, id, organism_accession, tf_locus_tag, tg_locus_tag, regulatory_function, evidence, experiment, experimental_condition, pubmedid, publication):
        self.id = id
        self.organism_accession = organism_accession
        self.tf_locus_tag = tf_locus_tag
        self.tg_locus_tag = tg_locus_tag
        self.regulatory_function = regulatory_function
        self.evidence = evidence
        self.experiment = experiment
        self.experimental_condition = experimental_condition
        self.pubmedid = pubmedid
        self.publication = publication

"""
    Class defining a new Regulatory Interaction
"""
class RegulatoryInteraction:
    def __init__(self, id, assembly_accession, tf_locus_tag, tg_locus_tag, regulatory_function, pubmedid_source):
        self.id = id
        self.assembly_accession = assembly_accession
        self.tf_locus_tag = tf_locus_tag
        self.tg_locus_tag = tg_locus_tag
        self.regulatory_function = regulatory_function
        self.pubmedid_source = pubmedid_source

"""
    Class defining a new TFBS prediction
"""
class TFBS:
    def __init__(self, id, regulatory_interaction_id, assembly_accession, tf_locus_tag, tg_locus_tag, pwm_id, strand, start, end, sequence, weight, pval, ln_pval, sig):
        self.id = id
        self.regulatory_interaction_id = regulatory_interaction_id
        self.assembly_accession = assembly_accession
        self.tf_locus_tag = tf_locus_tag
        self.tg_locus_tag = tg_locus_tag
        self.pwm_id = pwm_id
        self.strand = strand
        self.start = start
        self.end = end
        self.sequence = sequence
        self.weight = weight
        self.pval = pval
        self.ln_pval = ln_pval
        self.sig = sig