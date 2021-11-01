"""
    Class defining a Model Regulatory Interaction
"""
class ModelRegulatory:
    def __init__(self, id, tf_locus_tag, tg_locus_tag, regulatory_function, evidence, experiment, experimental_condition, pubmedid, publication):
        self.id = id
        self.tf_locus_tag = tf_locus_tag
        self.tg_locus_tag = tg_locus_tag
        self.regulatory_function = regulatory_function
        self.evidence = evidence
        self.experiment = experiment
        self.experimental_condition = experimental_condition
        self.pubmedid = pubmedid
        self.publication = publication
