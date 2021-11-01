"""
    Class defining a new Regulatory Interaction
"""
class RegulatoryInteraction:
    def __init__(self, id, tf_locus_tag, tg_locus_tag, regulatory_function, pubmedid_source):
        self.id = id
        self.tf_locus_tag = tf_locus_tag
        self.tg_locus_tag = tg_locus_tag
        self.regulatory_function = regulatory_function
        self.pubmedid_source = pubmedid_source
