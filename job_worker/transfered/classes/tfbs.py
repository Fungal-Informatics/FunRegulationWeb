"""
    Class defining a new TFBS prediction
"""
class TFBS:
    def __init__(self, id, regulatory_interaction_id, pwm_id, strand, start, end, sequence, weight, pval, ln_pval, sig):
        self.id = id
        self.regulatory_interaction_id = regulatory_interaction_id
        self.pwm_id = pwm_id
        self.strand = strand
        self.start = start
        self.end = end
        self.sequence = sequence
        self.weight = weight
        self.pval = pval
        self.ln_pval = ln_pval
        self.sig = sig