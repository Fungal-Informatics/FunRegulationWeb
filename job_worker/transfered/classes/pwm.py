"""
    Class defining a PWM
"""
class Pwm:
    def __init__(self, id, locus_tag, motif_id, status, tf_family, motif_type, msource_author, msource, pubmedid):
        self.id = id
        self.locus_tag = locus_tag
        self.motif_id = motif_id
        self.status = status
        self.tf_family = tf_family
        self.motif_type = motif_type
        self.msource_author = msource_author
        self.msource = msource
        self.pubmedid = pubmedid