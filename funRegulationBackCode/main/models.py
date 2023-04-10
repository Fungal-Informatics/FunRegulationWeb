# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#   * Rearrange models' order
#   * Make sure each model has one field with primary_key=True
#   * Make sure each ForeignKey and OneToOneField has `on_delete` set to the desired behavior
#   * Remove `managed = False` lines if you wish to allow Django to create, modify, and delete the table
# Feel free to rename the models, but don't rename db_table values or field names.
from django.db import models

class Gene(models.Model):
    organism_accession = models.ForeignKey('Organism', models.DO_NOTHING, db_column='organism_accession', blank=True, null=True)
    locus_tag = models.CharField(primary_key=True, max_length=-1)
    symbol_gene = models.CharField(max_length=-1, blank=True, null=True)
    description = models.CharField(max_length=-1, blank=True, null=True)
    is_tf = models.BooleanField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'gene'


class ModelRegulatory(models.Model):
    tf_locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, db_column='tf_locus_tag', blank=True, null=True)
    tg_locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, db_column='tg_locus_tag', blank=True, null=True)
    regulatory_function = models.CharField(max_length=-1, blank=True, null=True)
    evidence = models.CharField(max_length=-1, blank=True, null=True)
    experiment = models.CharField(max_length=-1, blank=True, null=True)
    experimental_condition = models.CharField(max_length=-1, blank=True, null=True)
    pubmedid = models.CharField(max_length=-1, blank=True, null=True)
    publication = models.CharField(max_length=-1, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'model_regulatory'


class Organism(models.Model):
    accession = models.CharField(primary_key=True, max_length=-1)
    order = models.CharField(max_length=-1, blank=True, null=True)
    genus = models.CharField(max_length=-1, blank=True, null=True)
    species = models.CharField(max_length=-1, blank=True, null=True)
    strain = models.CharField(max_length=-1, blank=True, null=True)
    is_model = models.BooleanField(blank=True, null=True)
    cis_bp = models.BooleanField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'organism'


class Orthology(models.Model):
    model_protein = models.ForeignKey('Protein', models.DO_NOTHING, db_column='model_protein', blank=True, null=True)
    target_protein = models.ForeignKey('Protein', models.DO_NOTHING, db_column='target_protein', blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'orthology'


class Promoter(models.Model):
    locus_tag = models.OneToOneField(Gene, models.DO_NOTHING, db_column='locus_tag', blank=True, null=True)
    strand = models.CharField(max_length=1, blank=True, null=True)
    source = models.CharField(max_length=-1, blank=True, null=True)
    start = models.IntegerField(blank=True, null=True)
    stop = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'promoter'


class Protein(models.Model):
    locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, db_column='locus_tag', blank=True, null=True)
    id = models.CharField(primary_key=True, max_length=-1)
    interpro = models.CharField(max_length=-1, blank=True, null=True)
    pfam = models.CharField(max_length=-1, blank=True, null=True)
    go = models.CharField(max_length=-1, blank=True, null=True)
    gene3d = models.CharField(max_length=-1, blank=True, null=True)
    reactome = models.CharField(max_length=-1, blank=True, null=True)
    panther = models.CharField(max_length=-1, blank=True, null=True)
    uniprot = models.CharField(max_length=-1, blank=True, null=True)
    kegg_enzyme = models.CharField(max_length=-1, blank=True, null=True)
    cazy = models.CharField(max_length=-1, blank=True, null=True)
    uniparc = models.CharField(max_length=-1, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'protein'


class Pwm(models.Model):
    locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, db_column='locus_tag', blank=True, null=True)
    motif_id = models.CharField(max_length=-1, blank=True, null=True)
    status = models.CharField(max_length=1, blank=True, null=True)
    tf_family = models.CharField(max_length=-1, blank=True, null=True)
    motif_type = models.CharField(max_length=-1, blank=True, null=True)
    msource_author = models.CharField(max_length=-1, blank=True, null=True)
    msource = models.CharField(max_length=-1, blank=True, null=True)
    pubmedid = models.CharField(max_length=-1, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'pwm'


class RegulatoryInteraction(models.Model):
    tf_locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, db_column='tf_locus_tag', blank=True, null=True)
    tg_locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, db_column='tg_locus_tag', blank=True, null=True)
    regulatory_function = models.CharField(max_length=-1, blank=True, null=True)
    pubmedid_source = models.CharField(max_length=-1, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'regulatory_interaction'


class Tfbs(models.Model):
    regulatory_interaction = models.ForeignKey(RegulatoryInteraction, models.DO_NOTHING, blank=True, null=True)
    pwm = models.ForeignKey(Pwm, models.DO_NOTHING, blank=True, null=True)
    strand = models.CharField(max_length=-1, blank=True, null=True)
    start = models.CharField(max_length=-1, blank=True, null=True)
    end = models.CharField(max_length=-1, blank=True, null=True)
    sequence = models.CharField(max_length=-1, blank=True, null=True)
    weight = models.CharField(max_length=-1, blank=True, null=True)
    pval = models.CharField(max_length=-1, blank=True, null=True)
    ln_pval = models.CharField(max_length=-1, blank=True, null=True)
    sig = models.CharField(max_length=-1, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'tfbs'
