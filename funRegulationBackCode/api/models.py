import enum

from django.db import models
from django.db.models import Q, Case, When, IntegerField, Sum
from django.db.models.signals import post_save
from django.contrib.auth.models import User
from django_celery_results.models import TaskResult
from django.dispatch import receiver

from funRegulationTool import task_utils
from root.utils import file_utils, task_status_utils, funregulationtools_utils

class Gene(models.Model):
    organism_accession = models.ForeignKey('Organism', models.DO_NOTHING, db_column='organism_accession', blank=True, null=True, related_name='org_accession')
    locus_tag = models.CharField(primary_key=True, max_length=255, db_column='locus_tag')
    symbol_gene = models.CharField(max_length=255, blank=True, null=True)
    is_tf = models.BooleanField(blank=True, null=True)

    class Meta:
        constraints = [models.UniqueConstraint(fields=['organism_accession','locus_tag'], name='unique_gene')]
        db_table = 'gene'

class ModelRegulatory(models.Model):
    organism_accession = models.ForeignKey('Organism', models.DO_NOTHING, db_column='organism_accession', blank=True, null=True, related_name='mr_gene_organism_accession')
    tf_locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, db_column='tf_locus_tag', blank=True, null=True, related_name='mr_gene_tf_locus_tag')
    tg_locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, db_column='tg_locus_tag', blank=True, null=True, related_name='mr_gene_tg_locus_tag')
    regulatory_function = models.CharField(max_length=255, blank=True, null=True)
    evidence = models.CharField(max_length=255, blank=True, null=True)
    experiment = models.CharField(max_length=255, blank=True, null=True)
    experimental_condition = models.CharField(max_length=255, blank=True, null=True)
    pubmedid = models.CharField(max_length=255, blank=True, null=True)
    publication = models.CharField(max_length=255, blank=True, null=True)

    class Meta:
        constraints = [models.UniqueConstraint(fields=['organism_accession','tf_locus_tag','tg_locus_tag'], 
                                               name='unique_model_regulatory')]
        db_table = 'model_regulatory'

class Organism(models.Model):
    accession = models.CharField(primary_key=True, max_length=255)
    order = models.CharField(max_length=255, blank=True, null=True)
    genus = models.CharField(max_length=255, blank=True, null=True)
    species = models.CharField(max_length=255, blank=True, null=True)
    strain = models.CharField(max_length=255, blank=True, null=True)
    is_model = models.BooleanField(blank=True, null=True)
    cis_bp = models.BooleanField(blank=True, null=True)

    def __str__(self):
        return self.accession

class Protein(models.Model):
    organism_accession = models.ForeignKey(Organism, models.DO_NOTHING, db_column='organism_accession', blank=True, null=True, related_name='prot_gene_organism_accession')
    locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, db_column='locus_tag', blank=True, null=True, related_name='prot_gene_locus_tag')
    id = models.CharField(primary_key=True, max_length=255)
    product = models.CharField(max_length=255, blank=True, null=True)
    interpro = models.CharField(max_length=255, blank=True, null=True)
    pfam = models.CharField(max_length=255, blank=True, null=True)
    go = models.CharField(max_length=255, blank=True, null=True)
    gene3d = models.CharField(max_length=255, blank=True, null=True)
    reactome = models.CharField(max_length=255, blank=True, null=True)
    panther = models.CharField(max_length=255, blank=True, null=True)
    uniprot = models.CharField(max_length=255, blank=True, null=True)
    ec_number = models.CharField(max_length=255, blank=True, null=True)
    cazy = models.CharField(max_length=255, blank=True, null=True)

    class Meta:
        constraints = [models.UniqueConstraint(fields=['organism_accession','locus_tag', 'id'], name='unique_protein')]
        db_table = 'protein'

class Orthology(models.Model):
    model_organism_accession = models.ForeignKey(Organism, models.DO_NOTHING,db_column='model_organism_accession', null=True,related_name='ort_model_organism_accession', default='')
    model_locus_tag = models.ForeignKey(Gene,models.DO_NOTHING, db_column='model_locus_tag', null=True,related_name='ort_model_locus_tag', default='')
    model_protein = models.ForeignKey(Protein, models.DO_NOTHING, db_column='model_protein', related_name='model_protein', default='')
    target_organism_accession = models.ForeignKey(Organism, models.DO_NOTHING, db_column='target_organism_accession', null=True,related_name='ort_target_organism_accession', default='')
    target_locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, db_column='target_locus_tag', null=True,related_name='ort_target_locus_tag', default='')
    target_protein = models.ForeignKey(Protein, models.DO_NOTHING, db_column='target_protein', related_name='target_protein' , default='')

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['id','model_organism_accession','model_locus_tag',
                                            'target_organism_accession','target_locus_tag',
                                            'model_protein','target_protein'], 
                                            name='unique_orthology')]
        db_table = 'orthology'

class Promoter(models.Model):
    organism_accession = models.ForeignKey(Organism, models.DO_NOTHING,db_column='organism_accession', related_name='prot_organism_accession', blank=True)
    locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, db_column='locus_tag', blank=True, primary_key=True)
    strand = models.CharField(max_length=3, blank=True, null=True)
    source = models.CharField(max_length=255, blank=True, null=True)
    start = models.IntegerField(blank=True, null=True)
    stop = models.IntegerField(blank=True, null=True)
    promoter_seq = models.TextField(blank=True, null=True)

    class Meta:
        constraints = [models.UniqueConstraint(fields=['organism_accession','locus_tag'], name='unique_promoter')]
        db_table = 'promoter'

class Pwm(models.Model):
    organism_accession = models.ForeignKey(Organism, models.DO_NOTHING,db_column='organism_accession', default='')
    locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, db_column='locus_tag', blank=True, null=True)
    motif_id = models.CharField(max_length=255, blank=True, null=True)
    status = models.CharField(max_length=1, blank=True, null=True)
    tf_family = models.CharField(max_length=255, blank=True, null=True)
    motif_type = models.CharField(max_length=255, blank=True, null=True)
    msource_author = models.CharField(max_length=255, blank=True, null=True)
    msource = models.CharField(max_length=255, blank=True, null=True)
    pubmedid = models.CharField(max_length=255, blank=True, null=True)

    class Meta:
        constraints = [models.UniqueConstraint(fields=['organism_accession','locus_tag'], name='unique_pwm')]
        db_table = 'pwm'

class RegulatoryInteraction(models.Model):
    organism_accession = models.ForeignKey(Organism, models.DO_NOTHING,db_column='organism_accession', default='')
    tf_locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, db_column='tf_locus_tag', blank=True, null=True, related_name='reg_tf_locus_tag')
    tg_locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, db_column='tg_locus_tag', blank=True, null=True, related_name='reg_tg_locus_tag')
    regulatory_function = models.CharField(max_length=255, blank=True, null=True)
    pubmedid_source = models.CharField(max_length=255, blank=True, null=True)

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['id','organism_accession','tf_locus_tag','tg_locus_tag'], 
                                            name='unique_regulatory_interaction')]
        db_table = 'regulatory_interaction'

class Tfbs(models.Model):
    regulatory_interaction = models.ForeignKey(RegulatoryInteraction, models.DO_NOTHING, blank=True, null=True)
    organism_accession = models.ForeignKey(Organism, models.DO_NOTHING, db_column='organism_accession', blank=True, null=True, related_name='tfbs_gene_organism_accession')
    tf_locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, db_column='tf_locus_tag', blank=True, null=True, related_name='tfbs_gene_tf_locus_tag')
    tg_locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, db_column='tg_locus_tag', blank=True, null=True, related_name='tfbs_gene_tg_locus_tag')
    pwm = models.ForeignKey(Pwm, models.DO_NOTHING, blank=True, null=True)
    strand = models.CharField(max_length=255, blank=True, null=True)
    start = models.CharField(max_length=255, blank=True, null=True)
    end = models.CharField(max_length=255, blank=True, null=True)
    sequence = models.CharField(max_length=255, blank=True, null=True)
    weight = models.CharField(max_length=255, blank=True, null=True)
    pval = models.CharField(max_length=255, blank=True, null=True)
    ln_pval = models.CharField(max_length=255, blank=True, null=True)
    sig = models.CharField(max_length=255, blank=True, null=True)

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['regulatory_interaction','organism_accession',
                                            'tf_locus_tag','tg_locus_tag'], 
                                            name='unique_tfbs')]
        db_table = 'tfbs'
  
class ProjectAnalysisRegistry(models.Model):
    active = models.BooleanField(default=True)
    date_created = models.DateTimeField(auto_now_add=True)
    created_by = models.ForeignKey(User, related_name='analysis_registries_created', on_delete=models.PROTECT)
    date_inactive = models.DateTimeField(null=True, blank=True)
    organism_accession = models.CharField(max_length=100, null=False, default="vazio")
    download_organism = models.BooleanField(default=False)
    proteinortho_analyse = models.BooleanField(default=True)
    rsat_analyse = models.BooleanField(default=False)
    proteinortho_analysed = models.BooleanField(default=False)
    proteinortho_error = models.IntegerField(null=True, blank=True)
    task_proteinortho = models.OneToOneField(TaskResult, null=True, blank=True,
                                         related_name='proteinortho_analysis_registry_item',
                                         on_delete=models.SET_NULL)
    rsat_analysed = models.BooleanField(default=False)
    rsat_error = models.IntegerField(null=True, blank=True)
    task_rsat = models.OneToOneField(TaskResult, null=True, blank=True,
                                         related_name='rsat_analysis_registry_item',
                                         on_delete=models.SET_NULL)
    download_completed = models.BooleanField(default=False)
    task_download_organism = models.OneToOneField(TaskResult, null=True, blank=True,on_delete=models.SET_NULL)
    task_create_graph = models.OneToOneField(TaskResult, null=True, blank=True, 
                                             related_name='create_graph_analysis_registry_item',
                                             on_delete=models.SET_NULL)
    create_graph_error = models.IntegerField(null=True, blank=True)
    task_calculate_degree = models.OneToOneField(TaskResult, null=True, blank=True, 
                                                 related_name='calculate_degree_analysis_registry_item',
                                                 on_delete=models.SET_NULL)
    degree_error = models.IntegerField(null=True, blank=True)
    task_calculate_closeness = models.OneToOneField(TaskResult, null=True, blank=True,
                                                    related_name='calculate_closeness_analysis_registry_item',
                                                    on_delete=models.SET_NULL)
    closeness_error = models.IntegerField(null=True, blank=True)
    task_calculate_betweenness = models.OneToOneField(TaskResult, null=True, blank=True,
                                                      related_name='calculate_betweenness_analysis_registry_item',
                                                      on_delete=models.SET_NULL)
    betweenness_error = models.IntegerField(null=True, blank=True)
    task_calculate_eigenvector = models.OneToOneField(TaskResult, null=True, blank=True,
                                                      related_name='calculate_eigenvector_analysis_registry_item',
                                                      on_delete=models.SET_NULL)
    eigenvector_error = models.IntegerField(null=True, blank=True)
    task_calculate_harmonic = models.OneToOneField(TaskResult, null=True, blank=True,
                                                   related_name='calculate_harmonic_analysis_registry_item',
                                                   on_delete=models.SET_NULL)
    harmonic_error = models.IntegerField(null=True, blank=True)
    task = models.OneToOneField(TaskResult, null=True, blank=True,
                                related_name='analysis_registry', on_delete=models.SET_NULL)

    def __str__(self):
        return str(self.pk)

class CalculateCentralityRegistry(models.Model):
    locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, primary_key=True, db_column='locus_tag', related_name='calculate_centrality_gene_locus_tag')
    degree = models.FloatField(null=True, blank=True)
    closeness = models.FloatField(null=True, blank=True)
    betweenness = models.FloatField(null=True, blank=True)
    eigenvector = models.FloatField(null=True, blank=True)
    harmonic = models.FloatField(null=True, blank=True)

    def __str__(self):
        return str(self.pk)

class ProteinOrthoErrorType(enum.Enum):
    COMMAND_ERROR = 1

class RsatErrorType(enum.Enum):
    COMMAND_ERROR = 1

class CalculateCentralityErrorType(enum.Enum):
    COMMAND_ERROR = 1

class SystemPreferenceType(enum.Enum):
    STRING = 1
    INTEGER = 2
    JSON = 3
    BOOLEAN = 4
    FLOAT = 5

class Profile(models.Model):
    user = models.OneToOneField(User, related_name='profile', on_delete=models.PROTECT)
    organization = models.CharField(max_length=200, null=True, blank=True)
    country = models.CharField(max_length=100, null=True, blank=True)
    BrazilianState = models.CharField(max_length=2, null=True, blank=True)
    is_admin = models.BooleanField(default=False)
    account_confirmation = models.BooleanField(default=False)
    account_confirmation_date = models.DateTimeField(null=True, blank=True)
    first_login_date = models.DateTimeField(null=True, blank=True)

    def __str__(self):
        return self.user.username
