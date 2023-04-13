import enum

from django.db import models
from django.db.models import Q, Case, When, IntegerField, Sum
from django.db.models.signals import post_save
from django.contrib.auth.models import User
from django_celery_results.models import TaskResult
from django.dispatch import receiver

from funRegulationTool import task_utils

from root.utils import file_utils, task_status_utils


class Gene(models.Model):
    organism_accession = models.ForeignKey('Organism', models.DO_NOTHING, db_column='organism_accession', blank=True, null=True, related_name='org_accession')
    locus_tag = models.CharField(primary_key=True, max_length=255, db_column='locus_tag')
    symbol_gene = models.CharField(max_length=255, blank=True, null=True)
    description = models.CharField(max_length=500, blank=True, null=True)
    is_tf = models.BooleanField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'gene'


class ModelRegulatory(models.Model):
    tf_locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, db_column='tf_locus_tag', blank=True, null=True, related_name='gene_tf_locus_tag')
    tg_locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, db_column='tg_locus_tag', blank=True, null=True, related_name='gene_tg_locus_tag')
    regulatory_function = models.CharField(max_length=255, blank=True, null=True)
    evidence = models.CharField(max_length=255, blank=True, null=True)
    experiment = models.CharField(max_length=255, blank=True, null=True)
    experimental_condition = models.CharField(max_length=255, blank=True, null=True)
    pubmedid = models.CharField(max_length=255, blank=True, null=True)
    publication = models.CharField(max_length=255, blank=True, null=True)

    class Meta:
        managed = False
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
    
    class Meta:
        managed = False
        db_table = 'organism'


class Orthology(models.Model):
    model_protein = models.ForeignKey('Protein', models.DO_NOTHING, db_column='model_protein', blank=True, null=True, related_name='protein_model_protein')
    target_protein = models.ForeignKey('Protein', models.DO_NOTHING, db_column='target_protein', blank=True, null=True, related_name='protein_target_protein')

    class Meta:
        managed = False
        db_table = 'orthology'


class Promoter(models.Model):
    locus_tag = models.OneToOneField(Gene, models.DO_NOTHING, db_column='locus_tag', blank=True, null=True)
    strand = models.CharField(max_length=1, blank=True, null=True)
    source = models.CharField(max_length=255, blank=True, null=True)
    start = models.IntegerField(blank=True, null=True)
    stop = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'promoter'


class Protein(models.Model):
    locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, db_column='locus_tag', blank=True, null=True, related_name='gene_locus_tag')
    id = models.CharField(primary_key=True, max_length=255)
    interpro = models.CharField(max_length=255, blank=True, null=True)
    pfam = models.CharField(max_length=255, blank=True, null=True)
    go = models.CharField(max_length=255, blank=True, null=True)
    gene3d = models.CharField(max_length=255, blank=True, null=True)
    reactome = models.CharField(max_length=255, blank=True, null=True)
    panther = models.CharField(max_length=255, blank=True, null=True)
    uniprot = models.CharField(max_length=255, blank=True, null=True)
    kegg_enzyme = models.CharField(max_length=255, blank=True, null=True)
    cazy = models.CharField(max_length=255, blank=True, null=True)
    uniparc = models.CharField(max_length=255, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'protein'


class Pwm(models.Model):
    id = models.CharField(primary_key=True, max_length=255)
    locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, db_column='locus_tag', blank=True, null=True)
    motif_id = models.CharField(max_length=255, blank=True, null=True)
    status = models.CharField(max_length=1, blank=True, null=True)
    tf_family = models.CharField(max_length=255, blank=True, null=True)
    motif_type = models.CharField(max_length=255, blank=True, null=True)
    msource_author = models.CharField(max_length=255, blank=True, null=True)
    msource = models.CharField(max_length=255, blank=True, null=True)
    pubmedid = models.CharField(max_length=255, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'pwm'


class RegulatoryInteraction(models.Model):
    tf_locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, db_column='tf_locus_tag', blank=True, null=True, related_name='tf_locus_tag')
    tg_locus_tag = models.ForeignKey(Gene, models.DO_NOTHING, db_column='tg_locus_tag', blank=True, null=True, related_name='tg_locus_tag')
    regulatory_function = models.CharField(max_length=255, blank=True, null=True)
    pubmedid_source = models.CharField(max_length=255, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'regulatory_interaction'


class Tfbs(models.Model):
    regulatory_interaction = models.ForeignKey(RegulatoryInteraction, models.DO_NOTHING, blank=True, null=True)
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
        managed = False
        db_table = 'tfbs'

class Project(models.Model):
    name = models.CharField(max_length=100)
    strain = models.CharField(max_length=10)
    description = models.CharField(max_length=1000, null=True, blank=True)
    locus_prefix = models.CharField(max_length=10)
    transcript_suffix_type = models.IntegerField()
    contig_prefix = models.CharField(max_length=20)
    date_created = models.DateTimeField(auto_now_add=True)
    created_by = models.ForeignKey(User, related_name='projects', on_delete=models.PROTECT)
    removed = models.BooleanField(default=False)
    date_removed = models.DateTimeField(null=True, blank=True)
    removed_by = models.ForeignKey(User, null=True, blank=True,
                                   related_name='projects_removed', on_delete=models.PROTECT)    

    def __str__(self):
        return self.name

    def get_qtd_pending_analysis(self, distinct_feature=False):
        fields = task_status_utils.get_analysis_fields()
        #features = GeneFeature.objects.filter(gentool_utils.get_valid_feature_filters(self))
        features = task_status_utils.get_feature_status_statistic(features)\
            .filter(task_status_utils.has_analysis_in_progress_filter())
        if distinct_feature:
            return features.count()

        aggregations = {}
        for field in fields:
            features = features.annotate(**{
                'qtd_%s' % field: Case(When(**{'last_%s_analysis_in_progress' % field: True}, then=1),
                                       default=0, output_field=IntegerField())
            })
            aggregations['total_%s' % field] = Sum('qtd_%s' % field)
        totals = features.aggregate(**aggregations)
        total = 0
        for field in fields:
            t = totals.get('total_%s' % field, 0)
            if t is None:
                t = 0
            total += t
        return total

    def has_pending_analysis(self):
        return self.get_qtd_pending_analysis() > 0

    def has_pending_imports(self):
        return self.gene_imports.filter(~Q(task__status__in=task_utils.get_task_status_finished())).count() > 0

    def has_pending_exports(self):
        return self.export_registries.filter(~Q(task__status__in=task_utils.get_task_status_finished())).count() > 0
    
class ProjectAnalysisRegistry(models.Model):
    project = models.ForeignKey(Project, related_name='analysis_registries', on_delete=models.PROTECT)
    date_created = models.DateTimeField(auto_now_add=True)
    created_by = models.ForeignKey(User, related_name='analysis_registries_created', on_delete=models.PROTECT)
    date_inactive = models.DateTimeField(null=True, blank=True)
    pfam_analyse = models.BooleanField(default=False)
    interpro_analyse = models.BooleanField(default=False)
    task = models.OneToOneField(TaskResult, null=True, blank=True,
                                related_name='analysis_registry', on_delete=models.SET_NULL)

    def __str__(self):
        return str(self.pk)

class ProjectAnalysisRegistryItem(models.Model):
    registry = models.ForeignKey(ProjectAnalysisRegistry, related_name='items', on_delete=models.PROTECT)
    #feature = models.ForeignKey(GeneFeature, related_name='analysis_registries_items', on_delete=models.PROTECT)
    active = models.BooleanField(default=True)
    date_inactive = models.DateTimeField(null=True, blank=True)
    #pfam_analysed = models.BooleanField(default=False)
    #pfam_error = models.IntegerField(null=True, blank=True)
    proteinortho_analysed = models.BooleanField(default=False)
    proteinortho_error = models.IntegerField(null=True, blank=True)
    #task_pfam = models.OneToOneField(TaskResult, null=True, blank=True,
    #                                 related_name='pfam_analysis_registry_item', on_delete=models.SET_NULL)
    task_proteinortho = models.OneToOneField(TaskResult, null=True, blank=True,
                                         related_name='proteinortho_analysis_registry_item',
                                         on_delete=models.SET_NULL)

    def __str__(self):
        return str(self.pk)

class ProteinOrthoErrorType(enum.Enum):
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
    is_admin = models.BooleanField(default=False)
    account_confirmation = models.BooleanField(default=False)
    account_confirmation_date = models.DateTimeField(null=True, blank=True)
    first_login_date = models.DateTimeField(null=True, blank=True)

    def __str__(self):
        return self.user.username


@receiver(post_save, sender=User)
def create_user_profile(sender, instance, created, **kwargs):
    if created:
        Profile.objects.create(user=instance)


@receiver(post_save, sender=User)
def save_user_profile(sender, instance, **kwargs):
    instance.profile.save()
