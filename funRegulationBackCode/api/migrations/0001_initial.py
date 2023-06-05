# Generated by Django 4.2 on 2023-04-16 23:04

from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        ('django_celery_results', '0003_auto_20181106_1101'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='Gene',
            fields=[
                ('locus_tag', models.CharField(db_column='locus_tag', max_length=255, primary_key=True, serialize=False)),
                ('symbol_gene', models.CharField(blank=True, max_length=255, null=True)),
                ('description', models.CharField(blank=True, max_length=500, null=True)),
                ('is_tf', models.BooleanField(blank=True, null=True)),
            ],
            options={
                'db_table': 'gene',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='ModelRegulatory',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('regulatory_function', models.CharField(blank=True, max_length=255, null=True)),
                ('evidence', models.CharField(blank=True, max_length=255, null=True)),
                ('experiment', models.CharField(blank=True, max_length=255, null=True)),
                ('experimental_condition', models.CharField(blank=True, max_length=255, null=True)),
                ('pubmedid', models.CharField(blank=True, max_length=255, null=True)),
                ('publication', models.CharField(blank=True, max_length=255, null=True)),
            ],
            options={
                'db_table': 'model_regulatory',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='Organism',
            fields=[
                ('accession', models.CharField(max_length=255, primary_key=True, serialize=False)),
                ('order', models.CharField(blank=True, max_length=255, null=True)),
                ('genus', models.CharField(blank=True, max_length=255, null=True)),
                ('species', models.CharField(blank=True, max_length=255, null=True)),
                ('strain', models.CharField(blank=True, max_length=255, null=True)),
                ('is_model', models.BooleanField(blank=True, null=True)),
                ('cis_bp', models.BooleanField(blank=True, null=True)),
            ],
            options={
                'db_table': 'organism',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='Orthology',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
            ],
            options={
                'db_table': 'orthology',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='Promoter',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('strand', models.CharField(blank=True, max_length=1, null=True)),
                ('source', models.CharField(blank=True, max_length=255, null=True)),
                ('start', models.IntegerField(blank=True, null=True)),
                ('stop', models.IntegerField(blank=True, null=True)),
            ],
            options={
                'db_table': 'promoter',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='Protein',
            fields=[
                ('id', models.CharField(max_length=255, primary_key=True, serialize=False)),
                ('interpro', models.CharField(blank=True, max_length=255, null=True)),
                ('pfam', models.CharField(blank=True, max_length=255, null=True)),
                ('go', models.CharField(blank=True, max_length=255, null=True)),
                ('gene3d', models.CharField(blank=True, max_length=255, null=True)),
                ('reactome', models.CharField(blank=True, max_length=255, null=True)),
                ('panther', models.CharField(blank=True, max_length=255, null=True)),
                ('uniprot', models.CharField(blank=True, max_length=255, null=True)),
                ('kegg_enzyme', models.CharField(blank=True, max_length=255, null=True)),
                ('cazy', models.CharField(blank=True, max_length=255, null=True)),
                ('uniparc', models.CharField(blank=True, max_length=255, null=True)),
            ],
            options={
                'db_table': 'protein',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='Pwm',
            fields=[
                ('id', models.CharField(max_length=255, primary_key=True, serialize=False)),
                ('motif_id', models.CharField(blank=True, max_length=255, null=True)),
                ('status', models.CharField(blank=True, max_length=1, null=True)),
                ('tf_family', models.CharField(blank=True, max_length=255, null=True)),
                ('motif_type', models.CharField(blank=True, max_length=255, null=True)),
                ('msource_author', models.CharField(blank=True, max_length=255, null=True)),
                ('msource', models.CharField(blank=True, max_length=255, null=True)),
                ('pubmedid', models.CharField(blank=True, max_length=255, null=True)),
            ],
            options={
                'db_table': 'pwm',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='RegulatoryInteraction',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('regulatory_function', models.CharField(blank=True, max_length=255, null=True)),
                ('pubmedid_source', models.CharField(blank=True, max_length=255, null=True)),
            ],
            options={
                'db_table': 'regulatory_interaction',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='Tfbs',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('strand', models.CharField(blank=True, max_length=255, null=True)),
                ('start', models.CharField(blank=True, max_length=255, null=True)),
                ('end', models.CharField(blank=True, max_length=255, null=True)),
                ('sequence', models.CharField(blank=True, max_length=255, null=True)),
                ('weight', models.CharField(blank=True, max_length=255, null=True)),
                ('pval', models.CharField(blank=True, max_length=255, null=True)),
                ('ln_pval', models.CharField(blank=True, max_length=255, null=True)),
                ('sig', models.CharField(blank=True, max_length=255, null=True)),
            ],
            options={
                'db_table': 'tfbs',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='Project',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=100)),
                ('strain', models.CharField(max_length=10)),
                ('description', models.CharField(blank=True, max_length=1000, null=True)),
                ('locus_prefix', models.CharField(max_length=10)),
                ('transcript_suffix_type', models.IntegerField()),
                ('contig_prefix', models.CharField(max_length=20)),
                ('date_created', models.DateTimeField(auto_now_add=True)),
                ('removed', models.BooleanField(default=False)),
                ('date_removed', models.DateTimeField(blank=True, null=True)),
                ('created_by', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, related_name='projects', to=settings.AUTH_USER_MODEL)),
                ('removed_by', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.PROTECT, related_name='projects_removed', to=settings.AUTH_USER_MODEL)),
            ],
        ),
        migrations.CreateModel(
            name='ProjectAnalysisRegistry',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('date_created', models.DateTimeField(auto_now_add=True)),
                ('date_inactive', models.DateTimeField(blank=True, null=True)),
                ('pfam_analyse', models.BooleanField(default=False)),
                ('interpro_analyse', models.BooleanField(default=False)),
                ('created_by', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, related_name='analysis_registries_created', to=settings.AUTH_USER_MODEL)),
                ('project', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, related_name='analysis_registries', to='api.project')),
                ('task', models.OneToOneField(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='analysis_registry', to='django_celery_results.taskresult')),
            ],
        ),
        migrations.CreateModel(
            name='ProjectAnalysisRegistryItem',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('active', models.BooleanField(default=True)),
                ('date_inactive', models.DateTimeField(blank=True, null=True)),
                ('proteinortho_analysed', models.BooleanField(default=False)),
                ('proteinortho_error', models.IntegerField(blank=True, null=True)),
                ('registry', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, related_name='items', to='api.projectanalysisregistry')),
                ('task_proteinortho', models.OneToOneField(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='proteinortho_analysis_registry_item', to='django_celery_results.taskresult')),
            ],
        ),
        migrations.CreateModel(
            name='Profile',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('organization', models.CharField(blank=True, max_length=200, null=True)),
                ('is_admin', models.BooleanField(default=False)),
                ('account_confirmation', models.BooleanField(default=False)),
                ('account_confirmation_date', models.DateTimeField(blank=True, null=True)),
                ('first_login_date', models.DateTimeField(blank=True, null=True)),
                ('user', models.OneToOneField(on_delete=django.db.models.deletion.PROTECT, related_name='profile', to=settings.AUTH_USER_MODEL)),
            ],
        ),
    ]
