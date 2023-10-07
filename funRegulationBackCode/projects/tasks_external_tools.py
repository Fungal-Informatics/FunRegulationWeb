import logging

from celery import shared_task, chain
from django.conf import settings
from django.db import transaction
from django_celery_results.models import TaskResult

from funRegulationTool.task_utils import FunRegulationBaseTask
from root.engine.proteinortho_analyse_engine import ProteinOrthoAnalyseEngine
from root.engine.rsat_analyse_engine import RsatAnalyseEngine
from api.models import ProjectAnalysisRegistry

@shared_task(bind=True, name='run_proteinortho', base=FunRegulationBaseTask)
def task_run_proteinortho(self, registry_id, organism_accession):
    registry = ProjectAnalysisRegistry.objects.get(pk=registry_id)
    registry.task_proteinortho = TaskResult.objects.get(task_id=self.request.id)
    registry.save()
    logging.info('Running proteinOrtho for organism %s' % organism_accession)
    save_path_organism = settings.PROTEINORTHO_SAVE_PATH+organism_accession
    engine = ProteinOrthoAnalyseEngine(proteinOrtho_path=settings.PROTEINORTHO_PATH, work_folder=save_path_organism)
    engine.analyse_items(registry_id, organism_accession)
    logging.info('proteinOrtho analysis finished for organism %s' % organism_accession)
    
@shared_task(bind=True, name='run_rsat', base=FunRegulationBaseTask)
def task_run_rsat(self, registry_id, organism_accession):
    registry = ProjectAnalysisRegistry.objects.get(pk=registry_id)
    registry.task_rsat = TaskResult.objects.get(task_id=self.request.id)
    registry.save()
    logging.info('Running RSAT for organism %s' % organism_accession)
    engine = RsatAnalyseEngine(work_folder=settings.RSAT_SAVE_PATH+organism_accession, pwms_folder=settings.PWMS_FILES_PATH)
    engine.analyse_items(registry_id, organism_accession)
    logging.info('RSAT analysis finished for organism %s' % organism_accession)




    