import logging

from celery import shared_task
from django.conf import settings
from django_celery_results.models import TaskResult

from funRegulationTool.task_utils import FunRegulationBaseTask
from root.engine.proteinortho_analyse_engine import ProteinOrthoAnalyseEngine
from root.engine.rsat_analyse_engine import RsatAnalyseEngine
from api.models import ProjectAnalysisRegistry

logger = logging.getLogger('main')

@shared_task(bind=True, name='run_proteinortho', base=FunRegulationBaseTask)
def task_run_proteinortho(self, registry_id, organism_accession):
    registry = ProjectAnalysisRegistry.objects.get(pk=registry_id)
    registry.task_proteinortho = TaskResult.objects.get(task_id=self.request.id)
    registry.save()
    logger.info(f'Running proteinOrtho for organism {organism_accession}')
    save_path_organism = settings.PROTEINORTHO_SAVE_PATH+organism_accession
    engine = ProteinOrthoAnalyseEngine(proteinOrtho_path=settings.PROTEINORTHO_PATH, work_folder=save_path_organism)
    engine.analyse_items(registry_id, organism_accession)
    logger.info(f'proteinOrtho analysis finished for organism {organism_accession}')
    
@shared_task(bind=True, name='run_rsat', base=FunRegulationBaseTask)
def task_run_rsat(self, registry_id, organism_accession):
    registry = ProjectAnalysisRegistry.objects.get(pk=registry_id)
    registry.task_rsat = TaskResult.objects.get(task_id=self.request.id)
    registry.save()
    logger.info(f'Running RSAT for organism {organism_accession}')
    engine = RsatAnalyseEngine(work_folder=settings.RSAT_SAVE_PATH+organism_accession, pwms_folder=settings.PWMS_FILES_PATH)
    engine.analyse_items(registry_id, organism_accession)
    logger.info(f'RSAT analysis finished for organism {organism_accession}')




    