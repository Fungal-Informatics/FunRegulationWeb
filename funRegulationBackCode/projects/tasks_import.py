import logging
import root.lib.library as lib
from django.conf import settings
import os.path
from zipfile import ZipFile
from celery import shared_task
from django_celery_results.models import TaskResult
from funRegulationTool.task_utils import FunRegulationBaseTask
from api.models import ProjectAnalysisRegistry

import sys
from typing import List

from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
from ncbi.datasets.openapi import ApiException as DatasetsApiException
from ncbi.datasets import GenomeApi as DatasetsGenomeApi
from ncbi.datasets.metadata.genome import get_assembly_metadata_by_bioproject_accessions

from ncbi.datasets.package import dataset

@shared_task(bind=True, name='import_genes', base=FunRegulationBaseTask)
def task_import_genes(self, registry_id, organism_accession):
    registry = ProjectAnalysisRegistry.objects.get(pk=registry_id)
    registry.task_download_organism = TaskResult.objects.get(task_id=self.request.id)
    registry.save()
    download_path = settings.NCBI_DOWNLOAD_PATH
    accessions: List[str] = [organism_accession]
    zipfile_name = organism_accession+".zip"
    
    if not os.path.exists(download_path+organism_accession):
        os.makedirs(download_path+organism_accession)
    
    save_path = os.path.join(download_path+organism_accession, zipfile_name)

    with DatasetsApiClient() as api_client:
        genome_api = DatasetsGenomeApi(api_client)
        try:
            print("Begin download of genome data package ...")
            lib.log.info(f"Begin download of genome {organism_accession} data package ...")
            genome_ds_download = genome_api.download_assembly_package(
                accessions,
                include_annotation_type=["PROT_FASTA", "GENOME_GFF", "GENOME_GBFF"],
                _preload_content=False,
            )

            with open(save_path, "wb") as f:
                f.write(genome_ds_download.data)
            print(f"Download completed -- see {zipfile_name}")
            lib.log.info(f"Download completed of genome {organism_accession}")
            work_dir = download_path+organism_accession
            os.chdir(work_dir)
            with ZipFile(save_path, 'r') as zip:
                zip.extract(zip.namelist()[2]) #genome
                zip.extract(zip.namelist()[3]) #GFF3
                zip.extract(zip.namelist()[4]) #protein
                zip.extract(zip.namelist()[5]) #sequence
        except DatasetsApiException as e:
            registry.download_completed = False
            registry.save()
            sys.exit(f"Exception when calling download_assembly_package: {e}\n")
    registry.download_completed = True
    registry.save()