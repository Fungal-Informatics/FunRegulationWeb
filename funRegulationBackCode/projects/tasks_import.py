import logging
from django.conf import settings
import os.path
from zipfile import ZipFile
from celery import shared_task
from django_celery_results.models import TaskResult
from funRegulationTool.task_utils import FunRegulationBaseTask   

import sys
from typing import List

from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
from ncbi.datasets.openapi import ApiException as DatasetsApiException
from ncbi.datasets import GenomeApi as DatasetsGenomeApi
from ncbi.datasets.metadata.genome import get_assembly_metadata_by_bioproject_accessions

from ncbi.datasets.package import dataset

def import_genes():
    #result = task_import_genes.apply_async((import_registry.pk,))
    result = task_import_genes.delay()
    #import_registry.task = TaskResult.objects.get(task_id=result.task_id)
    #import_registry.save()

#@shared_task(bind=True, name='import_genes', base=FunRegulationBaseTask)
@shared_task(bind=True, name='import_genes')
def task_import_genes(self):
    organism_accession = "GCA_003184765.1"
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
            genome_ds_download = genome_api.download_assembly_package(
                accessions,
                include_annotation_type=["PROT_FASTA"],
                _preload_content=False,
            )

            with open(save_path, "wb") as f:
                f.write(genome_ds_download.data)
            print(f"Download completed -- see {zipfile_name}")
            work_dir = download_path+organism_accession
            os.chdir(work_dir)
            with ZipFile(save_path, 'r') as zip:
                zip.extract(zip.namelist()[2]) #genome
                zip.extract(zip.namelist()[3]) #protein
                zip.extract(zip.namelist()[4]) #sequence
        except DatasetsApiException as e:
            sys.exit(f"Exception when calling download_assembly_package: {e}\n")
