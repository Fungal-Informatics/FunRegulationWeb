import logging
from subprocess import Popen, PIPE
from django.db import transaction
from api.models import *
import urllib.parse
import os.path
from django.conf import settings
from root.engine.proteinOrtho_functions import select_protein_by_id, gbff_handler, construct_grn_orthology

class ProteinOrthoAnalyseEngine:
    def __init__(self, proteinOrtho_path=None, work_folder=None, timeout=None):
        self.proteinOrtho_path = proteinOrtho_path
        self.work_folder = work_folder
        self.timeout = timeout

    def analyse_items(self, registry_id, organism_accession):
        model_organism = []
        model_organism = self.__get_model_organism(organism_accession)
        target_organism_protein_file = settings.NCBI_DOWNLOAD_PATH+organism_accession+"/ncbi_dataset/data/"+organism_accession+"/protein.faa"
        target_organism_gbff_file = settings.NCBI_DOWNLOAD_PATH+organism_accession+"/ncbi_dataset/data/"+organism_accession+"/genomic.gbff"
                    
        command = ["perl", self.proteinOrtho_path, model_organism[0], target_organism_protein_file]

        item = ProjectAnalysisRegistry.objects.filter(pk=registry_id).first()
        organism_analysed = ProjectAnalysisRegistry.objects.filter(organism_accession=organism_accession)\
                .filter(proteinortho_analysed=True).count()
        
        if organism_analysed <= 0:
            if not os.path.exists(self.work_folder+"/proteinOrtho"):
                gbff_handler(organism_accession, target_organism_gbff_file)
                os.makedirs(self.work_folder+"/proteinOrtho")

            os.chdir(self.work_folder+"/proteinOrtho")

            proc = Popen(command, stdout=PIPE, stderr=PIPE)
            output, error = proc.communicate()
            
            filename = self.work_folder+"/proteinOrtho" + "/myproject.proteinortho.tsv"
            ret = proc.returncode

            if ret != 0:
                logging.info('proteinOrtho analysis error for organism %s, Error: ' % organism_accession, output )
                self.__set_error(item, ProteinOrthoErrorType.COMMAND_ERROR.value)
            else:
                with open(filename) as in_file:
                    for line in in_file:
                        if line.startswith("#"): 
                            continue
                        line_parts = line.strip().split("\t")
                        
                        model = urllib.parse.unquote(line_parts[3])
                        target = urllib.parse.unquote(line_parts[4])
                        
                        model_parts = model.strip().split(",")
                        target_parts = target.strip().split(",")

                        for record_model in model_parts:
                            for record_target in target_parts:
                                if (record_model != '*' and record_target != '*'):
                                    model_protein = select_protein_by_id(record_model)
                                    target_protein = select_protein_by_id(record_target)
                                    orthology = Orthology(model_protein=Protein.objects.get(locus_tag=model_protein),target_protein=Protein.objects.get(locus_tag=target_protein))

                                    if(orthology != None):
                                        orthology.save()

                in_file.close()
                construct_grn_orthology(model_organism[1], organism_accession)
                self.__set_analysed(item)
        else:
            logging.info('The organism {} already has the proteinOrtho data in database'.format(organism_accession))

    @staticmethod
    def __get_model_organism(organism_accession):
        results = []
        order_organism_user = Organism.objects.filter(accession=organism_accession).values('order')

        if(len(order_organism_user) > 0):
            for order_value in order_organism_user:
                order = order_value['order']

            model_organism = Organism.objects.filter(order=order).filter(is_model=True).values('order')

            if(len(model_organism) > 0):
                for order_model_value in model_organism:
                    order_model = order_model_value['order']
                
                if(order_model == 'Saccharomycetales'):
                    results.append(settings.ORGANISM_MODEL_SACCHAROMYCES_CEREVISIAE_PROTEIN_PATH)
                    results.append('R64-3-1-SGD') 
                elif(order_model == 'Eurotiales'):
                    results.append(settings.ORGANISM_MODEL_A_NIDULANS_PROTEIN_PATH)
                    results.append('s10-m04-r16-AspGD') 
                elif(order_model == 'Sordariales'):
                    results.append(settings.ORGANISM_MODEL_NEUROSPORA_CRASSA_PROTEIN_PATH)
                    results.append('GCA_000182925.2') 
                elif(order_model == 'Hypocreales'):
                    results.append(settings.ORGANISM_MODEL_FUSARIUM_GRAMINEARUM_PROTEIN_PATH)
                    results.append('GCA_000240135.3')
                return results
            else:
                return 'ERROR - NOT FOUND A MODEL ORGANISM FOR UPLOADED ORGANISM'
        else:
            return 'ERROR - NOT FOUND ORGANISM WITH THIS ACCESSION'

    @staticmethod
    def __set_error(item, error_type):
        item.proteinortho_error = error_type
        item.save(update_fields=['proteinortho_error'])

    @staticmethod
    def __set_analysed(item):
        item.proteinortho_analysed = True
        item.save(update_fields=['proteinortho_analysed'])