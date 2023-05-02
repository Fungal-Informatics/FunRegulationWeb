from subprocess import Popen, PIPE
from django.db import transaction
from api.models import *
import urllib.parse
import os.path
import pandas as pd
from django.conf import settings
from root.engine.proteinOrtho_functions import select_protein_by_id, insert_orthology, construct_grn_orthology
from root.engine.proteinOrtho_functions import parse_protein_file, gff3_handler

class ProteinOrthoAnalyseEngine:
    def __init__(self, proteinOrtho_path=None, work_folder=None, timeout=None):
        self.proteinOrtho_path = proteinOrtho_path
        self.work_folder = work_folder
        self.timeout = timeout

    def analyse_items(self, organism_accession):
        # VERIFICAR SE JÁ EXISTE
        model_organism = self.__get_model_organism(organism_accession)
        target_organism_protein_file = settings.NCBI_DOWNLOAD_PATH+organism_accession+"/ncbi_dataset/data/"+organism_accession+"/protein.faa"
        target_organism_gff_file = settings.NCBI_DOWNLOAD_PATH+organism_accession+"/ncbi_dataset/data/"+organism_accession+"/genomic.gff"
        
        # item = ProjectAnalysisRegistryItem.objects.select_related('feature')\
        #     .filter(feature__organism__accession=organism_accession, active=True,
        #             feature__removed=False, feature__organism__removed=False)
                    
        command = ["perl", self.proteinOrtho_path, model_organism, target_organism_protein_file]

        if not os.path.exists(self.work_folder+"/proteinOrtho"):
            gff3_handler(target_organism_gff_file)
            parse_protein_file(target_organism_protein_file)
            os.makedirs(self.work_folder+"/proteinOrtho")

        os.chdir(self.work_folder+"/proteinOrtho")

        proc = Popen(command, stdout=PIPE, stderr=PIPE)
        output, error = proc.communicate()
        proc.communicate()
        filename = self.work_folder+"/proteinOrtho" + "/myproject.proteinortho.tsv"
        ret = proc.returncode
        if ret != 0:
            print(output)
            print('END OUTPUT')
            print(error)
            #self.__set_error(item, ProteinOrthoErrorType.COMMAND_ERROR.value)
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
                                print(target_protein)
                                # orthology = Orthology(model_protein,target_protein)
                                # if(orthology != None):
                                #     print(orthology)
                                #     print(orthology.target_protein)
                                    #insert_orthology(orthology)
            in_file.close()
            #construct_grn_orthology()

    @staticmethod
    def __get_model_organism(organism_accession):
        df = pd.read_csv(settings.ALL_ORGANISMS_FILE_PATH,sep = '\t')

        order_organism_user = df[['#order']][df['accession'] == organism_accession].squeeze()
        if(len(order_organism_user) > 0):
            model_organism = df[(df['#order'] == order_organism_user) & (df['model'] == True)]

            if(len(model_organism) > 0):
                order_model_organism = model_organism['#order'].squeeze()
                if(order_model_organism == 'Saccharomycetales'):
                    return settings.ORGANISM_MODEL_SACCHAROMYCES_CEREVISIAE_PROTEIN_PATH
                elif(order_model_organism == 'Eurotiales'):
                    return settings.ORGANISM_MODEL_A_NIDULANS_PROTEIN_PATH
                elif(order_model_organism == 'Sordariales'):
                    return settings.ORGANISM_MODEL_NEUROSPORA_CRASSA_PROTEIN_PATH
                elif(order_model_organism == 'Hypocreales'):
                    return settings.ORGANISM_MODEL_FUSARIUM_GRAMINEARUM_PROTEIN_PATH
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