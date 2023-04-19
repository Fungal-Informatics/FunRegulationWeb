from subprocess import Popen, PIPE
from django.db import transaction
from api.models import *
import urllib.parse
import os.path
from root.engine.proteinOrtho_functions import select_protein_by_id, insert_orthology, construct_grn_orthology

class ProteinOrthoAnalyseEngine:
    def __init__(self, proteinOrtho_path=None, work_folder=None, timeout=None):
        self.proteinOrtho_path = proteinOrtho_path
        self.work_folder = work_folder
        self.timeout = timeout

    def analyse_items(self):
        model_organism = "/home/gabriel/Downloads/TCC I/Software/ProteinOrtho/proteinortho-master/test/C.faa"
        target_organism = "/home/gabriel/Downloads/TCC I/Software/ProteinOrtho/proteinortho-master/test/E.faa"
    
        command = ["perl", self.proteinOrtho_path, model_organism, target_organism]

        if not os.path.exists(self.work_folder):
            os.makedirs(self.work_folder)

        os.chdir(self.work_folder)

        proc = Popen(command, stdout=PIPE, stderr=PIPE)
        output, error = proc.communicate()
        filename = self.work_folder + "/myproject.proteinortho.tsv"
        ret = proc.returncode
        if ret != 0:
            pass
            #self.__set_error("ERROR EXECUTING PROTEINORTHO", ProteinOrthoErrorType.COMMAND_ERROR.value)
        else:
            with open(filename) as in_file:
                for line in in_file:
                    if line.startswith("#"): continue
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
                                orthology = Orthology(model_protein,target_protein)
                                if(orthology != None):
                                    insert_orthology(orthology)
            in_file.close()
            construct_grn_orthology()




    @staticmethod
    def __set_error(item, error_type):
        item.proteinortho_error = error_type
        item.save(update_fields=['proteinortho_error'])

    @staticmethod
    def __set_analysed(item):
        item.proteinortho_analysed = True
        item.save(update_fields=['proteinortho_analysed'])