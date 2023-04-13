from subprocess import Popen, PIPE
from django.db import transaction
from api.models import ProjectAnalysisRegistryItem, ProteinOrthoErrorType
import os.path
from general_functions import parse_orthology_file

class ProteinOrthoAnalyseEngine:
    def __init__(self, proteinOrtho_path=None, work_folder=None, timeout=None):
        self.proteinOrtho_path = proteinOrtho_path
        self.work_folder = work_folder
        self.timeout = timeout

    def analyse_items(self, id_items):
        item = ProjectAnalysisRegistryItem.objects.select_related('feature')\
            .filter(pk__in=id_items, active=True,
                    feature__removed=False, feature__gene__removed=False,
                    feature__gene__project__removed=False)
        organism1 = "/home/gabriel/Downloads/TCC I/Software/ProteinOrtho/proteinortho-master/test/C.faa"
        organism2 = "/home/gabriel/Downloads/TCC I/Software/ProteinOrtho/proteinortho-master/test/E.faa"
    
        command = ["perl", self.proteinOrtho_path, organism1, organism2]

        if not os.path.exists(self.work_folder):
            os.makedirs(self.work_folder)

        os.chdir(self.work_folder)

        proc = Popen(command, stdout=PIPE, stderr=PIPE)
        output, error = proc.communicate()
        filename = "/home/gabriel/Desktop/FunRegulationBack-end/FunRegulationAPI/funRegulationBackCode/projects/proteinOrtho_Results/myproject.proteinortho.tsv"
        ret = proc.returncode
        if ret != 0:
            self.__set_error(item, ProteinOrthoErrorType.COMMAND_ERROR.value)
            
        else:
            parse_orthology_file(filename)




    @staticmethod
    def __set_error(item, error_type):
        item.proteinOrtho_error = error_type
        item.save(update_fields=['proteinOrtho_error'])

    @staticmethod
    def __set_analysed(item):
        item.proteinOrtho_analysed = True
        item.save(update_fields=['proteinOrtho_analysed'])