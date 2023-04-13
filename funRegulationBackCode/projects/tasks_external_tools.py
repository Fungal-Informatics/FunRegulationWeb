import logging

from celery import shared_task
from django.conf import settings
from django.db import transaction
from django_celery_results.models import TaskResult

from funRegulationTool.task_utils import FunRegulationBaseTask
from root.engine.proteinortho_analyse_engine import ProteinOrthoAnalyseEngine
from asyncio.subprocess import STDOUT
from subprocess import Popen, PIPE
import urllib.parse
import sys
import os.path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from api.models import *
import root.lib.library as lib
import psycopg2
import root.general_functions as LenzFunctions
#from root.general_functions import select_tfs_by_organism

def run_proteinortho(import_registry_id):
    result = task_run_proteinortho.delay()
    # result = task_run_proteinortho.apply_async((import_registry.pk,))
    # import_registry.task = TaskResult.objects.get(task_id=result.task_id)
    # import_registry.save()
    
def run_rsat():
    result = task_run_rsat.delay()

@shared_task(bind=True, name='run_rsat', base=FunRegulationBaseTask)
def task_run_rsat(self):
    tf_list = list()
    pwm_list = list()
    regulatory_interactions_list = list()
    out_folder = settings.RSAT_SAVE_PATH
    organism_accession = "GCA_000240135.3"

    genome = SeqIO.to_dict(SeqIO.parse(open(settings.ORGANISM_MODEL_FUSARIUM_GRAMINEARUM_PATH), 'fasta'))

    #rsat_pwms = "/home/gabriel/packages/rsat/perl-scripts/M01890_2.00.txt"
    #rsat_promoter = "/home/gabriel/packages/rsat/perl-scripts/FOXG_10713_promoter.fasta"
    tf_list = LenzFunctions.select_tfs_by_organism(organism_accession)
    for tf in tf_list:
        pwm_list = LenzFunctions.select_pwms_by_locus_tag(tf.locus_tag)
        for pwm in pwm_list:
            #print(pwm.locus_tag.locus_tag)
            regulatory_interactions_list = LenzFunctions.select_regulatory_interactions_by_tf_locus_tag(pwm.locus_tag.locus_tag)
            for regulatory_interaction in regulatory_interactions_list:
                prediction_file = out_folder + "tfbs_predictions/" + organism_accession + "/" +regulatory_interaction.tf_locus_tag+"-"+regulatory_interaction.tg_locus_tag+"-"+pwm.motif_id+ ".txt"
                # this condition verify if exists a TFBS prediction already performed for this pwm
                # to avoid duplicated predictions with the same tg_promoter and pwm
                # note that pwms in Cis-Bp could be inferred by similarity tf_status = I
                if not os.path.isfile(prediction_file):
                    #Extract Promoter
                    promoter_sequence = LenzFunctions.extract_promoter(genome, regulatory_interaction.tg_locus_tag)
                    
                    ## Prepare matrix
                    in_matrix = settings.PWMS_FILES_PATH + pwm.motif_id + ".txt"

                    # Input matrix (file content)
                    in_matrix_file = open(in_matrix, "r")
                    matrix = in_matrix_file.read()
                    in_matrix_file.close()

                    cmd = "matrix-scan -v 1 -matrix_format cis-bp -m "+in_matrix+" -pseudo 1 -decimals 1 -2str -origin end -bginput -markov 1 -bg_pseudo 0.01 -return limits -return sites -return pval -return rank -lth score 1 -uth pval 1e-4 -i "+promoter_sequence+" -seq_format fasta -n score"
    
                    rsat_call = Popen(cmd, shell=True,stdout=PIPE,stderr=PIPE)
                    out, error = rsat_call.communicate()

                    #out_file = out_folder + "temp.txt"
                    # lib.log.info("Call RSAT: "+'\n'+
                    #          " tf_locus_tag: "+ regulatory_interaction.tf_locus_tag+'\n'+
                    #          " tg_locus_tag: "+ regulatory_interaction.tg_locus_tag+'\n'+
                    #          " pwm_id: "+ pwm.motif_id)
                    
                    ret = rsat_call.returncode
                    if ret != 0:
                        print("RSAT failed %d %s %s" % (rsat_call.returncode, out, error))
                    else:
                        if matrix == ("Pos	A	C	G	T"+'\n'):
                            lib.log.info("Null Matrix: " + pwm.motif_id)
                        else:
                            result = out
                            # Write result in output file
                            with open(prediction_file, 'w') as out_file:
                                out_file.write(result)
                            out_file.close()

                            with open (prediction_file) as in_file:
                                for line in in_file:
                                    if line.startswith("#"): 
                                        continue
                                    parts = line.strip().split("\\t")
                                    #Create new TFBS prediction for each RSAT prediction result
                                    strand = urllib.parse.unquote(parts[76])
                                    start = urllib.parse.unquote(parts[77])
                                    end = urllib.parse.unquote(parts[78])
                                    sequence = urllib.parse.unquote(parts[79])
                                    weight = urllib.parse.unquote(parts[80])
                                    pval = urllib.parse.unquote(parts[81])
                                    ln_pval = urllib.parse.unquote(parts[82])
                                    sig = urllib.parse.unquote(parts[83])
                                    tfbs = Tfbs(0, regulatory_interaction.id, pwm.id, strand, start, end, sequence, weight, pval, ln_pval, sig)
                                    #insert_tfbs_prediction(tfbs)
                                in_file.close()
                else:
                    lib.log.info("TFBS Prediction File already exists: " +regulatory_interaction.tf_locus_tag+"-"+regulatory_interaction.tg_locus_tag+"-"+pwm.motif_id+ ".txt")

@shared_task(bind=True, name='run_proteinortho', base=FunRegulationBaseTask)
def task_run_proteinortho(self, items):
    logging.info('importing genes based on registry %s' % items)
    organism_accession = "GCA_003184765.3"
    save_path_organism = settings.PROTEINORTHO_SAVE_PATH+organism_accession
    engine = ProteinOrthoAnalyseEngine(proteinOrtho_path=settings.PROTEINORTHO_PATH, work_folder=save_path_organism)
    engine.analyse_items(items)
    logging.info('proteinOrtho analysis finished for items %s' % items)
    