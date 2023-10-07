import logging
from subprocess import Popen, PIPE
from django.db import transaction
from api.models import *
import urllib.parse
import os.path
from root.engine.rsat_functions import *
from django.conf import settings
import root.lib.library as lib
from Bio import SeqIO

class RsatAnalyseEngine:
    def __init__(self, work_folder=None, pwms_folder=None, timeout=None):
        self.work_folder = work_folder
        self.pwms_folder = pwms_folder
        self.timeout = timeout

    def analyse_items(self, registry_id, organism_accession):
        tf_list = list()
        pwm_list = list()
        regulatory_interactions_list = list()
        prediction_list = list()
        target_pwm_file = settings.TARGET_PWMS_FILES_PATH+organism_accession+'/TF_Information.txt'

        item = ProjectAnalysisRegistry.objects.filter(pk=registry_id).first()

        if os.path.exists(target_pwm_file):
            parse_pwm_file(target_pwm_file, organism_accession)

            tf_list = select_tfs_by_organism(organism_accession)
            for tf in tf_list:
                pwm_list = select_pwms_by_locus_tag(tf.locus_tag)
                for pwm in pwm_list:
                    regulatory_interactions_list = select_regulatory_interactions_by_tf_locus_tag(pwm.locus_tag.locus_tag)
                    for regulatory_interaction in regulatory_interactions_list:
                        prediction = str(regulatory_interaction.tf_locus_tag.locus_tag+"-"+regulatory_interaction.tg_locus_tag.locus_tag+"-"+pwm.motif_id)
                        # this condition verify if exists a TFBS prediction already performed for this pwm
                        # to avoid duplicated predictions with the same tg_promoter and pwm
                        # note that pwms in Cis-Bp could be inferred by similarity tf_status = I
                        if(prediction not in prediction_list):
                            prediction_list.append(prediction)
                            promoter = select_promoter_by_locus_tag(regulatory_interaction.tg_locus_tag.locus_tag)
                            promoter_sequence = ">" + regulatory_interaction.tg_locus_tag.locus_tag + "\n" + promoter.promoter_seq
                            temp_file = self.work_folder+"/test.txt"

                            with open(temp_file, 'w') as f:
                                for line in promoter_sequence:
                                    f.write(line)
                            
                            # Prepare matrix
                            in_matrix = self.pwms_folder + pwm.motif_id + ".txt"
                            
                            # Input matrix (file content)
                            in_matrix_file = open(in_matrix, "r")
                            matrix = in_matrix_file.read()
                            in_matrix_file.close()
                            if(matrix == ("Pos	A	C	G	T"+'\n')):
                                lib.log.info("Null Matrix: " + pwm.motif_id)
                            else:
                                cmd = "matrix-scan -v 1 -matrix_format cis-bp -m "+in_matrix+" -pseudo 1 -decimals 1 -2str -origin start -bginput -markov 1 -bg_pseudo 0.01 -return limits -return sites -return pval -return rank -lth score 1 -uth pval 1e-4 -i "+temp_file+" -seq_format fasta -n score"
                
                                rsat_call = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
                                output, error = rsat_call.communicate()
                                ret = rsat_call.returncode
                                
                                if ret != 0:
                                    logging.info('Rsat analysis error for organism %s, Error: ' % organism_accession, output )
                                    self.__set_error(item, RsatErrorType.COMMAND_ERROR.value)
                                else:
                                    result = str(output)
                                    parts = result.strip().split("\\t")
                                    # Create new TFBS prediction for each RSAT prediction result
                                    strand = urllib.parse.unquote(parts[76])
                                    start = urllib.parse.unquote(parts[77])
                                    end = urllib.parse.unquote(parts[78])
                                    sequence = urllib.parse.unquote(parts[79])
                                    weight = urllib.parse.unquote(parts[80])
                                    pval = urllib.parse.unquote(parts[81])
                                    ln_pval = urllib.parse.unquote(parts[82])
                                    sig = urllib.parse.unquote(parts[83])

                                    tfbs = Tfbs(regulatory_interaction = RegulatoryInteraction.objects.get(pk=regulatory_interaction.id),
                                                organism_accession = Organism.objects.get(accession=organism_accession),
                                                tf_locus_tag = Gene.objects.get(locus_tag=regulatory_interaction.tf_locus_tag.locus_tag),
                                                tg_locus_tag = Gene.objects.get(locus_tag=regulatory_interaction.tg_locus_tag.locus_tag),
                                                pwm = Pwm.objects.get(pk=pwm.id), strand = strand, start = start, end = end, 
                                                sequence = sequence, weight = weight, pval = pval, ln_pval = ln_pval, sig = sig)
                                    if strand == 'ncol' and start == 'nrow' and end == 'pseudo' and sequence == 'Wmin' and weight == 'Wmax' and pval == 'Wrange\\n;':
                                        pass
                                    else:
                                        tfbs.save()
                        else:
                            lib.log.info("TFBS Prediction File already exists: " +regulatory_interaction.tf_locus_tag.locus_tag+"-"+regulatory_interaction.tg_locus_tag.locus_tag+"-"+pwm.motif_id+ ".txt")
            self.__set_analysed(item)
        else:
            lib.log.info('The tf information file for '+organism_accession+' doesnt exist')

    @staticmethod
    def __set_error(item, error_type):
        item.rsat_error = error_type
        item.rsat_analysed = False
        item.save(update_fields=['rsat_error', 'rsat_analysed'])

    @staticmethod
    def __set_analysed(item):
        item.rsat_analysed = True
        item.save(update_fields=['rsat_analysed'])