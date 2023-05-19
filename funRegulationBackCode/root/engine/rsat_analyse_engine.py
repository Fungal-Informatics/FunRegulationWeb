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
    def __init__(self, work_folder=None, pwms_folder = None, timeout=None):
        self.work_folder = work_folder
        self.pwms_folder = pwms_folder
        self.timeout = timeout

    def analyse_items(self, organism_accession):
        tf_list = list()
        pwm_list = list()
        regulatory_interactions_list = list()
        target_pwm_file = settings.TARGET_PWMS_FILES_PATH+organism_accession+'/TF_Information.txt'
        out_folder = self.work_folder

        if os.path.exists(target_pwm_file):
            parse_pwm_file(target_pwm_file)
            dna_model_file = self.__get_dna_model_organism(organism_accession)

            if(dna_model_file.find('ERROR') == -1):
                genome = SeqIO.to_dict(SeqIO.parse(open(dna_model_file), 'fasta'))

                tf_list = select_tfs_by_organism(organism_accession)
                for tf in tf_list:
                    pwm_list = select_pwms_by_locus_tag(tf.locus_tag)
                    print(pwm_list)
                    for pwm in pwm_list:
                        regulatory_interactions_list = select_regulatory_interactions_by_tf_locus_tag(pwm.locus_tag.locus_tag)
                        for regulatory_interaction in regulatory_interactions_list:
                            print('TEST REGULATORY_INTERACTION')
                            prediction_file = out_folder + "tfbs_predictions/" + organism_accession + "/" +regulatory_interaction.tf_locus_tag+"-"+regulatory_interaction.tg_locus_tag+"-"+pwm.motif_id+ ".txt"
                            # this condition verify if exists a TFBS prediction already performed for this pwm
                            # to avoid duplicated predictions with the same tg_promoter and pwm
                            # note that pwms in Cis-Bp could be inferred by similarity tf_status = I
                            if not os.path.isfile(prediction_file):
                                promoter_sequence = extract_promoter(genome, regulatory_interaction.tg_locus_tag)
                                
                                ## Prepare matrix
                                in_matrix = self.pwms_folder + pwm.motif_id + ".txt"

                                # Input matrix (file content)
                                in_matrix_file = open(in_matrix, "r")
                                matrix = in_matrix_file.read()
                                in_matrix_file.close()

                                cmd = "matrix-scan -v 1 -matrix_format cis-bp -m "+in_matrix+" -pseudo 1 -decimals 1 -2str -origin end -bginput -markov 1 -bg_pseudo 0.01 -return limits -return sites -return pval -return rank -lth score 1 -uth pval 1e-4 -i "+promoter_sequence+" -seq_format fasta -n score"
                
                                rsat_call = Popen(cmd, shell=True,stdout=PIPE,stderr=PIPE)
                                out, error = rsat_call.communicate()
                                print("ENTERED IN RSAT CALL SOFTWARE!!")
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
            else:
                lib.log.info(dna_model_file)
        else:
            lib.log.info('The tf information file for '+organism_accession+' doesnt exist')

    @staticmethod
    def __get_dna_model_organism(organism_accession):
        order_organism_user = Organism.objects.filter(accession=organism_accession).values('order')

        if(len(order_organism_user) > 0):
            for order_value in order_organism_user:
                order = order_value['order']

            model_organism = Organism.objects.filter(order=order).filter(is_model=True).values('order')

            if(len(model_organism) > 0):
                for order_model_value in model_organism:
                    order_model = order_model_value['order']
                
                if(order_model == 'Saccharomycetales'):
                    return settings.ORGANISM_MODEL_SACCHAROMYCES_CEREVISIAE_PROTEIN_PATH
                elif(order_model == 'Eurotiales'):
                    return settings.ORGANISM_MODEL_SACCHAROMYCES_CEREVISIAE_PATH
                elif(order_model == 'Sordariales'):
                    return settings.ORGANISM_MODEL_NEUROSPORA_CRASSA_PATH 
                elif(order_model == 'Hypocreales'):
                    return settings.ORGANISM_MODEL_FUSARIUM_GRAMINEARUM_PATH
            else:
                return 'ERROR - NOT FOUND A MODEL ORGANISM FOR UPLOADED ORGANISM'
        else:
            return 'ERROR - NOT FOUND ORGANISM WITH THIS ACCESSION'

    @staticmethod
    def __set_error(item, error_type):
        item.rsat_error = error_type
        item.save(update_fields=['rsat_error'])

    @staticmethod
    def __set_analysed(item):
        item.rsat_analysed = True
        item.save(update_fields=['rsat_analysed'])