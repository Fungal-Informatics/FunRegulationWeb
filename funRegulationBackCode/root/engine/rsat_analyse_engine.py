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

    def analyse_items(self):
        tf_list = list()
        pwm_list = list()
        regulatory_interactions_list = list()
        out_folder = self.work_folder
        organism_accession = "GCA_000240135.3"

        genome = SeqIO.to_dict(SeqIO.parse(open(settings.ORGANISM_MODEL_FUSARIUM_GRAMINEARUM_PATH), 'fasta'))

        #rsat_pwms = "/home/gabriel/packages/rsat/perl-scripts/M01890_2.00.txt"
        #rsat_promoter = "/home/gabriel/packages/rsat/perl-scripts/FOXG_10713_promoter.fasta"

        tf_list = select_tfs_by_organism(organism_accession)
        for tf in tf_list:
            pwm_list = select_pwms_by_locus_tag(tf.locus_tag)
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




    @staticmethod
    def __set_error(item, error_type):
        item.rsat_error = error_type
        item.save(update_fields=['rsat_error'])

    @staticmethod
    def __set_analysed(item):
        item.rsat_analysed = True
        item.save(update_fields=['rsat_analysed'])