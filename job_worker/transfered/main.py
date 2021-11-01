"""
###########################################################
#
# load_data.py
#
# FunRegulation: Network Inference (Alexandre Rafael Lenz)
# Universidade do Estado da Bahia (UNEB)
# https://github.com/Fungal-Informatics/funregulation
# Last update: 21.09.2021 (Lenz, A. R.)
#
#   INPUT DATA:
#
# a) Genomic data:
#    Data were collected and organized in tab-delimited files.
#    i) Organisms
#    ii) Genes
#    iii) Proteins
#    iv) Annotations
#
# b) Cis-BP:
#    Data were collected and organized in tab-delimited files for each species.
#    i) TF and PWM data
#
###########################################################
"""
GRNInferenceVersion = '1.0 - 04/10/2021'

import psycopg2
import sys, re
import os, platform
import urllib.parse
import lib.library as lib
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import namedtuple
from suds.client import Client

#Initialize Folder Paths
in_folder = "/Users/arlenz/Documents/Alexandre/Trabalho/Pesquisa/FunWorld/Tools/FunRegulation/Software/Database/"
out_folder = "/Users/arlenz/Documents/Alexandre/Trabalho/Pesquisa/FunWorld/Tools/FunRegulation/Software/Database/Results/"

#Initialize file Paths
in_file_organisms = os.path.join(in_folder,'Ensembl_Species.tsv')
in_file_genes = os.path.join(in_folder,'Model/Aspergillus_nidulans/A_nidulans_FGSC_A4_current_features.gff3')
in_file_proteins = os.path.join(in_folder,'Model/Aspergillus_nidulans/A_nidulans_FGSC_A4_current_orf_trans_all.fasta')
in_file_pwms = os.path.join(in_folder,'Cis-BP/Species/Model/Aspergillus_nidulans_2021_09_17_6_25_pm/TF_Information.txt')
in_file_model_regulatory = os.path.join(in_folder,'Model/Aspergillus_nidulans.tsv')
in_file_orthology = os.path.join(in_folder,'Model/Fusarium_graminearum/Anidulans.proteinortho.tsv')
in_file_genome = os.path.join(in_folder,'Model/Fusarium_graminearum/Fusarium_graminearum_gca_000240135.ASM24013v3.dna.toplevel.fa')
#Initialize Output File Paths
#...

#DB definitions
dbConnection = None;

"""
################ PROMOTER LENGHT DEFINITION ##############
"""
upstream = -1000
downstream = 0

"""
##############   NAMEDTUPLE  DEFINITIONS  #################
"""
# Initialized GeneInfo named tuple to handle with GFF3 annotation. Note: namedtuple is immutable
gffInfoFields = ["seqid", "source", "ltype", "start", "end", "score", "strand", "phase", "attributes"]
GFFRecord = namedtuple("GFFRecord", gffInfoFields)

#create log file
log_name = os.path.join(out_folder, 'logfiles', 'funregulation.log')
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)






                     

                    








"""
    Main function of this program
""" 
if __name__ == '__main__':
    
    #Create a database connection
    dbConnection = create_db_connection()
    
    #Load Input Files
    #parse_organism_file(in_file_organisms)
    #gff3_handler(in_file_genes)
    #parse_protein_file(in_file_proteins)
    #parse_pwm_file(in_file_pwms)
    #parse_model_regulatory_file(in_file_model_regulatory)
    
    #Construct GRN
    #####Orthology
    #execute ProteinOrtho *****
    #parse_orthology_file(in_file_orthology)
    #construct_grn_orthology('2', '4')
    #####TFBS Predictions
    construct_grn_tfbs_predictions('4')
    
    if (dbConnection):
        dbConnection.close()
        lib.log.info("The DB connection is closed")
