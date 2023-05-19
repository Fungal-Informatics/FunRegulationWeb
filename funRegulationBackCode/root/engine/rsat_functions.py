import psycopg2
import root.lib.library as lib
import urllib.parse
from api.models import *
import os
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from django.conf import settings

#create log file
log_name = os.path.join(settings.LOG_FILE_PATH)
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)

def parse_pwm_file(in_file_pwm):
    lib.log.info("Parsing "+ in_file_pwm)
    with open(in_file_pwm) as in_file:
        for line in in_file:
            if line.startswith("#"): 
                continue
            line_parts = line.strip().split("\t")
            
            motif_id = urllib.parse.unquote(line_parts[3]) 
            if motif_id != '.' and motif_id != 'Motif_ID':
                locus_tag = urllib.parse.unquote(line_parts[5])
                status = urllib.parse.unquote(line_parts[8])
                tf_family = urllib.parse.unquote(line_parts[9])
                motif_type = urllib.parse.unquote(line_parts[14])
                msource = urllib.parse.unquote(line_parts[16])
                msource_author = urllib.parse.unquote(line_parts[17])
                pubmedid = urllib.parse.unquote(line_parts[19])
                if pubmedid == 'NULL':
                    pubmedid = ''
                pwm = Pwm(locus_tag=Gene.objects.get(locus_tag=locus_tag), motif_id=motif_id, status=status, tf_family=tf_family, motif_type=motif_type, msource_author=msource_author, msource=msource, pubmedid=pubmedid)
                #insert_pwm(pwm)
    in_file.close()
    lib.log.info(in_file_pwm + " parsed correctly")

def select_tfs_by_organism(organism_accession):
    genes = []
    genes = Gene.objects.filter(is_tf=True).filter(organism_accession=organism_accession)
    return genes
    # dbConnection = create_db_connection()
    # tf_list = list()
    # try:
    #     cursor = dbConnection.cursor()
    #     cursor.execute("SELECT * from gene WHERE is_tf = 'True' and organism_accession = %s", (organism_accession,))
    #     rec = cursor.fetchall()
    #     #queryset = Gene.objects.filter(is_tf = 'true').filter(organism_accession = organism_accession)
    #     for row in rec:
    #         tf_list.append(Gene(row[0],row[1],row[2],row[3],row[4]))
    #     cursor.close()
    #     return tf_list
    # except (Exception, psycopg2.Error) as error:
    #     lib.log.info("Failed to execute the select into table gene", error)
    #     lib.log.info(organism_accession)

def select_pwms_by_locus_tag(locus_tag):
    pwms = list()
    try:
        pwms = Pwm.objects.filter(locus_tag=locus_tag)
        return pwms
    except(Exception) as error:
        lib.log.info("Failed to execute the select into table pwm", error)
        lib.log.info(locus_tag)

def select_regulatory_interactions_by_tf_locus_tag(tf_locus_tag):
    dbConnection = create_db_connection()
    regulatory_interactions = list()
    try:
        cursor = dbConnection.cursor()
        #postgreSQL_select_Query = "SELECT * from regulatory_interaction WHERE tf_locus_tag = %s"
        cursor.execute("SELECT * from regulatory_interaction WHERE tf_locus_tag = %s", (tf_locus_tag,))
        rec = cursor.fetchall()
        for row in rec:
            regulatory_interactions.append(RegulatoryInteraction(row[0],row[1],row[2],row[3],row[4]))
        cursor.close()
        return regulatory_interactions
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table regulatory_interaction", error)
        lib.log.info(tf_locus_tag)

def extract_promoter(genome, tg_locus_tag):
    promoter = select_promoter_by_locus_tag(tg_locus_tag)
    
    long_seq_record = genome[promoter.source]
    long_seq = long_seq_record.seq
    short_seq = str(long_seq)[promoter.start:promoter.stop]
    if promoter.strand == '-':
        short_seq = str(long_seq)[promoter.stop:promoter.start]
        my_dna = Seq(short_seq)
        my_dna = my_dna.reverse_complement()
        short_seq=str(my_dna)
        if promoter.stop > len(my_dna)+promoter.start :
            lib.log.info("Promoter of gene " + promoter.locus_tag + " can't be fully extracted")
    
    short_seq_record = SeqRecord(Seq(short_seq), id=promoter.locus_tag, name=promoter.locus_tag, description=promoter.source+'_'+str(promoter.start)+':'+str(promoter.stop))
    return short_seq_record

def select_promoter_by_locus_tag(tg_locus_tag):
    dbConnection = create_db_connection()
    try:
        cursor = dbConnection.cursor()
        postgreSQL_select_Query = "SELECT * from promoter WHERE locus_tag = %s"
        cursor.execute(postgreSQL_select_Query, (tg_locus_tag,))
        rec = cursor.fetchone()
        promoter = Promoter(rec[0],rec[1],rec[2],rec[3],rec[4])
        cursor.close()
        return promoter
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table promoter", error)
        lib.log.info(tg_locus_tag)

def insert_tfbs_prediction(tfbs):
    dbConnection = create_db_connection()
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO tfbs VALUES (default, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
                                                            (tfbs.regulatory_interaction_id, 
                                                            tfbs.pwm_id, 
                                                            tfbs.strand, 
                                                            tfbs.start, 
                                                            tfbs.end, 
                                                            tfbs.sequence, 
                                                            tfbs.weight, 
                                                            tfbs.pval, 
                                                            tfbs.ln_pval, 
                                                            tfbs.sig))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE tfbs")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE tfbs ")
        lib.log.info(str(tfbs.regulatory_interaction_id) + " " +
                    str(tfbs.pwm_id) + " " + 
                    str(tfbs.strand) + " " +
                    str(tfbs.sequence))
        
def create_db_connection():
    try:
        con = psycopg2.connect(host='localhost', database='funregulationtcc',
        user='postgres', password='postgres')
        lib.log.info("Successfully Connected to PostgreSQL")
        return con
    except (Exception, psycopg2.Error) as error:
        lib.log.info(error)