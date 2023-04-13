import psycopg2
import os, platform
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import urllib.parse
from collections import namedtuple
import lib.library as lib
#from modelsLenz.models import *
from api.models import *

#Initialize Folder Paths
in_folder = ""
out_folder = ""

#Initialize file Paths
in_file_orthology = os.path.join('Model/Fusarium_graminearum/Anidulans.proteinortho.tsv')
in_file_genome = os.path.join('Model/Fusarium_graminearum/Fusarium_graminearum_gca_000240135.ASM24013v3.dna.toplevel.fa')

upstream = -1000
downstream = 0

dbConnection = None

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
        
def update_gene(gene):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("UPDATE gene SET symbol_gene = %s, description = %s, is_tf = %s WHERE organism = %s AND locus_tag = %s",
                                (gene.symbol_gene, 
                                gene.description,
                                gene.is_tf,
                                gene.organism,
                                gene.locus_tag))
        dbConnection.commit()
        lib.log.info("Record updated successfully into TABLE gene")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to update data into TABLE gene", error)
        lib.log.info(str(gene.organism) + " " +
                     str(gene.locus_tag) + " " + 
                     str(gene.symbol_gene) + " " +
                     str(gene.description) + " " +
                     str(gene.is_tf))
        
def select_gene_by_locus_tag(locus_tag):
    gene = None
    try:
        cursor = dbConnection.cursor()
        postgreSQL_select_Query = "SELECT * FROM gene WHERE locus_tag = %s"
        cursor.execute(postgreSQL_select_Query, (locus_tag,))
        rec = cursor.fetchall()
        for row in rec:
            gene = Gene(row[0],row[1],row[2],row[3],row[4])
            return gene
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table gene", error)

def select_protein_by_id(protein_id):
    dbConnection = create_db_connection()
    protein = None
    try:
        cursor = dbConnection.cursor()
        postgreSQL_select_Query = "SELECT * FROM protein WHERE id = %s"
        cursor.execute(postgreSQL_select_Query, (protein_id,))
        rec = cursor.fetchall()
        for row in rec:
            protein = Protein(row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8],row[9],row[10],row[11])
            return protein
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table protein", error)
        #lib.log.info(source)
        lib.log.info(protein_id)

def insert_orthology(orthology):
    dbConnection = create_db_connection()
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO orthology VALUES (%s, %s)",
                                                                    (orthology.model_protein.id, 
                                                                    orthology.target_protein.id))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE orthology")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE orthology", error)
        lib.log.info(str(orthology.model_protein.id) + " " +
                    str(orthology.target_protein.id))
        
def parse_orthology_file(filename):
    lib.log.info("Parsing "+ filename)
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
                        #insert_orthology(orthology)
                    
    in_file.close()
    lib.log.info(filename + " parsed correctly")

def select_model_regulatory_by_organism_id(organism_id):
    model_regulatory_interactions = list()
    try:
        cursor = dbConnection.cursor()
        cursor.execute("SELECT DISTINCT * from model_regulatory model right join gene gen on model.tf_locus_tag = gen.locus_tag AND gen.organism = %s WHERE model.tf_locus_tag IS NOT NULL", (organism_id))
        rec = cursor.fetchall()
        for row in rec:
            model_regulatory = ModelRegulatory(row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8])
            model_regulatory_interactions.append(model_regulatory)
        cursor.close()
        return model_regulatory_interactions
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table model_regulatory", error)
        lib.log.info(organism_id)

def select_orthologs_by_target_organism(model_locus_tag, target_organism_id):
    orthology_list = list()
    try:
        cursor = dbConnection.cursor()
        cursor.execute("SELECT DISTINCT model_protein, target_protein from orthology ortho join protein prot on ortho.model_protein = prot.id join gene gen on prot.locus_tag = %s AND gen.organism = %s WHERE ortho.model_protein IS NOT NULL", (model_locus_tag, target_organism_id))
        rec = cursor.fetchall()
        for row in rec:
            ortho = Orthology(select_protein_by_id(row[0]),select_protein_by_id(row[1]))
            orthology_list.append(ortho)
        cursor.close()
        return orthology_list
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table orthology", error)
        lib.log.info(model_locus_tag)

def select_locus_tag_by_protein_id(protein_id):
    try:
        cursor = dbConnection.cursor()
        postgreSQL_select_Query = "SELECT locus_tag from protein WHERE id = %s"
        cursor.execute(postgreSQL_select_Query, (protein_id,))
        rec = cursor.fetchone()
        print(rec)
        cursor.close()
        return rec
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table protein", error)

def insert_regulatory_interaction(regulatory_interaction):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO regulatory_interaction VALUES (default, %s, %s, %s, %s)",
                                                                    (regulatory_interaction.tf_locus_tag, 
                                                                    regulatory_interaction.tg_locus_tag, 
                                                                    regulatory_interaction.regulatory_function,
                                                                    regulatory_interaction.pubmedid_source
                                                                    ))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE regulatory_interaction")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE regulatory_interaction", error)
        lib.log.info(str(regulatory_interaction.tf_locus_tag) + " " +
                    str(regulatory_interaction.tg_locus_tag) + " " +
                    str(regulatory_interaction.regulatory_function))
        
def construct_grn_orthology(model_organism_id, target_organism_id):
    
    model_regulatory_interactions = select_model_regulatory_by_organism_id(model_organism_id)
    
    for model_regulatory in model_regulatory_interactions:
        tf_orthologs = select_orthologs_by_target_organism(model_regulatory.tf_locus_tag, target_organism_id)
        tg_orthologs = select_orthologs_by_target_organism(model_regulatory.tg_locus_tag, target_organism_id)
        if len(tf_orthologs)!=0 and len(tg_orthologs)!=0:
            for ortholog_tf in tf_orthologs:
                for ortholog_tg in tg_orthologs:
                    regulatory_interaction = RegulatoryInteraction(0,ortholog_tf.target_protein.locus_tag, ortholog_tg.target_protein.locus_tag, model_regulatory.regulatory_function,model_regulatory.pubmedid)
                    
                    insert_regulatory_interaction(regulatory_interaction)
                    
                    #update gene as TF
                    tf = select_gene_by_locus_tag(ortholog_tf.target_protein.locus_tag)
                    tf.is_tf = 'True'
                    update_gene(tf)

def select_tfs_by_organism(organism_accession):
    dbConnection = create_db_connection()
    tf_list = list()
    try:
        cursor = dbConnection.cursor()
        cursor.execute("SELECT * from gene WHERE is_tf = 'True' and organism_accession = %s", (organism_accession,))
        rec = cursor.fetchall()
        #queryset = Gene.objects.filter(is_tf = 'true').filter(organism_accession = organism_accession)
        for row in rec:
            tf_list.append(Gene(row[0],row[1],row[2],row[3],row[4]))
        cursor.close()
        return tf_list
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table gene", error)
        lib.log.info(organism_accession)

def select_pwms_by_locus_tag(locus_tag):
    dbConnection = create_db_connection()
    pwm_list = list()
    try:
        cursor = dbConnection.cursor()
        cursor.execute("SELECT * from pwm WHERE locus_tag = %s", (locus_tag,))
        rec = cursor.fetchall()
        #rec = Pwm.objects.filter(locus_tag = locus_tag)
        
        for row in rec:
            pwm_list.append(Pwm(row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8]))
        cursor.close()
        return pwm_list
    except (Exception, psycopg2.Error) as error:
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

def select_promoter_by_locus_tag(tg_locus_tag):
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

def insert_tfbs_prediction(tfbs):
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


if __name__ == '__main__':
    
    #Create a database connection
    dbConnection = create_db_connection()
    
    if (dbConnection):
        dbConnection.close()
        lib.log.info("The DB connection is closed")