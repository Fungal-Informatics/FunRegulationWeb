import psycopg2
import root.lib.library as lib
from api.models import *
import os
import sys
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

def select_model_regulatory_by_organism_id(organism_id):
    dbConnection = create_db_connection()
    model_regulatory_interactions = list()
    try:
        cursor = dbConnection.cursor()
        cursor.execute("SELECT DISTINCT * from model_regulatory model right join gene gen on model.tf_locus_tag = gen.locus_tag AND gen.organism_accession = %s WHERE model.tf_locus_tag IS NOT NULL", (organism_id))
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
    dbConnection = create_db_connection()
    orthology_list = list()
    try:
        cursor = dbConnection.cursor()
        cursor.execute("SELECT DISTINCT model_protein, target_protein from orthology ortho join protein prot on ortho.model_protein = prot.id join gene gen on prot.locus_tag = %s AND gen.organism_accession = %s WHERE ortho.model_protein IS NOT NULL", (model_locus_tag, target_organism_id))
        rec = cursor.fetchall()
        for row in rec:
            ortho = Orthology(select_protein_by_id(row[0]),select_protein_by_id(row[1]))
            orthology_list.append(ortho)
        cursor.close()
        return orthology_list
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table orthology", error)
        lib.log.info(model_locus_tag)

def insert_regulatory_interaction(regulatory_interaction):
    dbConnection = create_db_connection()
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
        
def select_gene_by_locus_tag(locus_tag):
    dbConnection = create_db_connection()
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

def update_gene(gene):
    dbConnection = create_db_connection()
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("UPDATE gene SET symbol_gene = %s, description = %s, is_tf = %s WHERE organism_accession = %s AND locus_tag = %s",
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
        
def create_db_connection():
    try:
        con = psycopg2.connect(host='localhost', database='funregulationtcc',
        user='postgres', password='postgres')
        lib.log.info("Successfully Connected to PostgreSQL")
        return con
    except (Exception, psycopg2.Error) as error:
        lib.log.info(error)