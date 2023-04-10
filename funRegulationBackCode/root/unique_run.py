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
from modelsLenz.models import *

#Initialize Folder Paths
in_folder = ""
out_folder = ""

#Initialize file Paths
in_file_organisms = os.path.join("/home/gabriel/Desktop/FunRegulationBack-end/funregulation/FunRegulationWeb/funRegulationBackCode/funRegulation/projects/Database/Database/genbank_dataset_04-10-22.tsv")
in_file_genes = os.path.join('/home/gabriel/Downloads/ArquivosTCC/Database/Database/genbank/R64-3-1-SGD/Saccharomyces_cerevisiae_gene_set_R64-3-1_20210421.gff3')
in_file_proteins = os.path.join('/home/gabriel/Downloads/ArquivosTCC/Database/Database/genbank/GCA_000182925.2/Neurospora_crassa.NC12111.pep.all.fa')
in_file_pwms = os.path.join('/home/gabriel/Downloads/ArquivosTCC/Database/Database/cis-bp/species/model/s10-m04-r16-AspGD/TF_Information.txt')
in_file_model_regulatory = os.path.join('/home/gabriel/Downloads/ArquivosTCC/Database/Database/model-interactions/Saccharomyces_cerevisiae.tsv')
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

def create_db_connection():
    try:
        con = psycopg2.connect(host='localhost', database='funregulationtcc',
        user='postgres', password='postgres')
        lib.log.info("Successfully Connected to PostgreSQL")
        return con
    except (Exception, psycopg2.Error) as error:
        lib.log.info(error)
        

def insert_organism(organism):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO organism VALUES (%s, %s, %s, %s, %s, %s, %s)",
                                                            (
                                                            organism.accession,
                                                            organism.order, 
                                                            organism.genus, 
                                                            organism.species, 
                                                            organism.strain, 
                                                            organism.is_model, 
                                                            organism.cis_bp))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE organism")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE organism ")
        lib.log.info(error)
        lib.log.info(str(organism.order) + " " +
                    str(organism.genus) + " " + 
                    str(organism.species) + " " +
                    str(organism.strain) + " " +
                    str(organism.is_model) + " " +
                    str(organism.cis_bp))

"""
    Parse Ensembl Species input file
"""
def parse_organism_file(filename):
    lib.log.info("Parsing "+ filename)
    with open(filename) as in_file:
        for line in in_file:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            accession = urllib.parse.unquote(parts[0])
            order = urllib.parse.unquote(parts[1])
            genus = urllib.parse.unquote(parts[2])
            species = urllib.parse.unquote(parts[3])
            strain = urllib.parse.unquote(parts[4])
            is_model = urllib.parse.unquote(parts[5])
            cis_bp = urllib.parse.unquote(parts[6])
            organism = Organism(accession,order,genus,species,strain,is_model,cis_bp)
            insert_organism(organism)
    in_file.close()
    lib.log.info(filename + " parsed correctly")

"""
    Parse the GFF3 attribute column and return a dict
"""
def parse_gff_attributes(attributeString):
    if attributeString == ".": return {}
    ret ={}
    if ";" not in attributeString: 
        ret=attributeString
    if ";" in attributeString:
        for attribute in attributeString.split(";"):
            key, value = attribute.split("=")
            ret[urllib.parse.unquote(key)] = urllib.parse.unquote(value)
    return ret

"""
    A minimalistic GFF3 format parser.
    Yields objects that contain info about a single GFF3 feature.
"""
def parse_gff3_file(filename):
    lib.log.info("Parsing "+ filename)
    with open(filename) as infile:
        for line in infile:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            #If this fails, the file format is not standard-compatible
            assert len(parts) == len(gffInfoFields)
            #Normalize data
            normalizedInfo = {
                "seqid": None if parts[0] == "." else urllib.parse.unquote(parts[0]),
                "source": None if parts[1] == "." else urllib.parse.unquote(parts[1]),
                "ltype": None if parts[2] == "." else urllib.parse.unquote(parts[2]),
                "start": None if parts[3] == "." else int(parts[3]),
                "end": None if parts[4] == "." else int(parts[4]),
                "score": None if parts[5] == "." else float(parts[5]),
                "strand": None if parts[6] == "." else urllib.parse.unquote(parts[6]),
                "phase": None if parts[7] == "." else urllib.parse.unquote(parts[7]),
                "attributes": parse_gff_attributes(parts[8])
            }
            #Alternatively, you can emit the dictionary here, if you need mutability:
            #    yield normalizedInfo
            yield GFFRecord(**normalizedInfo)
    lib.log.info("GFF3 File parsed correctly")

def select_organism_by_assembly_name(source):
    organism = 0
    try:
        cursor = dbConnection.cursor()
        postgreSQL_select_Query = "SELECT accession FROM organism WHERE accession = %s"
        cursor.execute(postgreSQL_select_Query, (source,))
        records = cursor.fetchall()
        for row in records:
            organism = row[0]
        return organism
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table organism", error)
        lib.log.info(source)

def insert_gene(gene):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO gene VALUES (%s, %s, %s, %s, %s)",
                                                        (gene.organism, 
                                                        gene.locus_tag, 
                                                        gene.symbol_gene, 
                                                        gene.description, 
                                                        gene.is_tf))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE gene")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE gene", error)
        lib.log.info(str(gene.organism) + " " +
                     str(gene.locus_tag) + " " + 
                     str(gene.symbol_gene) + " " +
                     str(gene.description) + " " +
                     str(gene.is_tf))
                     
def update_gene(gene):
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

def insert_promoter(promoter):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO promoter VALUES (%s, %s, %s, %s, %s)",
                                                            (promoter.locus_tag, 
                                                            promoter.strand, 
                                                            promoter.source, 
                                                            promoter.start, 
                                                            promoter.stop))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE promoter")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE promoter", error)
        lib.log.info(str(promoter.locus_tag) + " " +
                    str(promoter.strand) + " " + 
                    str(promoter.source) + " " +
                    str(promoter.start) + " " +
                    str(promoter.stop))

def gff3_handler(in_file_genes):
    lib.log.info("Parsing "+ in_file_genes)
    recordCount = 0
    promoters_partially_extracted = 0
    organism_id = 0
    
    record_list = list()
    for record in parse_gff3_file(in_file_genes):
        record_list.append(record)
        
    pos = 0
    while (pos<len(record_list)):
        record = record_list[pos]
        if record.ltype == 'chromosome' or record.ltype == 'supercontig':
            source_size = record.end
            source = record.source
            organism_id = select_organism_by_assembly_name(source)
            #organism_id = select_organism_by_assembly_name("R64-3-1-SGD")
        else:
            if (record.ltype == 'gene' or 
                record.ltype == 'pseudogene' or 
                record.ltype == 'transposable_element_gene' or 
                record.ltype == 'blocked_reading_frame'):
                #Access attributes like this: my_strand = record.strand
                #gene = str(record.attributes)
                locus_tag = record.attributes.get("ID")
                #locus_tag = record.attributes.get("gene_id")
                #description = record.attributes.get("description")
                description = None #record.attributes.get("description")
                symbol_gene = ''
                if record.attributes.get("gene") is not None:
                    symbol_gene = record.attributes.get("gene")
                is_tf = False
                print(organism_id,locus_tag,symbol_gene,description,is_tf)
                gene = Gene(organism_id,locus_tag,symbol_gene,description,is_tf)
                insert_gene(gene)

                promoter = None
                if record.strand == '+':
                    if record.start+upstream > 0 :
                        promoter = Promoter(locus_tag, record.strand, record.seqid, record.start+upstream, record.start+downstream)
                    else:
                        # incomplete promoters
                        promoter = Promoter(locus_tag, record.strand, record.seqid, 1, record.start+downstream)
                        lib.log.info("Promoter of gene " + locus_tag + " can't be fully indentified")
                        promoters_partially_extracted += 1
                else:
                    if record.end-source_size <= 0 :
                        promoter = Promoter(locus_tag, record.strand, record.seqid, record.end-upstream, record.end-downstream)
                    else:
                        # incomplete promoters
                        promoter = Promoter(locus_tag, record.strand, record.seqid, source_size, record.end-downstream)
                        lib.log.info("Promoter of gene " + locus_tag + " can't be fully indentified")
                        promoters_partially_extracted += 1
                recordCount += 1
                insert_promoter(promoter)
        pos=pos+1
    
    lib.log.info("%d genes were found" % recordCount)
    lib.log.info("Promoters partially identified: %d" % promoters_partially_extracted)
    lib.log.info("GFF3 file successfully parsed")

def insert_protein(protein):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO protein VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
                                                            (protein.locus_tag, 
                                                            protein.id, 
                                                            protein.interpro, 
                                                            protein.pfam, 
                                                            protein.go,
                                                            protein.gene3d,
                                                            protein.reactome,
                                                            protein.panther,
                                                            protein.uniprot,
                                                            protein.kegg_enzyme,
                                                            protein.cazy,
                                                            protein.uniparc
                                                            ))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE protein")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE protein", error)
        lib.log.info(str(protein.locus_tag) + " " +
                    str(protein.id))

def parse_protein_file(in_file_proteins):
    lib.log.info("Parsing "+ in_file_proteins)
    for rec in SeqIO.parse(in_file_proteins, 'fasta'):
        
        #when locus_tag != protein_id
        rec.description = re.search(r'gene:(.*?) transcript:', rec.description).group(1)
        protein = Protein(rec.description, rec.id,'','','','','','','','','','')
        
        #when locus_tag == protein_id
        #protein = Protein(rec.id, rec.id,'','','','','','','','','','')
        #print(protein.locus_tag, protein.id)
        insert_protein(protein)

def insert_pwm(pwm):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO pwm VALUES (default, %s, %s, %s, %s, %s, %s, %s, %s)",
                                                        (pwm.locus_tag, 
                                                        pwm.motif_id, 
                                                        pwm.status, 
                                                        pwm.tf_family, 
                                                        pwm.motif_type,
                                                        pwm.msource_author,
                                                        pwm.msource,
                                                        pwm.pubmedid
                                                        ))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE pwm")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE pwm", error)
        lib.log.info(str(pwm.locus_tag) + " " +
                    str(pwm.motif_id) + " " +
                    str(pwm.status) + " " +
                    str(pwm.tf_family) + " " +
                    str(pwm.motif_type) + " " +
                    str(pwm.msource_author) + " " +
                    str(pwm.msource) + " " +
                    str(pwm.pubmedid))

def parse_pwm_file(in_file_pwm):
    lib.log.info("Parsing "+ in_file_pwm)
    with open(in_file_pwm) as in_file:
        for line in in_file:
            if line.startswith("#"): continue
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
                pwm = Pwm(0,locus_tag, motif_id, status, tf_family, motif_type, msource_author, msource, pubmedid)
                #print(pwm.locus_tag, pwm.motif_id, pwm.status, pwm.tf_family, pwm.motif_type, pwm.msource_author, pwm.msource, pwm.pubmedid)
                insert_pwm(pwm)
    in_file.close()
    lib.log.info(in_file_pwm + " parsed correctly")

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

def parse_model_regulatory_file(filename):
    lib.log.info("Parsing " + filename)
    with open(filename) as in_file:
        for line in in_file:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            
            #update gene as TF
            tf_locus_tag = urllib.parse.unquote(parts[0])
            tf = select_gene_by_locus_tag(tf_locus_tag)
            tf.is_tf = 'True'
            if tf.symbol_gene is None or tf.symbol_gene == '':
                tf.symbol_gene = urllib.parse.unquote(parts[1])
            update_gene(tf)
            
            #update gene as TG
            tg_locus_tag = urllib.parse.unquote(parts[2])
            tg = select_gene_by_locus_tag(tg_locus_tag)
            if tg is not None:
                if tg.symbol_gene is None or tg.symbol_gene == '':
                    tg.symbol_gene = urllib.parse.unquote(parts[3])
                update_gene(tg)
                regulatory_function = urllib.parse.unquote(parts[4])
                evidence = urllib.parse.unquote(parts[5])
                experiment = urllib.parse.unquote(parts[6])
                experimental_condition = urllib.parse.unquote(parts[7])
                pubmedid = urllib.parse.unquote(parts[8])
                publication = urllib.parse.unquote(parts[9])
                mr = ModelRegulatory(0,tf_locus_tag,tg_locus_tag,regulatory_function,evidence,experiment,experimental_condition,pubmedid,publication)
                #print(mr.tf_locus_tag,mr.tg_locus_tag,mr.regulatory_function,mr.evidence,mr.experiment,mr.experimental_condition,mr.pubmedid,mr.publication)
                insert_model_regulatory(mr)
    in_file.close()
    lib.log.info(filename + " parsed correctly")

def insert_model_regulatory(model_regulatory):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO model_regulatory VALUES (default, %s, %s, %s, %s, %s, %s, %s, %s)",
                                                                    (model_regulatory.tf_locus_tag, 
                                                                    model_regulatory.tg_locus_tag, 
                                                                    model_regulatory.regulatory_function, 
                                                                    model_regulatory.evidence, 
                                                                    model_regulatory.experiment,
                                                                    model_regulatory.experimental_condition,
                                                                    model_regulatory.pubmedid,
                                                                    model_regulatory.publication
                                                                    ))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE model_regulatory")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE model_regulatory", error)
        lib.log.info(str(model_regulatory.tf_locus_tag) + " " +
                    str(model_regulatory.tg_locus_tag) + " " +
                    str(model_regulatory.regulatory_function))               


"""
    Main function of this program
    This program should run just "one time", to populate the database with the model organisms
""" 
if __name__ == '__main__':
    
    #Create a database connection
    dbConnection = create_db_connection()
    
    #Load Input Files
    #1-parse_organism_file(in_file_organisms)
    #2-gff3_handler(in_file_genes)
    #3-parse_protein_file(in_file_proteins)
    #4-parse_pwm_file(in_file_pwms)
    #5-parse_model_regulatory_file(in_file_model_regulatory)
    
    if (dbConnection):
        dbConnection.close()
        lib.log.info("The DB connection is closed")