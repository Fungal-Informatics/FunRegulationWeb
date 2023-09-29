import os
import psycopg2
from collections import namedtuple
import sys, re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import urllib.parse
from OldModels.models import *
from time import sleep

gffInfoFields = ["seqid", "source", "ltype", "start", "end", "score", "strand", "phase", "attributes"]
GFFRecord = namedtuple("GFFRecord", gffInfoFields)
modelsOrganisms = ['GCA_000182925.2', 'GCA_000240135.3', 'R64-3-1-SGD', 's10-m04-r16-AspGD']
modelsAccessions = {'NEUROSPORA': 'GCA_000182925.2', 
                    'FUSARIUM': 'GCA_000240135.3', 
                    'SACCHAROMYCES':'R64-3-1-SGD',
                    'A.NIDULANS':'s10-m04-r16-AspGD'}

#GENES E PROMOTERS
#1 - in_file_genes = os.path.join('/home/gabriel/Downloads/ArquivosTCC/Database/Database/genbank/GCA_000182925.2/Neurospora_crassa.NC12.51.chr.gff3')
#2 - in_file_genes = os.path.join('/home/gabriel/Downloads/ArquivosTCC/Database/Database/genbank/GCA_000240135.3/Fusarium_graminearum_gca_000240135.ASM24013v3.51.gff3')
#3 - in_file_genes = os.path.join('/home/gabriel/Downloads/ArquivosTCC/Database/Database/genbank/R64-3-1-SGD/Saccharomyces_cerevisiae_gene_set_R64-3-1_20210421.gff3')
#4 - in_file_genes = os.path.join('/home/gabriel/Downloads/ArquivosTCC/Database/Database/genbank/s10-m04-r16-AspGD/A_nidulans_FGSC_A4_current_features.gff3')

#PROTEINS
#1 - in_file_proteins = os.path.join('/home/gabriel/Downloads/ArquivosTCC/Database/Database/genbank/GCA_000182925.2/Neurospora_crassa.NC12.pep.all.fa')
#2 - in_file_proteins = os.path.join('/home/gabriel/Downloads/ArquivosTCC/Database/Database/genbank/GCA_000240135.3/Fusarium_graminearum_gca_000240135.ASM24013v3.pep.all.fa')
#3 - in_file_proteins = os.path.join('/home/gabriel/Downloads/ArquivosTCC/Database/Database/genbank/R64-3-1-SGD/Saccharomyces_cerevisiae_proteins_R64-3-1_20210421.fasta')
#4 - in_file_proteins = os.path.join('/home/gabriel/Downloads/ArquivosTCC/Database/Database/genbank/s10-m04-r16-AspGD/A_nidulans_FGSC_A4_current_orf_trans_all.fasta')

#PWMS
#in_file_pwms = os.path.join('/home/gabriel/Downloads/ArquivosTCC/Database/Database/cis-bp/species/model/GCA_000182925.2/TF_Information.txt')

#MODELS_INTERACTIONS
#1 - in_file_model_regulatory = os.path.join('/home/gabriel/Downloads/ArquivosTCC/Database/Database/model-interactions/Aspergillus_nidulans.tsv')
#2 - in_file_model_regulatory = os.path.join('/home/gabriel/Downloads/ArquivosTCC/Database/Database/model-interactions/Fusarium_graminearum.tsv')
#3 - in_file_model_regulatory = os.path.join('/home/gabriel/Downloads/ArquivosTCC/Database/Database/model-interactions/Neurospora_crassa.tsv')
in_file_model_regulatory = os.path.join('/home/gabriel/Downloads/ArquivosTCC/Database/Database/model-interactions/Saccharomyces_cerevisiae.tsv')
dbConnection = None;
upstream = -1000
downstream = 0

def create_db_connection():
    try:
        con = psycopg2.connect(host='localhost', database='funregulationtcc',
        user='postgres', password='postgres')
        print("Successfully Connected to PostgreSQL")
        return con
    except (Exception, psycopg2.Error) as error:
        pass
        #lib.log.info(error)

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

def parse_gff3_file(filename):
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
    #lib.log.info("GFF3 File parsed correctly")

def insert_gene(gene):
    dbConnection = create_db_connection()
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO gene VALUES (%s, %s, %s, %s)",
                                                        (gene.locus_tag,
                                                         gene.symbol_gene, 
                                                         gene.is_tf, 
                                                         gene.organism_accession))
        dbConnection.commit()
        print("Record inserted successfully into TABLE gene")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        print(error)

def insert_promoter(promoter):
    dbConnection = create_db_connection()
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO promoter VALUES (%s, %s, %s, %s, %s, %s, %s)",
                                                            (promoter.locus_tag,
                                                            promoter.strand, 
                                                            promoter.source, 
                                                            promoter.start, 
                                                            promoter.stop, 
                                                            promoter.promoter_seq,
                                                            promoter.organism_accession))
        dbConnection.commit()
        print("Record inserted successfully into TABLE promoter")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        print("Failed to insert data into TABLE promoter", error)
         
def gff3_handler(in_file_genes):
    print("Parsing "+ in_file_genes)
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
            #organism_id = select_organism_by_assembly_name(source)
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
                if record.attributes.get("Gene") is not None:
                    symbol_gene = record.attributes.get("Gene")
                is_tf = False
                #print(symbol_gene)
                organsim_accession = 's10-m04-r16-AspGD'
                gene = Gene(organsim_accession,locus_tag,symbol_gene,is_tf)
                insert_gene(gene)
                
                promoter = None
                if record.strand == '+':
                    if record.start+upstream > 0 :
                        promoter = Promoter(organsim_accession,locus_tag, record.strand, record.seqid, record.start+upstream, record.start+downstream, None)
                        #promoter = Promoter(locus_tag=Gene.objects.get(locus_tag=locus_tag), strand=record.strand, source=record.seqid, start=record.start+upstream, stop=record.start+downstream, promoter_seq=None)
                    else:
                        # incomplete promoters
                        promoter = Promoter(organsim_accession,locus_tag, record.strand, record.seqid, 1, record.start+downstream, None)
                        #promoter = Promoter(locus_tag=Gene.objects.get(locus_tag=locus_tag), strand=record.strand, source=record.seqid, start=1, stop=record.start+downstream, promoter_seq=None)
                        print("Promoter of gene " + locus_tag + " can't be fully indentified")
                        promoters_partially_extracted += 1
                else:
                    if record.end-source_size <= 0 :
                        promoter = Promoter(organsim_accession,locus_tag, record.strand, record.seqid, record.end-upstream, record.end-downstream, None)
                        #promoter = Promoter(locus_tag=Gene.objects.get(locus_tag=locus_tag), strand=record.strand, source=record.seqid, start=record.end-upstream, stop=record.end-downstream,promoter_seq=None)
                    else:
                        # incomplete promoters
                        promoter = Promoter(organsim_accession,locus_tag, record.strand, record.seqid, source_size, record.end-downstream, None)
                        #promoter = Promoter(locus_tag=Gene.objects.get(locus_tag=locus_tag), strand=record.strand, source=record.seqid, start=source_size, stop=record.end-downstream, promoter_seq=None)
                        print("Promoter of gene " + locus_tag + " can't be fully indentified")
                        promoters_partially_extracted += 1
                recordCount += 1
                #promoter.save()
                insert_promoter(promoter)
        pos=pos+1

def insert_protein(protein):
    dbConnection = create_db_connection()
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO protein VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
                                                            (protein.id,
                                                            protein.product,
                                                            protein.interpro,
                                                            protein.pfam,
                                                            protein.go, 
                                                            protein.gene3d, 
                                                            protein.reactome,
                                                            protein.panther,
                                                            protein.uniprot,
                                                            protein.ec_number,
                                                            protein.cazy,
                                                            protein.locus_tag,
                                                            protein.organism_accession
                                                            ))
        dbConnection.commit()
        print("Record inserted successfully into TABLE protein")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        print("Failed to insert data into TABLE protein", error)
        print(str(protein.locus_tag) + " " +
                    str(protein.id))
        
def parse_protein_file(in_file_proteins, organism_accession):
    print("Parsing "+ in_file_proteins)
    for rec in SeqIO.parse(in_file_proteins, 'fasta'):
        
        #when locus_tag != protein_id
        #rec.description = re.search(r'gene:(.*?) transcript:', rec.description).group(1)
        #protein = Protein(organism_accession, rec.description, rec.id,'','','','','','','','','','')
        #print(rec.description, rec.id)

        #when locus_tag == protein_id
        protein = Protein(organism_accession, rec.id, rec.id,'','','','','','','','','','')
        #print(rec.id)
        insert_protein(protein)

def insert_pwm(pwm):
    dbConnection = create_db_connection()
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO pwm VALUES (default, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
                                                        (pwm.motif_id,
                                                        pwm.status, 
                                                        pwm.tf_family, 
                                                        pwm.motif_type, 
                                                        pwm.msource_author, 
                                                        pwm.msource,
                                                        pwm.pubmedid,
                                                        pwm.locus_tag,
                                                        pwm.organism_accession
                                                        ))
        dbConnection.commit()
        print("Record inserted successfully into TABLE pwm")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        print("Failed to insert data into TABLE pwm", error)
        print(str(pwm.organism_accession) + " " +
                    str(pwm.locus_tag) + " " +
                    str(pwm.motif_id) + " " +
                    str(pwm.status) + " " +
                    str(pwm.tf_family) + " " +
                    str(pwm.motif_type) + " " +
                    str(pwm.msource_author) + " " +
                    str(pwm.msource) + " " +
                    str(pwm.pubmedid))
        
def parse_pwm_file(organism_accession):
    in_file_pwm = os.path.join('/home/gabriel/Downloads/ArquivosTCC/Database/Database/cis-bp/species/model/'+organism_accession+'/TF_Information.txt')
    print("Parsing "+ in_file_pwm)
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
                pwm = Pwm(0, organism_accession, locus_tag, motif_id, status, tf_family, motif_type, msource_author, msource, pubmedid)
                insert_pwm(pwm)
    in_file.close()
    print(in_file_pwm + " parsed correctly")

def select_gene_by_locus_tag(locus_tag):
    dbConnection = create_db_connection()
    gene = None
    try:
        cursor = dbConnection.cursor()
        postgreSQL_select_Query = "SELECT * FROM gene WHERE locus_tag = %s"
        cursor.execute(postgreSQL_select_Query, (locus_tag,))
        rec = cursor.fetchall()
        for row in rec:
            gene = Gene(locus_tag=row[0],symbol_gene=row[1],is_tf=row[2],organism_accession=row[3])
            return gene
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        print("Failed to execute the select into table gene", error)

def update_gene(gene):
    dbConnection = create_db_connection()
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("UPDATE gene SET symbol_gene = %s, is_tf = %s WHERE organism_accession = %s AND locus_tag = %s",
                                (gene.symbol_gene,
                                gene.is_tf,
                                gene.organism_accession,
                                gene.locus_tag))
        dbConnection.commit()
        print("Record updated successfully into TABLE gene")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        print("Failed to update data into TABLE gene", error)
        print(str(gene.organism_accession) + " " +
                     str(gene.locus_tag) + " " + 
                     str(gene.symbol_gene) + " " +
                     str(gene.is_tf))

def insert_model_regulatory(model_regulatory):
    dbConnection = create_db_connection()
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO model_regulatory VALUES (default, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
                                                                    (model_regulatory.regulatory_function,
                                                                    model_regulatory.evidence, 
                                                                    model_regulatory.experiment, 
                                                                    model_regulatory.experimental_condition, 
                                                                    model_regulatory.pubmedid, 
                                                                    model_regulatory.publication,
                                                                    model_regulatory.organism_accession,
                                                                    model_regulatory.tf_locus_tag,
                                                                    model_regulatory.tg_locus_tag
                                                                    ))
        dbConnection.commit()
        print("Record inserted successfully into TABLE model_regulatory")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        print("Failed to insert data into TABLE model_regulatory", error)
        print(str(model_regulatory.tf_locus_tag) + " " +
                    str(model_regulatory.tg_locus_tag) + " " +
                    str(model_regulatory.regulatory_function))
            
def parse_model_regulatory_file(filename, organism_accession):
    print("Parsing " + filename)
    with open(filename) as in_file:
        for line in in_file:
            if line.startswith("#"): 
                continue
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
                model_regulatory = ModelRegulatory(0, organism_accession, tf_locus_tag,tg_locus_tag,regulatory_function,evidence,experiment,experimental_condition,pubmedid,publication)
                insert_model_regulatory(model_regulatory)
    in_file.close()
    print(filename + " parsed correctly")

#gff3_handler(in_file_genes)
#parse_protein_file(in_file_proteins, 's10-m04-r16-AspGD')
# for organism in modelsOrganisms:
#     parse_pwm_file(organism)
#parse_model_regulatory_file(in_file_model_regulatory, 's10-m04-r16-AspGD')
#parse_model_regulatory_file(in_file_model_regulatory, 'GCA_000240135.3')
#parse_model_regulatory_file(in_file_model_regulatory, 'GCA_000182925.2')
parse_model_regulatory_file(in_file_model_regulatory, 'R64-3-1-SGD')