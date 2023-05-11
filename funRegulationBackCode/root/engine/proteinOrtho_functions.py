import psycopg2
import root.lib.library as lib
from api.models import *
import os
import sys, re
import urllib.parse
from collections import namedtuple
from Bio import SeqIO
from django.conf import settings
from BCBio import GFF
from time import sleep

upstream = -1000
downstream = 0

#create log file
log_name = os.path.join(settings.LOG_FILE_PATH)
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)

gffInfoFields = ["seqid", "source", "ltype", "start", "end", "score", "strand", "phase", "attributes"]
GFFRecord = namedtuple("GFFRecord", gffInfoFields)

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

def select_organism_by_assembly_name(source):
    dbConnection = create_db_connection()
    organism = 0
    try:
        cursor = dbConnection.cursor()
        postgreSQL_select_Query = "SELECT accession FROM api_organism WHERE accession = %s"
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
    dbConnection = create_db_connection()
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO gene VALUES (%s, %s, %s, %s, %s)",
                                                        (gene.organism_accession.accession, 
                                                        gene.locus_tag, 
                                                        gene.symbol_gene, 
                                                        gene.description, 
                                                        gene.is_tf))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE gene")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE gene", error)
        lib.log.info(str(gene.organism_accession.accession) + " " +
                     str(gene.locus_tag) + " " + 
                     str(gene.symbol_gene) + " " +
                     str(gene.description) + " " +
                     str(gene.is_tf))

def insert_promoter(promoter):
    dbConnection = create_db_connection()
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
                #locus_tag = record.attributes.get("ID")
                locus_tag = record.attributes.get("gene_id")
                description = record.attributes.get("description")
                #description = None #record.attributes.get("description")
                symbol_gene = ''
                if record.attributes.get("gene") is not None:
                    symbol_gene = record.attributes.get("gene")
                is_tf = False
                gene = Gene(organism_id,locus_tag,symbol_gene,description,is_tf)
                print(gene)
                #insert_gene(gene)

                # promoter = None
                # if record.strand == '+':
                #     if record.start+upstream > 0 :
                #         promoter = Promoter(locus_tag, record.strand, record.seqid, record.start+upstream, record.start+downstream)
                #     else:
                #         # incomplete promoters
                #         promoter = Promoter(locus_tag, record.strand, record.seqid, 1, record.start+downstream)
                #         lib.log.info("Promoter of gene " + locus_tag + " can't be fully indentified")
                #         promoters_partially_extracted += 1
                # else:
                #     if record.end-source_size <= 0 :
                #         promoter = Promoter(locus_tag, record.strand, record.seqid, record.end-upstream, record.end-downstream)
                #     else:
                #         # incomplete promoters
                #         promoter = Promoter(locus_tag, record.strand, record.seqid, source_size, record.end-downstream)
                #         lib.log.info("Promoter of gene " + locus_tag + " can't be fully indentified")
                #         promoters_partially_extracted += 1
                # recordCount += 1
                #insert_promoter(promoter)
        pos=pos+1
    
    # lib.log.info("%d genes were found" % recordCount)
    # lib.log.info("Promoters partially identified: %d" % promoters_partially_extracted)
    # lib.log.info("GFF3 file successfully parsed")

def gff3_handler2(in_file_genes, organism_accession):
    lib.log.info("Parsing "+ in_file_genes)
    recordCount = 0
    promoters_partially_extracted = 0
    organism_id = 0
    
    record_list = list()
    for record in parse_gff3_file(in_file_genes):
        record_list.append(record)
    
    pos = 0
    organism_id = organism_accession
    
    while (pos<len(record_list)):
        record = record_list[pos]
        if record.ltype == 'chromosome' or record.ltype == 'supercontig' or record.ltype == 'region':
            source_size = record.end
        elif (record.ltype == 'gene' or record.ltype == 'pseudogene' or 
              record.ltype == 'transposable_element_gene' or 
              record.ltype == 'blocked_reading_frame'):
            #Access attributes like this: my_strand = record.strand
            #gene = str(record.attributes)
            locus_tag = record.attributes.get("locus_tag")
            description = record.attributes.get("product")
            symbol_gene = ''
            if record.attributes.get("gene") is not None:
                symbol_gene = record.attributes.get("gene")
            is_tf = False
            
            gene = Gene(organism_accession=Organism.objects.get(accession=organism_id),locus_tag=locus_tag,symbol_gene=symbol_gene,description=description,is_tf=is_tf)
            gene.save()

            promoter = None
            if record.strand == '+':
                if record.start+upstream > 0 :
                    promoter = Promoter(locus_tag=Gene.objects.get(locus_tag=locus_tag), strand=record.strand, source=record.seqid, start=record.start+upstream, stop=record.start+downstream)
                else:
                    # incomplete promoters
                    promoter = Promoter(locus_tag=Gene.objects.get(locus_tag=locus_tag), strand=record.strand, source=record.seqid, start=1, stop=record.start+downstream)
                    lib.log.info("Promoter of gene " + locus_tag + " can't be fully indentified")
                    promoters_partially_extracted += 1
            else:
                if record.end-source_size <= 0 :
                    promoter = Promoter(locus_tag=Gene.objects.get(locus_tag=locus_tag), strand=record.strand, source=record.seqid, start=record.end-upstream, stop=record.end-downstream)
                else:
                    # incomplete promoters
                    promoter = Promoter(locus_tag=Gene.objects.get(locus_tag=locus_tag), strand=record.strand, source=record.seqid, start=source_size, stop=record.end-downstream)
                    lib.log.info("Promoter of gene " + locus_tag + " can't be fully indentified")
                    promoters_partially_extracted += 1
            recordCount += 1
            
            promoter.save()
        pos=pos+1
    
    lib.log.info("%d genes were found" % recordCount)
    lib.log.info("Promoters partially identified: %d" % promoters_partially_extracted)
    lib.log.info("GFF3 file successfully parsed")

def insert_protein(protein):
    dbConnection = create_db_connection()
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
        #rec.description = re.search(r'gene:(.*?) transcript:', rec.description).group(1)
        locus_tag = rec.description.split()[3].strip()
        
        #when locus_tag == protein_id
        #protein = Protein(rec.id, rec.id,'','','','','','','','','','')
        protein = Protein(locus_tag=Gene.objects.get(locus_tag=locus_tag), id=rec.id,interpro='',pfam='',go='',gene3d='',reactome='',panther='',uniprot='',kegg_enzyme='',cazy='',uniparc='')
        protein.save()
        
def select_protein_by_id(protein_id):
    protein = Protein.objects.filter(id=protein_id).values('locus_tag')
    
    if(len(protein) > 0):
        for locus_tag_value in protein:
            locus_tag = locus_tag_value['locus_tag']

    return locus_tag
    # dbConnection = create_db_connection()
    # protein = None
    # try:
    #     cursor = dbConnection.cursor()
    #     postgreSQL_select_Query = "SELECT * FROM protein WHERE id = %s"
    #     cursor.execute(postgreSQL_select_Query, (protein_id,))
    #     rec = cursor.fetchall()
    #     for row in rec:
    #         protein = Protein(row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8],row[9],row[10],row[11])
    #         return protein
    #     cursor.close()
    # except (Exception, psycopg2.Error) as error:
    #     lib.log.info("Failed to execute the select into table protein", error)
    #     #lib.log.info(source)
    #     lib.log.info(protein_id)

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