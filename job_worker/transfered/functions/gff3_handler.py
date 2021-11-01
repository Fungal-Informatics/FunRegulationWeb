"""
    GFF3 handler
"""
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
        else:
            if (record.ltype == 'gene' or 
                record.ltype == 'pseudogene' or 
                record.ltype == 'transposable_element_gene' or 
                record.ltype == 'blocked_reading_frame'):
                #Access attributes like this: my_strand = record.strand
                #gene = str(record.attributes)
                locus_tag = record.attributes.get("ID")
                #description = record.attributes.get("description")
                description = None #record.attributes.get("description")
                symbol_gene = ''
                if record.attributes.get("Name") is not None:
                    symbol_gene = record.attributes.get("Name")
                is_tf = False
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
