"""
    Parse Protein Fasta file
"""
def parse_protein_file(in_file_proteins):
    lib.log.info("Parsing "+ in_file_proteins)
    for rec in SeqIO.parse(in_file_proteins, 'fasta'):
        
        #when locus_tag != protein_id
        #rec.description = re.search(r'gene:(.*?) transcript:', rec.description).group(1)
        #protein = Protein(rec.description, rec.id,'','','','','','','','','','')
        
        #when locus_tag == protein_id
        protein = Protein(rec.id, rec.id,'','','','','','','','','','')
        insert_protein(protein)