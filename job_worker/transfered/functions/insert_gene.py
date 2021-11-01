"""
    Insert gene
"""
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