"""
    Insert protein
"""
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
