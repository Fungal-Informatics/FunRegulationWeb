"""
    Insert organism
"""
def insert_organism(organism):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO organism VALUES (default, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
                                                            (organism.order, 
                                                            organism.genus, 
                                                            organism.species, 
                                                            organism.strain, 
                                                            organism.taxon, 
                                                            organism.assembly_name, 
                                                            organism.assembly_accession, 
                                                            organism.is_model, 
                                                            organism.cis_bp))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE organism")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE organism ")
        lib.log.info(str(organism.order) + " " +
                    str(organism.genus) + " " + 
                    str(organism.species) + " " +
                    str(organism.strain))
