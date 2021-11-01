"""
    Insert Orthology
"""
def insert_orthology(orthology):
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