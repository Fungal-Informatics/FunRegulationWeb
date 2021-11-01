"""
    Select protein by id
"""
def select_protein_by_id(protein_id):
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
        lib.log.info(source)