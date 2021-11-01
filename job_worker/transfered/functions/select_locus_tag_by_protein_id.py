"""
    Select locus_tag by protein id
"""
def select_locus_tag_by_protein_id(protein_id):
    try:
        cursor = dbConnection.cursor()
        postgreSQL_select_Query = "SELECT locus_tag from protein WHERE id = %s"
        cursor.execute(postgreSQL_select_Query, (protein_id,))
        rec = cursor.fetchone()
        print(rec)
        cursor.close()
        return rec
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table protein", error)