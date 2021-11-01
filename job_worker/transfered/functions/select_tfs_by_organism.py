"""
    Select tfs by organism id
"""
def select_tfs_by_organism(organism_id):
    tf_list = list()
    try:
        cursor = dbConnection.cursor()
        cursor.execute("SELECT * from gene WHERE is_tf = 'True' and organism = %s", (organism_id))
        rec = cursor.fetchall()
        for row in rec:
            tf_list.append(Gene(row[0],row[1],row[2],row[3],row[4]))
        cursor.close()
        return tf_list
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table gene", error)
        lib.log.info(organism_id)