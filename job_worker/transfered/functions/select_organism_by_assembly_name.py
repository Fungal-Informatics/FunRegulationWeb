"""
    Select organism by assembly
"""
def select_organism_by_assembly_name(source):
    organism = 0
    try:
        cursor = dbConnection.cursor()
        postgreSQL_select_Query = "SELECT id FROM organism WHERE assembly_name = %s"
        cursor.execute(postgreSQL_select_Query, (source,))
        records = cursor.fetchall()
        for row in records:
            organism = row[0]
        return organism
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table organism", error)
        lib.log.info(source)