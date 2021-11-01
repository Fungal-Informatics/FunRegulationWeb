"""
    Select gene by locus_tag
"""
def select_gene_by_locus_tag(locus_tag):
    gene = None
    try:
        cursor = dbConnection.cursor()
        postgreSQL_select_Query = "SELECT * FROM gene WHERE locus_tag = %s"
        cursor.execute(postgreSQL_select_Query, (locus_tag,))
        rec = cursor.fetchall()
        for row in rec:
            gene = Gene(row[0],row[1],row[2],row[3],row[4])
            return gene
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table gene", error)

