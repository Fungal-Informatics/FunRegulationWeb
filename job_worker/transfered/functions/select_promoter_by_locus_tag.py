"""
    Select promoter by tg_locus_tag
"""
def select_promoter_by_locus_tag(tg_locus_tag):
    try:
        cursor = dbConnection.cursor()
        postgreSQL_select_Query = "SELECT * from promoter WHERE locus_tag = %s"
        cursor.execute(postgreSQL_select_Query, (tg_locus_tag,))
        rec = cursor.fetchone()
        promoter = Promoter(rec[0],rec[1],rec[2],rec[3],rec[4])
        cursor.close()
        return promoter
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table promoter", error)
        lib.log.info(tg_locus_tag)