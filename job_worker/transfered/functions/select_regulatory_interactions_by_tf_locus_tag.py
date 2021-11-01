"""
    Select regulatory interactions by tf_locus_tag
"""
def select_regulatory_interactions_by_tf_locus_tag(tf_locus_tag):
    regulatory_interactions = list()
    try:
        cursor = dbConnection.cursor()
        postgreSQL_select_Query = "SELECT * from regulatory_interaction WHERE tf_locus_tag = %s"
        cursor.execute(postgreSQL_select_Query, (tf_locus_tag,))
        rec = cursor.fetchall()
        for row in rec:
            regulatory_interactions.append(RegulatoryInteraction(row[0],row[1],row[2],row[3],row[4]))
        cursor.close()
        return regulatory_interactions
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table regulatory_interaction", error)
        lib.log.info(tf_locus_tag)