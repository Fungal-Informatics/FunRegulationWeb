"""
    Select model regulatory interactions by organism id
"""
def select_model_regulatory_by_organism_id(organism_id):
    model_regulatory_interactions = list()
    try:
        cursor = dbConnection.cursor()
        cursor.execute("SELECT DISTINCT * from model_regulatory model right join gene gen on model.tf_locus_tag = gen.locus_tag AND gen.organism = %s WHERE model.tf_locus_tag IS NOT NULL", (organism_id))
        rec = cursor.fetchall()
        for row in rec:
            model_regulatory = ModelRegulatory(row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8])
            model_regulatory_interactions.append(model_regulatory)
        cursor.close()
        return model_regulatory_interactions
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table model_regulatory", error)
        lib.log.info(organism_id)