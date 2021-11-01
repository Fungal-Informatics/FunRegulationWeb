"""
    Insert Model regulatory interaction
"""
def insert_model_regulatory(model_regulatory):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO model_regulatory VALUES (default, %s, %s, %s, %s, %s, %s, %s, %s)",
                                                                    (model_regulatory.tf_locus_tag, 
                                                                    model_regulatory.tg_locus_tag, 
                                                                    model_regulatory.regulatory_function, 
                                                                    model_regulatory.evidence, 
                                                                    model_regulatory.experiment,
                                                                    model_regulatory.experimental_condition,
                                                                    model_regulatory.pubmedid,
                                                                    model_regulatory.publication
                                                                    ))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE model_regulatory")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE model_regulatory", error)
        lib.log.info(str(model_regulatory.tf_locus_tag) + " " +
                    str(model_regulatory.tg_locus_tag) + " " +
                    str(model_regulatory.regulatory_function))