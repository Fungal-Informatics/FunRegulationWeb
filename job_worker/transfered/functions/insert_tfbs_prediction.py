"""
    Insert TFBS prediction
"""
def insert_tfbs_prediction(tfbs):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO tfbs VALUES (default, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
                                                            (tfbs.regulatory_interaction_id, 
                                                            tfbs.pwm_id, 
                                                            tfbs.strand, 
                                                            tfbs.start, 
                                                            tfbs.end, 
                                                            tfbs.sequence, 
                                                            tfbs.weight, 
                                                            tfbs.pval, 
                                                            tfbs.ln_pval, 
                                                            tfbs.sig))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE tfbs")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE tfbs ")
        lib.log.info(str(tfbs.regulatory_interaction_id) + " " +
                    str(tfbs.pwm_id) + " " + 
                    str(tfbs.strand) + " " +
                    str(tfbs.sequence))