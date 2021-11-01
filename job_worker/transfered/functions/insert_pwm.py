"""
    Insert pwm
"""
def insert_pwm(pwm):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO pwm VALUES (default, %s, %s, %s, %s, %s, %s, %s, %s)",
                                                        (pwm.locus_tag, 
                                                        pwm.motif_id, 
                                                        pwm.status, 
                                                        pwm.tf_family, 
                                                        pwm.motif_type,
                                                        pwm.msource_author,
                                                        pwm.msource,
                                                        pwm.pubmedid
                                                        ))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE pwm")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE pwm", error)
        lib.log.info(str(pwm.locus_tag) + " " +
                    str(pwm.motif_id) + " " +
                    str(pwm.status) + " " +
                    str(pwm.tf_family) + " " +
                    str(pwm.motif_type) + " " +
                    str(pwm.msource_author) + " " +
                    str(pwm.msource) + " " +
                    str(pwm.pubmedid))