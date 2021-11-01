
"""
    Select pwms by tf_locus_tag
"""
def select_pwms_by_locus_tag(locus_tag):
    pwm_list = list()
    try:
        cursor = dbConnection.cursor()
        postgreSQL_select_Query = "SELECT * from pwm WHERE locus_tag = %s"
        cursor.execute(postgreSQL_select_Query, (locus_tag,))
        rec = cursor.fetchall()
        for row in rec:
            pwm_list.append(Pwm(row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8]))
        cursor.close()
        return pwm_list
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to execute the select into table pwm", error)
        lib.log.info(locus_tag)