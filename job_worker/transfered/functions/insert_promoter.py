"""
    Insert promoter
"""
def insert_promoter(promoter):
    try:
        cursor = dbConnection.cursor()
        count = cursor.execute("INSERT INTO promoter VALUES (%s, %s, %s, %s, %s)",
                                                            (promoter.locus_tag, 
                                                            promoter.strand, 
                                                            promoter.source, 
                                                            promoter.start, 
                                                            promoter.stop))
        dbConnection.commit()
        lib.log.info("Record inserted successfully into TABLE promoter")
        cursor.close()
    except (Exception, psycopg2.Error) as error:
        lib.log.info("Failed to insert data into TABLE promoter", error)
        lib.log.info(str(promoter.locus_tag) + " " +
                    str(promoter.strand) + " " + 
                    str(promoter.source) + " " +
                    str(promoter.start) + " " +
                    str(promoter.stop))