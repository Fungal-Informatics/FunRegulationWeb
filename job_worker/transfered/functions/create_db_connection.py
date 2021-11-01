
""" create a database connection to a SQLite database
"""
def create_db_connection():
    try:
        con = psycopg2.connect(host='localhost', database='funregulation',
        user='postgres', password='123456')
        lib.log.info("Successfully Connected to PostgreSQL")
        return con
    except (Exception, psycopg2.Error) as error:
        lib.log.info(error)
        


