import os
import psycopg2

def database_setup():

    conn = psycopg2.connect(
        database = os.environ.get('PG_DB'),
        user = os.environ.get('PG_USER'),
        password = os.environ.get('PG_PASS'),
        host = os.environ.get('PG_HOST'),
        port = os.environ.get('PG_PORT'),
    )

    print ("database_setup: opened database connection successfully")


