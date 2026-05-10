'''
extracts some data on PFAM domain length from mysql database
'''

import time
import mysql.connector
import sshtunnel
import pandas as pd
from dotenv import load_dotenv
import os


def initiate_ssh_tunnel():
    load_dotenv()

    # SSHTunnel user-specs
    ssht_address = ('128.196.195.140', 22)
    ssht_username = 'garretwilson'                      # user specific
    ssht_password = os.getenv("SSHT_PASSWORD")          # user specific
    remote_bind_addy = ("127.0.0.1", 3306)

    # Establish SSHTunnel
    server = sshtunnel.SSHTunnelForwarder(
        ssh_address_or_host=ssht_address,
        ssh_username=ssht_username,
        ssh_password=ssht_password,
        remote_bind_address=remote_bind_addy
    )
    server.start()

    return server


def initiate_mysql_connection(server):
    # MySQL user-specs
    User = 'garretwilson'
    Password = os.getenv("MYSQL_PASSWORD")              # user specific
    Database = 'PFAMphylostratigraphy'                  # user specific

    # Establish MySQL connection
    cnx = mysql.connector.connect(
        user=User,
        password=Password,
        database=Database,
        port=server.local_bind_port,
        connect_timeout=1000,
        use_pure=True
    )
    print("\nMade Connection\n\n")

    return cnx


def get_PfamUID_DomainLength_Ensembl(cursor, Ensembl_dict):
    query = "SELECT PfamUID, DomainLength FROM EnsemblGenomes_DomainMetrics_Complete"

    cursor.execute(query)

    # data is of format [(PfamUID, DomainLength), ...]
    data = cursor.fetchall()

    # put fetched data into data_Ensembl
    for i in data:
        Ensembl_dict["PfamUID"].append(i[0])
        Ensembl_dict["DomainLength"].append(i[1])

    print("\nEnsembl Domains Complete\n")
    return Ensembl_dict


def get_PfamUID_DomainLength_NCBI(cursor, NCBI_dict):
    query = "SELECT PfamUID, DomainLength FROM NCBIGenomes_DomainMetrics_Complete"

    cursor.execute(query)

    # data is of format [(PfamUID, DomainLength), ...]
    data = cursor.fetchall()

    # put fetched data into data_NCBI
    for i in data:
        NCBI_dict["PfamUID"].append(i[0])
        NCBI_dict["DomainLength"].append(i[1])

    print("\nNCBI Domains Complete\n")
    return NCBI_dict


def main():
    start = time.time()

    #initaite ssh tunnel and mysql connection
    server = initiate_ssh_tunnel()
    cnx = initiate_mysql_connection(server)

    # make object to curse through mySQL tables
    cursor = cnx.cursor(buffered=True)

    # data frame for Ensembl
    Ensembl_dict = {"PfamUID": [],
                    "DomainLength": []
                    }

    # data frame for NCBI
    NCBI_dict = {"PfamUID": [],
                 "DomainLength": []
                 }

    # get Ensembl PfamUIDs and Domain Lengths
    Ensembl_data = get_PfamUID_DomainLength_Ensembl(cursor, Ensembl_dict)

    # get NCBI PfamUIDs and Domain Lengths
    NCBI_data = get_PfamUID_DomainLength_NCBI(cursor, NCBI_dict)

    # build data frames and write csv
    Ensembl_df = pd.DataFrame(Ensembl_data)
    Ensembl_df.set_index("PfamUID", inplace=True)
    Ensembl_df.to_csv("data/DomainLength_Ensembl_PfamUID.csv")

    NCBI_df = pd.DataFrame(NCBI_data)
    NCBI_df.set_index("PfamUID", inplace=True)
    NCBI_df.to_csv("data/DomainLength_NCBI_PfamUID.csv")

    # end program
    cnx.commit()
    cnx.close()
    server.stop()
    end = time.time()
    print(" ")
    print("Process time: " + str((end-start)) + " s")

main()
