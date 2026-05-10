"""
Extracts unique PfamUIDs and their MLLossRates from PfamMLLikeBayesPipeline_MLLoss_AnimalOnly.txt
file

Extract MeanAAComp_AnimalSpecific for each unique PfamUID from
    PFAMphylostratigraphy.PfamUIDsTable_EnsemblAndNCBI

Write data to a csv file, Pfam_Metrics.csv


order of amino acids in MySQL tables:
['A','R','N','D','C','E','Q','G','H','O','I','L','K','M','F','P','U','S','T','W','Y','V']
"""

import time
import mysql.connector
import sshtunnel
import pandas as pd
import os
from dotenv import load_dotenv
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def initiate_ssh_tunnel():
    load_dotenv()

    # SSHTunnel user-specs
    ssht_address = ('128.196.195.140', 22)
    ssht_username = 'garretwilson'                      # user specific
    ssht_password = os.getenv("SSHT_PASSWORD")          # user specific
    remote_bind_address = ("127.0.0.1", 3306)

    # Establish SSHTunnel
    server = sshtunnel.SSHTunnelForwarder(
        ssh_address_or_host=ssht_address,
        ssh_username=ssht_username,
        ssh_password=ssht_password,
        remote_bind_address=remote_bind_address
    )

    server.start()

    return server


def initiate_mysql_connection(server):
    # MySQL user-specs
    user = 'garretwilson'
    password = os.getenv("MYSQL_PASSWORD")              # user specific
    database = 'PFAMphylostratigraphy'                  # user specific

    # Establish MySQL connection
    cnx = mysql.connector.connect(
        user=user,
        password=password,
        database=database,
        port=server.local_bind_port,
        connect_timeout=1000,
        use_pure=True
    )

    print("\nMade Connection\n\n")

    return cnx


def get_pfamsUIDs_MLLossRates_AAfreq(in_file, data, cursor):
    '''
    Get all unique PfamUIDs and their MLLossRates from the in_file.
    PfamUIDs are the first column (index  0)
    MLLossRates are the seventh column (index 6)
    :param cursor: object to execute SQL queries
    :param data: dictionary to hold extracted PfamUIDs, MLLossRates, and amino acid frequency data
    :param in_file: PfamMLLikeBayesPipeline_MLLoss_AnimalOnly.txt
    :return: dictionary with PfamUIDs, MLLossRates, and amino acid frequencies
    '''
    file = open(in_file).readlines()
    for line in file[1:]:  # in_file has header
        temp_line = line.split("\t")
        pfam = temp_line[0]
        MLLossRate = float(temp_line[6])
        AAfreq = get_pfam_AAfreq(pfam, cursor)          # AAfreq is a list

        data["PfamUID"].append(pfam)
        data["MLLossRate"].append(MLLossRate)

        data_keys = list(data.keys())
        index = 0
        for k in data_keys[2:]:         # only the 22 amino acid keys
            # if k != "PfamUID" and k != "MLLossRate":
            freq = float(AAfreq[index])
            data[k].append(freq)
            index += 1

    return data


def get_pfam_AAfreq(pfam, cursor):
    '''
    Get MeanAAComp_AnimalSpecific for each unique PfamUID from
    PFAMphylostratigraphy.PfamUIDsTable_EnsemblAndNCBI
    :param pfam: unique PfamUID
    :param cursor: object to execute SQL queries
    :return: a list with the average frequency of each amino acid for a Pfam
    '''
    # table to access
    data_table = "PfamUIDsTable_EnsemblAndNCBI"

    # SQL query
    select_data_query = "SELECT MeanAAComp_AnimalSpecific FROM " + data_table + " WHERE PfamUID " \
                        "LIKE " + "'%" + pfam + "%'"

    cursor.execute(select_data_query)

    # format of .fetchall() output is [('AAcomp, AAcomp, ...'),] a list with one tuple with one
    # string
    raw_AAfreq_data = cursor.fetchall()

    # when no AA freq data in PfamUIDsTable_EnsemblAndNCBI table
    # in PfamUIDsTable_EnsemblAndNCBI table, empty cells are handled as below
    if raw_AAfreq_data == [('',)] or raw_AAfreq_data == [(None,)]:          # 12 Pfams
        print(pfam)
        print(raw_AAfreq_data)
        print("No amino acid frequency data in PfamUIDsTable_EnsemblAndNCBI. Will attempt to "
              "source from PFAMphylostratigraphy.NCBIGenomes_DomainMetrics_Complete \n"
              "and PFAMphylostratigraphy.EnsemblGenomes_DomainMetrics_Complete\n")

        return alternate_AAcomp_data(pfam, cursor)

    # format of pfam_AAfreq_data is [AAcomp, AAcomp, ...], a list
    else:
        pfam_AAfreq_data = raw_AAfreq_data[0][0].split(",")
        if len(pfam_AAfreq_data) != 22:
            print(pfam)
            print(pfam_AAfreq_data)
            print("Pfam frequency data does not have entries for all of the 22 amino acids\n")
        return pfam_AAfreq_data


def alternate_AAcomp_data(pfam, cursor):
    '''
    PfamUID had no mean amino acid frequency data in PfamUIDsTable_EnsemblAndNCBI. Will attempt to
    pull instances amino acid frequencies from
    PFAMphylostratigraphy.EnsemblGenomes_DomainMetrics_Complete and
    PFAMphylostratigraphy.NCBIGenomes_DomainMetrics_Complete and average the frequencies of the Pfam
    instances. If no data is found in either of these tables, we move to
    PFAMphylostratigraphy.PfamAlignments, extract the peptide sequences of instances and calculate
    average amino acid frequencies from peptide instances
    :param pfam: a PfamUID
    :param cursor: cursor object to execute queries
    :return: a list with the average frequency of each amino acid for a Pfam
    '''
    # tables to access
    data_table_Ensembl = "EnsemblGenomes_DomainMetrics_Complete"
    data_table_NCBI = "NCBIGenomes_DomainMetrics_Complete"

    # SQL query
    query_Ensembl = "SELECT PercentAminoAcidComposition FROM " + data_table_Ensembl + \
                    " WHERE PfamUID LIKE " + "'%" + pfam + "%'"

    cursor.execute(query_Ensembl)
    raw_Ensembl_data = cursor.fetchall()
    # print("raw Ensembl data:")
    # print(raw_Ensembl_data)

    # SQl query
    query_NCBI = "SELECT PercentAminoAcidComposition FROM " + data_table_NCBI + " WHERE PfamUID " \
                 "LIKE " + "'%" + pfam + "%'"

    cursor.execute(query_NCBI)
    raw_NCBI_data = cursor.fetchall()
    # print("raw NCBI data:")
    # print(raw_NCBI_data)

    # combine data from NCBI and Ensembl
    raw_Ensembl_data.extend(raw_NCBI_data)

    # format of .fetchall() output is [('AAcomp, AAcomp, ...'), ('...')], a list of tuples
    # representing instances and each tuple contains a single string of instance frequencies
    all_AAfreq_data = raw_Ensembl_data
    # print("all AA freq data:")
    # print(all_AAfreq_data)
    # print(len(all_AAfreq_data))
    # print("\n")

    # if no comp data in either table, move to peptides
    if all_AAfreq_data == []:  # empty cells in NCBI and Ensembl tables are []
        print("No amino acid frequency data in alternate tables. Sourcing from peptide data in "
              "PFAMphylostratigraphy.PfamAlignments\n\n")

        return AAfreq_from_peptide(pfam, cursor)

    # cleaned_AAfreq_data format is a 2D list with floats [[AAcomp, AAcomp, AAcomp, ...], [...]]
    # required to make data frame
    cleaned_AAfreq_data = []
    for tup in all_AAfreq_data:
        temp = tup[0].split(",")
        temp_float = [float(i) for i in temp]
        cleaned_AAfreq_data.append(temp_float)
    # print(cleaned_AAfreq_data)
    # print(len(cleaned_AAfreq_data))

    # create data frame
    df_AAfreq = pd.DataFrame(cleaned_AAfreq_data, columns=['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G',
                                                           'H', 'O', 'I', 'L', 'K', 'M', 'F', 'P',
                                                           'U', 'S', 'T', 'W', 'Y', 'V'])

    # calculate mean of each column (each AA)
    AAfreq_means = df_AAfreq.mean()

    # return a list with average freq of each AA
    # print(AAfreq_means.tolist())
    return AAfreq_means.tolist()


def AAfreq_from_peptide(pfam, cursor):
    '''
    Extract Pfam peptides from PFAMphylostratigraphy.PfamAlignments, calculate amino acid frequency
    of each instance, and average the frequencies of the instances.
    Utilize the BioPython class ProteinAnalysis to calculate amino acid frequencies from peptides
    :param pfam: a PfamUID
    :param cursor: cursor object to execute queries
    :return: a list with the average frequency of each amino acid for a Pfam
    '''
    # table to access
    data_table = "PfamAlignments"

    # SQL query
    query_peptides = "SELECT AlignedPeptide FROM " + data_table + " WHERE PfamUID " \
                     "LIKE " + "'%" + pfam + "%'"

    cursor.execute(query_peptides)

    # format of .fetchall() is [('peptide',), ('peptide',), ...]
    raw_peptides = cursor.fetchall()
    # print(raw_peptides)
    # print(len(raw_peptides))

    # order of amino acids required in final data frame
    amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'O', 'I', 'L', 'K', 'M', 'F', 'P',
                   'U', 'S', 'T', 'W', 'Y', 'V']

    # initiate list to hold lists of amino acid frequency data of each instance
    all_AAfreq_data = []
    for tup in raw_peptides:
        peptide_string = ''.join(aa for aa in tup[0] if aa.isalpha())
        # print(peptide_string)

        # BioPython class, ProteinAnalysis, used to calculate amino acid frequencies from peptide
        # string
        peptide_object = ProteinAnalysis(peptide_string)
        peptide_freq_dict = peptide_object.amino_acids_percent
        # print(peptide_freq_dict)

        peptide_freq_lst = []
        for aa in amino_acids:
            if aa not in peptide_freq_dict.keys():      # O and U
                peptide_freq_lst.append(0.0)
            else:
                peptide_freq_lst.append(peptide_freq_dict[aa])
        # print(peptide_freq_lst)


        all_AAfreq_data.append(peptide_freq_lst)

    # print(all_AAfreq_data)
    # print(len(all_AAfreq_data))


    # create data frame
    df_AAfreq = pd.DataFrame(all_AAfreq_data, columns=['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H',
                                                       'O', 'I', 'L', 'K', 'M', 'F', 'P', 'U', 'S',
                                                       'T', 'W', 'Y', 'V'])

    # calculate mean of each column (each AA)
    AAfreq_means = df_AAfreq.mean()

    # print(AAfreq_means.tolist())
    return AAfreq_means.tolist()


def main():
    start = time.time()

    # initiate ssh tunnel and mysql connection
    server = initiate_ssh_tunnel()
    cnx = initiate_mysql_connection(server)

    # make object to curse through mySQL tables
    cursor = cnx.cursor(buffered=True)

    # dictionary object to store extracted data and build data frame
    data = {
        "PfamUID": [],
        "MLLossRate": [],
        "A": [],
        "R": [],
        "N": [],
        "D": [],
        "C": [],
        "E": [],
        "Q": [],
        "G": [],
        "H": [],
        "O": [],
        "I": [],
        "L": [],
        "K": [],
        "M": [],
        "F": [],
        "P": [],
        "U": [],
        "S": [],
        "T": [],
        "W": [],
        "Y": [],
        "V": []}

    # text file to pass to function with PfamUID and MLLossRate
    in_file = "data/PfamMLLikeBayesPipeline_MLLoss_AnimalOnly.txt"

    # get all unique PfamUIDs and MLLossRates from in_file and extract amino acid frequency data
    # from MySQL
    pfamUIDs_MLLossRates_AAfreqs = get_pfamsUIDs_MLLossRates_AAfreq(in_file, data, cursor)

    # build data frame
    df = pd.DataFrame(pfamUIDs_MLLossRates_AAfreqs)
    df.set_index("PfamUID", inplace=True)
    # pd.set_option("display.max_columns", None)
    # pd.set_option("display.max_rows", None)
    print(df)

    # write to csv
    df.to_csv("data/Pfam_Metrics_2.csv")

    # track process time
    end = time.time()
    process_time = end - start
    print("\n\n" + "run time: " + str(process_time))


main()

# look out for these in peptide sequences:
# X, J, B, Z
