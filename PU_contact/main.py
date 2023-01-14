import Bio.PDB
import numpy
import pandas
import os

#Create path for inputs folder
path_out = "OUT_PU"
path_pdb = "PDB"

# Exctract list of protein name from PDB folder
def list_pdb (path_pdb):
    pdb = []
    pdb_name = []
    for file in os.listdir(path_pdb):
        pdb.append(file)
        pdb_name.append(file.split('.')[0])
    return pdb_name

# Get all PU coordonates and names from Output file
def parserPU (file, path_out):
    with open(path_out + "/" + file,'r') as f:
        PU_list = []
        PU_name = []
        PU = []
        for line in f:
            l = line.split(' ')
            PU_list.append(l[4])
            PU_name.append(l[1])
        for i in PU_list:
            PU.append(i.split('-'))
    return PU, PU_name

# Extract residue for a PU list
def extract (chain, PU_coord, prot, path_pdb):
    structure = Bio.PDB.PDBParser().get_structure(prot, path_pdb + "/" + prot + ".pdb")
    model = structure[0]
    PU_list_res = []
    for index,value in enumerate(PU_coord):
        PU = []
        for residue in model[chain]:
           residu_id = residue.get_id()[1]
           if residu_id >= int(value[0]) and residu_id <= int(value[1]):
               PU.append(residue)
        PU_list_res.append(PU)
    return PU_list_res

# Calculate the CA distance between 2 residues
def calc_residue_dist(residue_one, residue_two) :
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    dist = numpy.sqrt(numpy.sum(diff_vector * diff_vector))
    return dist

# Calculate a matrix of CA distance between 2 PU
def calc_dist_matrix_PU(PU_1, PU_2) :
    matrix = numpy.zeros((len(PU_1), len(PU_2)), numpy.float)
    for row, residue_one in enumerate(PU_1) :
        for col, residue_two in enumerate(PU_2) :
            matrix[row, col] = calc_residue_dist(residue_one, residue_two)
    return matrix

# Calculate a matrix of number of contact between 2 PU for 2 lists of PU
def nb_contact_matrix (PU_list_1, PU_list_2) :
    matrix = numpy.zeros((len(PU_list_1), len(PU_list_2)), numpy.float)
    for row, i in enumerate(PU_list_1) :
        for col, j in enumerate(PU_list_2) :
            map = calc_dist_matrix_PU(i, j)
            contact_map = map < 8
            matrix[row, col] = numpy.count_nonzero(contact_map)
    return matrix

# Create a dataframe with contact between 2 lists of PU
def create_dataframe (nb_contact, PU_name_1, PU_name_2):
    df = pandas.DataFrame(data = nb_contact, index = PU_name_1, columns = PU_name_2)
    df = df.groupby(df.index).sum()
    df = df.groupby(df.columns, axis =1).sum()
    return df

# Concatenate 2 dataframe
def concatenate_sum_df (df_1, df_2):
    df_tot = pandas.concat([df_1,df_2],axis=1)
    df_tot = df_tot.groupby(df_tot.index).sum()
    df_tot = df_tot.groupby(df_tot.columns, axis = 1).sum()
    return df_tot

# Create Dataframe with all contacts beetween PU from the folder for proteins in the folder
def create_dataframe_out (path_out, path_pdb):
    list_prot = list_pdb(path_pdb)
    DF_prot = []
    for prot in list_prot:
        PU_chainA, PU_name_A = parserPU("Out_" + prot + "_A", path_out)
        PU_chainB, PU_name_B = parserPU("Out_" + prot + "_B", path_out)
        PU_list_res_A = extract ("A", PU_chainA, prot, path_pdb)
        PU_list_res_B = extract ("B", PU_chainB, prot, path_pdb)
        print("\n" + "protein = " + prot) # print the name of the protein to know where is the program for the calculation of the contact matrix
        print("Calcul contact AA") # print the 2 chain where the contact is calculate to know where is the program for the calculation of the contact matrix
        nb_contact_AA = nb_contact_matrix(PU_list_res_A, PU_list_res_A)
        print("Calcul contact AB") # print the 2 chain where the contact is calculate to know where is the program for the calculation of the contact matrix
        nb_contact_AB = nb_contact_matrix(PU_list_res_A, PU_list_res_B)
        print("Calcul contact BB") # print the 2 chain where the contact is calculate to know where is the program for the calculation of the contact matrix
        nb_contact_BB = nb_contact_matrix(PU_list_res_B, PU_list_res_B)
        df_A = create_dataframe (nb_contact_AA, PU_name_A, PU_name_A)
        df_B = create_dataframe (nb_contact_BB, PU_name_B, PU_name_B)
        df_AB = create_dataframe (nb_contact_AB, PU_name_A, PU_name_B)
        df_totA_B = concatenate_sum_df(df_A,df_B)
        df_tot = concatenate_sum_df(df_totA_B,df_AB)
        DF_prot.append(df_tot)
    for DF in DF_prot:
        for DFj in DF_prot:
            DF = pandas.concat([DF, DFj],axis=1)
            DF = DF.groupby(DF.index).sum()
            DF = DF.groupby(DF.columns, axis = 1).sum()
    return DF

DF = create_dataframe_out (path_out, path_pdb)
print(DF) # print the final dataframe
DF.to_csv("dataframe_PU.csv") # Save the output dataframe in csv
