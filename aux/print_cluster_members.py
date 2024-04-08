#!/bin/python3
import sys; sys.path.append("/hits/basement/mcm/niessnjb/sdas/sda-master/auxi")
import os
from fort55_lib_SDA7 import *
import numpy as np
import pdb

def print_usage():
    print("This script prints all members of  the cluster into a table")
    print("")
    print("The columns of the table include:")
    print("No: Line number of the member in the complexes file")
    print("Energy: The total energy of the cluster")
    print("Occur: The number of occurences of the complex")
    print("El-, HD-, ED-, LjEnergy: The Electrostatoc, Hydrophobic des, Electrophobic des and Lennard-Jones energies repectively")
    print("rmsd: The RMSD from the representative of the the cluster")
    print("rmsd_original: The RMSD from the original structure")

    print("")
    print("")
    print("")
    print("Usage:")
    print("thisprogram complexes_file pdb_file stats_file treecutting_file")
    print("complexes_file: The original complexes file from the simulations")
    print("pdb_file: The pdb file of protein 2")
    print("stats_file: The file containing the redirected output of the clust program")
    print("treecutting_file: The Treecutting.out file from the clust program")


def apply_roto_translation(coords, trans, rot):
    return (rot@coords.T).T + trans

def compute_rmsd(coords1, coords2):
    diff_squared = (coords1-coords2)**2
    return sqrt(1/diff_squared.shape[0] * np.sum(diff_squared))



def get_atom_coordinates_from_file(pdbfile):
 
    str_bf, coord, str_af = ReadPDBFile(pdbfile)
    atom_names = list(map(lambda s: s[12:15].strip(), str_bf))

    #sorting out only the backbone atoms
    coordinates_list = [coordinate for idx, coordinate in enumerate(coord) if(bool(set('CNO') & set(atom_names[idx])))]

    coordinates = np.array(coordinates_list)

    return coordinates

def get_rotation_matrices_from_data(lines):
    
    #rotation vectors
    rot_x = lines[:, 5:8]
    rot_y = lines[:, 8:11]
    rot_z = np.cross(rot_x, rot_y)

    #Merging vectors to rotation matrices
    rotation_matrices = np.transpose([rot_x, rot_y, rot_z], axes=(1,2,0))
    return rotation_matrices





def main(complexes_file, pdb_file, stats_file, treecutting_file):
    try:
        Fort55file = open(complexes_file, 'r')
    except IndexError:
        print_usage()
        return 

    if not os.path.isfile(pdb_file):
        print("PDB file does not exist")
        print_usage()
        return 


    coordinates = get_atom_coordinates_from_file(pdb_file)

    line1 = Fort55file.readline()
    line2 = Fort55file.readline()
    line3 = Fort55file.readline()
    line4 = Fort55file.readline()

    fields=line1.split()[1:]

    xc1 = np.array(line3.split()[1:4]).astype(float)
    xc2 = np.array(line4.split()[1:4]).astype(float)

    #translation to center of mass of p2
    xv2 = coordinates-xc2

    
    #Getting treecutting information
    cuts = np.loadtxt(treecutting_file, dtype=int, usecols=4)
    number_clusters = np.max(cuts)+1

    lines_data = np.loadtxt(Fort55file)[:len(cuts)]
    
    #extracting the position values
    trans_vec = lines_data[:,2:5]

    #getting rotation matrices
    rotation_matrices = get_rotation_matrices_from_data(lines_data)

    #Getting the energy values of the complexes
    energies = {}
    energies["TotEn"]=np.zeros(lines_data.shape[0])
    energies["El"]=np.zeros(lines_data.shape[0])
    energies["ED"]=np.zeros(lines_data.shape[0])
    energies["HD"]=np.zeros(lines_data.shape[0])
    energies["rLj"]=np.zeros(lines_data.shape[0])
    energies["Occur"] = np.zeros(lines_data.shape[0], dtype=int)

    for idx, field in enumerate(fields):
        if (field.startswith("TotEn")):
            energies["TotEn"] += lines_data[:, idx]
        elif (field.startswith("El")):
            energies["El"] += lines_data[:, idx]
        elif (field.startswith("ED")):
            energies["ED"] += lines_data[:, idx]
        elif (field.startswith("HD")):
            energies["HD"] += lines_data[:, idx]
        elif (field.startswith("rLJ")):
            energies["rLJ"] += lines_data[:,idx]
        elif (field.startswith("Occur")):
            energies["Occur"] += lines_data[:, idx].astype(int)



    #Reading the representatives from the stats.out file
    stats = open(stats_file).readlines()[8]
    representatives = np.array(stats.split()).astype(int)-1

    #Reading the dictionary representative cluster num
    lines_stats = open(stats_file).readlines()[13:]
    
    i=0
    rep_clustno = {}
    while(lines_stats[i] != "\n"):
        fields=lines_stats[i].split()
        rep_clustno[int(fields[3])]=int(fields[0])
        print(f"{fields[3]} {fields[0]}")
        i+=1



    #analyzing
    for idx, rep in enumerate(representatives):
        
        #Maybe put this into a function
        rot_matrices_in_clust = rotation_matrices[cuts==idx]
        trans_vec_in_clust = trans_vec[cuts==idx]
        
        coordinates_in_clust_rotated = (rot_matrices_in_clust@(xv2.T)).transpose((2,0,1))
        coordinates_in_clust = (coordinates_in_clust_rotated+trans_vec_in_clust).transpose(1,0,2)


        

        rotmatrix_rep = rotation_matrices[rep]
        trans_vec_rep = trans_vec[rep]

        coordinates_rep = (rotmatrix_rep@(xv2.T)).T+trans_vec_rep


        rmsds = np.array([compute_rmsd(coordinates_rep, coordinates_cur) for coordinates_cur in coordinates_in_clust])
        rmsds_original = np.array([compute_rmsd(coordinates, coordinates_cur+xc1) for coordinates_cur in coordinates_in_clust])
        Energies = energies["TotEn"][cuts==idx]
        Occurences = energies["Occur"][cuts==idx]
        ElEnergy = energies["El"][cuts==idx]
        HDEnergy = energies["ED"][cuts==idx]
        EDEnergy = energies["HD"][cuts==idx]
        LjEnergy = energies["rLj"][cuts==idx]
        numbers = np.argwhere(cuts==idx)[:,0]

       # breakpoint()

        array = np.array([numbers, Energies, Occurences, ElEnergy, HDEnergy, EDEnergy, LjEnergy, rmsds, rmsds_original]).T

       # breakpoint()

        np.savetxt(f"members_clusters{rep_clustno[rep+1]}.out", array, header = "   No        Energy        Occur      ElEnergy     HDEnergy      EDEnergy      LjEnergy       rmsd   rsmd_original",
               fmt = " %5i   %6.4e   %6.4e   %6.4e   %6.4e   %6.4e   %6.4e   %6.4e  %6.4e", comments="#")





if __name__=='__main__':
    if (len(sys.argv)<5 or sys.argv[1]=="--h" or sys.argv[1]=="-h"):
        print_usage()
        sys.exit()
    
    complexes_file = sys.argv[1]
    pdb_file = sys.argv[2]
    stats_file = sys.argv[3]
    treecutting_file = sys.argv[4]
    main(complexes_file, pdb_file, stats_file, treecutting_file)
