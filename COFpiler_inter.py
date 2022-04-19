# -*- coding: utf-8 -*-
# @Author  : Yingying Zhang
# @Time    : 4/19/2022 2:39 PM
# @Function: Build the statistial model with intercalated monomer
# @Usage: python COFpiler_inter.py -data data_file -i input_structure -o out[ut_format -path /file_path/ -L layer_number
# -M model_number --set set_number
# input structure of intercalated monomer, structures above and below the monomer should be saved separately, with name
# of stakcing_type1/2/3

from ase.io import read, write
import numpy as np
from ase.visualize import view
import math
import pandas as pd
import argparse


# calculate the probabilities from the relative energies
def get_probability(energy_data, temperature):
    if energy_data is None:
        p_out = energy_data
    else:
        energies_in = np.array(energy_data)
        p_exp_sum = 0
        p_exp = np.zeros(len(energies_in))
        p_out = np.zeros(len(energies_in))
        for i in range(len(energies_in)):
            p_exp[i] = np.exp(-energies_in[i] * 120.37 / temperature)
            p_exp_sum += p_exp[i]
        for n in range(len(energies_in)):
            p_out[n] = p_exp[n] / p_exp_sum
    p_sum = np.sum(p_out)
    if p_sum != 1:
        p_out[len(energies_in) - 1] += 1 - p_sum
    return p_out


# redefine the lattice parameter
def def_cell(structure, s_vectorSum, s):
    c_n = s_vectorSum + s
    l_a = structure.get_cell()[0]
    l_b = structure.get_cell()[1]

    while np.abs(c_n[1]) > 0.8 * np.abs(l_b[1]):
        if (c_n[1] > 0 and l_b[1] > 0) or (c_n[1] < 0 and l_b[1] < 0):
            c_n -= l_b
        elif (c_n[1] < 0 and l_b[1] > 0) or (c_n[1] > 0 and l_b[1] < 0):
            c_n += l_b
    while np.abs(c_n[0]) > 0.8 * np.abs(l_a[0]):
        if (c_n[0] > 0 and l_a[0] > 0) or (c_n[0] < 0 and l_a[0] < 0):
            c_n -= l_a
        elif (c_n[0] < 0 and l_a[0] > 0) or (c_n[0] > 0 and l_a[0] < 0):
            c_n += l_a
    if c_n[1] < 0 and np.arctan(np.abs(c_n[2] / c_n[1])) < np.radians(50):
        c_n += l_b
    if c_n[0] < 0 and np.arctan(np.abs(c_n[2] / c_n[0])) < np.radians(50):
        c_n += l_a
    return c_n


# rotate the structure, structure has iso-energetic shift vectors because of the symmetry
def rot(s, num_rot, rot_angle):
    rots = np.zeros([num_rot])
    for j in range(num_rot):
        rots[j] = 1 / num_rot
    rot_p = np.random.choice(num_rot, p=rots)
    rot_angle += rot_p * 2* np.pi / num_rot
    rot_matrix = ([math.cos(rot_angle), math.sin(rot_angle)], [-math.sin(rot_angle), math.cos(rot_angle)])
    s_vector_rot = np.dot(rot_matrix, ([s[0], s[1]]))  # rotate clockwise
    s_vector[0:2] = [s_vector_rot[0], s_vector_rot[1]]
    return s_vector, rot_angle


# get the mirror structure and slip direction, atom1 and atom2 are the diagonal atoms of the COF layer
def mirror(str_in, s_vector):
    center = str_in.get_positions().mean(axis=0)
    s_vector_x = s_vector[1]
    s_vector_y = s_vector[0]
    s_vector[0:2] = [s_vector_x, s_vector_y]
    for atom in str_in:
        x = atom.position[1] - center[1] + center[0]
        y = atom.position[0] - center[0] + center[1]
        atom.position[0:2] = [x,y]
    return str_in, s_vector


# ###################### Input parameters ######################
parser = argparse.ArgumentParser(description="Create a structure of a 2D layer which follows a parametric function")
parser.add_argument("-data", type=str, help="The input file include data needed, form as: form as: the 1st column "
                                             "is the stacking_type, the 2nd is the Erel (relative energies), the 3rd, "
                                             "4th and 5th columns are the x, y, z (vector of the shift_vectors)")
parser.add_argument("-i", "--instr", type=str, help="The input structure")
parser.add_argument("-T", "--tem", type=int, default=293, help="The experimental synthesize temperature")
parser.add_argument("-path", type=str, help="The path where the infile and instr are, the output structure will also "
                                             "save there")
parser.add_argument("-o", "--output_format", type=str, default="cif", help="The output structure format")
parser.add_argument("-m", "--mirror", type=bool, default=True, help="enable the consider of mirror shift or not")
parser.add_argument("-mp", "--mplane", type=list, default=[0.001, 1], help="the mirror plane for s2 shift")
parser.add_argument("-L", type=int, help="intercalated monomer every L layer")
parser.add_argument("-M", type=int, help="the model number")
parser.add_argument("--set", type=int, help="the number of sets")

args = parser.parse_args()
n_rot = 4

# ####################### READING DATA #######################
print('--------------------------------------')
print('Reading data...')
data = pd.read_excel(args.path + args.data, sheet_name='data')

stacking_types = data['stacking_type']
data = data.set_index('stacking_type')
energies = data['Erel']
s_vectors = data[['x', 'y', 'z']]

prob = get_probability(energies, args.tem)  # calculate probabilities from given energies
data['prob'] = prob
print(data)

f = open(args.path + "record_" + str(args.L) + '.txt', 'w+')  # record the stacking type for the randomize structure
print('stacking_types are:\n'+str(stacking_types)+'\n', file=f)
print('\nprobabilities are:\n' + str(prob)+'\n', file=f)
print('\ns_vectors are:\n' + str(s_vectors)+'\n', file=f)

for fnum in range(1, args.M+1):
    print('--------------------------------------')
    print('Building the %sst statistical model...' %fnum)
    print('\nfile ' + str(fnum), file=f)
    s_p = np.random.choice(len(prob), p=prob)  # AA_e,AA_n,AA_l,AB_e,AB_n,AB_l
    stacking_type = stacking_types[s_p]
    s_vector_sum = [0, 0, 0]
    angle_sum = 0
    m = 0

    file1 = read(args.path + stacking_type + '1.xyz')
    file3 = read(args.path + stacking_type + '3.xyz')
    center_base = (file1.get_positions().mean(axis=0) + file3.get_positions().mean(axis=0))/2
    structure = file1.copy()
    lattice = file1.get_cell()
    s_vector = s_vectors.loc[stacking_type].copy()

    for layer in range(2, args.set*(args.L+1)+1):    # generate the shift vector for the 2nd layer
        if (layer+1) % (args.L + 1) == 0:   # the intercalated monomer
            s_p = np.random.choice(len(prob), p=prob)
            stacking_type = stacking_types[s_p]
            s_vector = s_vectors.loc[stacking_type].copy()
            str_in = read(args.path + stacking_type + '2.xyz')
            center = str_in.get_positions().mean(axis=0)
            # move the intercalated monomer back to the center, which will have a new rotated slip vector
            for i in range(len(str_in)):
                str_in[i].position[0:2] = str_in[i].position[0:2] - s_vector[0:2]
            mirror_p = np.random.choice(2, p=(1/2, 1/2))
            if mirror_p == 0:
                str_in = mirror(str_in.copy(), s_vector)[0]
                s_vector = mirror(str_in.copy(), s_vector)[1]
        elif (layer) % (args.L + 1) == 0:  # the layers above
            str_in = read(args.path + stacking_type + '3.xyz')
        else:  # the layers below
            str_in = read(args.path + stacking_type + '1.xyz')
            if (layer+args.L) % (args.L + 1) == 0:
                pass
            else:
                m+=1

        # the input structures are from L=2, so for more L, the z position of the initial structure should be changed.
        for atom in str_in:
            atom.position[2] += 1/3 * lattice[2,2] * m

        if (layer+1) % (args.L + 1) == 0:   # the intercalated monomer
            s_vector = rot(s_vector, n_rot, angle_sum)[0]
            angle_sum = rot(s_vector, n_rot, angle_sum)[1]

            for i in range(len(str_in)):
                structure.append(str_in.get_chemical_symbols()[i])
                structure[-1].position[0:2] = str_in[i].position[0:2] + s_vector[0:2]
                structure[-1].position[2] = str_in[i].position[2] + (s_vector[2]) * math.floor((layer - 1) / (args.L + 1))
            print(layer, '\t', stacking_type, '\t', s_vector.values, '\t', s_vector_sum, file=f)
        else:
            for i in range(len(str_in)):
                structure.append(str_in.get_chemical_symbols()[i])
                structure[-1].position[0:2] = str_in[i].position[0:2]
                structure[-1].position[2] = str_in[i].position[2] + (s_vector[2]) * math.floor((layer - 1) / (args.L + 1))

        # to count the final lattice c
        if (layer) % (args.L + 1) == 0:
            s_vector_sum += s_vector.values + 1/3 * (args.L-2) * s_vector.values

    c = def_cell(structure, s_vector_sum, [s_vector.values[0], s_vector.values[1], 0])
    structure.set_cell([structure.get_cell()[0], structure.get_cell()[1], c], scale_atoms=False)
    view(structure)

    out_filename = str(stacking_types[s_p]) + '_' + str(fnum) + '_' + str(args.L) + 'l.' + args.output_format
    write(args.path + out_filename, structure, format=args.output_format)

f.close()
