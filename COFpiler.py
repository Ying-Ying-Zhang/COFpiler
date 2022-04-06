from ase.io import read, write
import numpy as np
from ase.visualize import view
import math
import pandas as pd
from scipy.optimize import fsolve
import os
import re
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
    if 'AB_' in stacking_type:
        c_n = s_vectorSum - s
        c_n[2] = s_vectorSum[2] + s[2]
    elif 'AA' or 'ABC' in stacking_type:
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


# obtain the mirrored shift vector
def mirror(i):
    k = (mirror_plane[1] / mirror_plane[0])
    x0, y0 = s_vector[0], s_vector[1]
    x, y = i[0], i[1]
    return [
        y0 + y - k * (x0 + x),
        k * (y - y0) + x - x0
    ]


# rotate the structure, structure has iso-energetic shift vectors because of the symmetry
def rot(s, num_rot, rot_angle):
    rots = np.zeros([num_rot])
    for j in range(num_rot):
        rots[j] = 1 / num_rot
    rot_p = np.random.choice(num_rot, p=rots)
    rot_angle += rot_p * np.pi / 3
    rot_matrix = ([math.cos(rot_angle), math.sin(rot_angle)], [-math.sin(rot_angle), math.cos(rot_angle)])
    s_vector_rot = np.dot(rot_matrix, ([s[0], s[1]]))  # rotate clockwise
    s_vector[0:2] = [s_vector_rot[0], s_vector_rot[1]]
    return s_vector


# ###################### Input parameters ######################
parser = argparse.ArgumentParser(description="Create a structure of a 2D layer which follows a parametric function")
parser.add_argument("-data", type=str, help="The input file include data needed, form as: form as: the 1st column "
                                             "is the stacking_type, the 2nd is the Erel (relative energies), the 3rd, "
                                             "4th and 5th columns are the x, y, z (vector of the shift_vectors)")
parser.add_argument("-i", "--instr", type=str, help="The input structure")
parser.add_argument("-T", "--tem", type=int, default=293, help="The experimental synthesize temperature")
parser.add_argument("-path", type=str, help="The path where the infile and instr are, the output structure will also "
                                             "save there")
parser.add_argument("-o", "--outstr_format",  type=str, default="cif", help="The output structure format")
parser.add_argument("-s", "--symmetry", type=str, default='C6', help="the symmetry for the structure, C3, C6 or C4")
parser.add_argument("-m", "--mirror", type=bool, default=True, help="enable the consider of mirror shift or not")
parser.add_argument("-mp", "--mplane", type=list, default=[0.001, 1], help="the mirror plane for s2 shift")
parser.add_argument("-L", type=int, help="the layer number")
parser.add_argument("-M", type=int, help="the model number")
output_folder = '/statistical_structures/'

args = parser.parse_args()

L = args.L  # the number of layer
M = args.M  # the number of model

# ####################### READING DATA #######################
print('READING DATA')
file_mono = read(args.path + args.instr)
atom_num = file_mono.get_number_of_atoms()

lattice = file_mono.get_cell()
data = pd.read_excel(args.path + args.data, sheet_name='data')
stacking_types = data['stacking_type']
data = data.set_index('stacking_type')

# the shift vector divides to two part, s_lattice and s_vector. Therefore, they can rotate separately and combined after
# rotation. The s_lattice for symmetry of 4 might differs, could be 1/2, or 1/3, it needs to be changed by user.
n_rot = int(re.sub("\D", "", args.symmetry))

if n_rot == 4:
    AB_s_type = input("which AB slip is included? AB_s_axis or AB_s_diag: ")
    if AB_s_type == 'AB_s_axis':
        data.loc['AB_s2_diag', 'Erel'] = 0
        s_AB_lattice = np.array([0, 1 / 2 * lattice[0, 0], 0])
    elif AB_s_type == 'AB_s_diag':
        data.loc['AB_s2_axis', 'Erel'] = 0
        s_AB_lattice = np.array([1 / 2 * lattice[0, 0], 1 / 2 * lattice[0, 0], 0])
else:
    s_AB_lattice = np.array([0, np.sqrt(3) / 3 * lattice[0, 0], 0])

energies = data['Erel']
prob = get_probability(energies, args.tem)  # calculate probabilities from given energies
s_vectors = data[['x', 'y', 'z']]
data['prob'] = prob

print(data)
# ####################### START ################################
angle_prob = []
for angle in range(n_rot):
    angle_prob.append(1 / n_rot)

if not os.path.exists(args.path + output_folder):
    os.mkdir(args.path + output_folder)
out_path = (args.path + output_folder)

# the shift for AA and AB/ABC respectively
s_AA_lattice = np.array([0, 0, 0])

f = open(out_path + "record_" + str(L) + '.txt', 'w+')  # record the stacking type for the randomize structure
f.write('stacking_types are:\n' + str(stacking_types) + '\n')
f.write('\np_2nd are:\n' + str(prob) + '\n')
f.write('\ns_vectors are:\n' + str(s_vectors) + '\n')

for fnum in range(M):
    print('Build the %dth model' % M)
    print('\nfile ' + str(fnum), file=f)
    print('layer_number', '\t', 'stacking_type', '\t', 's_vector', '\t\t\t\t', 's_vector_sum', file=f)
    print('layer_number', '\t', 'stacking_type', '\t', 's_vector', '\t\t\t\t', 's_vector_sum')

    str_out = file_mono.copy()
    angle = 0
    angle_lattice = 0
    s_vector_sum = [0, 0, 0]
    str_in = file_mono.copy()
    m = 1  # count the layer number is odd or even

    for layer in range(L - 1):
        s_p = np.random.choice(len(prob), p=prob)
        stacking_type = stacking_types[s_p]

        s2_p = np.random.choice(n_rot, p=angle_prob)
        angle += (s2_p + 1) * 2 * np.pi / n_rot

        s3_p = np.random.choice(n_rot, p=angle_prob)
        angle_lattice += (s3_p + 1) * 2 * np.pi / n_rot

        if 'AA' in stacking_type:
            s_lattice = s_AA_lattice.copy()
            s_vector = s_vectors.loc[stacking_type]
        elif 'AB_' in stacking_type:
            s_vector = s_vectors.loc[stacking_type].copy() - s_AB_lattice.copy()
            if (m % 2) == 0:
                s_lattice = -s_AB_lattice.copy()[0:2]
            else:
                s_lattice = s_AB_lattice.copy()[0:2]
            s_lattice = rot(s_lattice.copy().copy(), n_rot, angle_lattice)
        else:
            print('The stacking is neither AA, AB nor ABC')

        if 'ecl' in stacking_type:
            pass
        elif 's1' in stacking_type:
            s_vector = rot(s_vector, n_rot, angle)
        elif 's2' in stacking_type:
            if args.mirror:
                mirror_plane = args.mplane
                s_vector_new = fsolve(mirror, [0, 0])
                s_vector[0:2] = [s_vector_new[0], s_vector_new[1]]
            else:
                pass
            s_vector = rot(s_vector, n_rot, angle)
        else:
            print('the stacking type can not be recognized')

        if 'AB_' in stacking_type or 'ABC' in stacking_type:
            s_vector[0:2] += s_AB_lattice[0:2]
        else:
            pass
        m += 1

        s_vector_sum += s_vector.copy()
        for atom in range(atom_num):
            str_out.append(str_in.get_chemical_symbols()[atom])
            str_out.positions[-1] = str_in[atom].position + s_vector_sum
        print(layer + 2, '\t', stacking_type, '\t', s_vector.values, '\t', s_vector_sum.values, file=f)
        print(layer + 2, '\t\t', stacking_type, '\t', s_vector.values, '\t', s_vector_sum.values)

    # redefine the c lattice
    c = def_cell(str_out, s_vector_sum.values, s_vector.values)
    lattice_new = [lattice[0], lattice[1], c]
    str_out.set_cell(lattice_new, scale_atoms=False)
    out_filename = str(L) + '_' + str(fnum) + '.' + args.outstr_format
    write(out_path + out_filename, str_out, format=args.outstr_format)
f.close()
view(str_out)
