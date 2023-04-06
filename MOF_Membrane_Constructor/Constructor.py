import tkinter   #Version tcl/tk 8.6
from threading import Thread   #
import os
import re
import math
import numpy as np    #Version 1.23.0
from shapely import geometry  #Version 1.8.2
import time
from tkinter.filedialog import askdirectory


########settings#############
global slab_h   #m
global vacuum
file_dir = '.\\Example\\'#Input path
global outpath  # output path

facet_mode= "AUTO" # "AUTO" or "SET" to Control the manner to identify an exposed facet.
facet_index_input = "0 0 1"  #input facet index for the Set option to identify an exposed facet

grid_size = 0.8  #The input value of cell size to perform automatic identification of the exposed facet
output_files_type = ["cif"] #the list of the formats of output models
H_charge_set = 0.15   #if atomic charge is detected in the input, the surface saturated H will be charged with this value
global facet

separate_layer=0
protonation_yn = 1

min_X_set=0
min_Y_set=0

min_X_length=24
min_Y_length=24
########globale values#############
file_num = 0  # The number of input files of MOFs
log_lines=[]
error_lines = ["error_structure\n"]
#log_lines=["MOF-Membrane Constructor-Job Log\n\n"]

running = True

charge_all_N=[]
charge_all=[]
nonmetal_elemt_u=[]
metal_elemt_u=[]
metal_connect_structure_u=[]
global N2
f_con0=[]

Filename=''
facet_index="0 0 1"

# # the following are the radii of the atoms to determine their connectivity
atom_r = {"H": "0.38",  "Li": "0.86", "Be": "0.53", "B": "1.01", "C": "0.88", "N": "0.86", "O": "0.89",  "F": "0.82",
          "Na": "1.66", "Mg": "1.41", "Al": "1.53", "Si": "1.38","P": "1.28", "S": "1.20", "Cl": "1.17", "K": "2.03",
          "Ca": "1.76", "Sc": "1.70", "Ti": "1.65", "V": "1.51", "Cr": "1.53","Mn": "1.53","Fe": "1.43", "Co": "1.31",
          "Ni": "1.33", "Cu": "1.31", "Zn": "1.41", "Ga": "1.40","Ge": "1.35","As": "1.39","Se": "1.40", "Br": "1.39",
          "Rb": "2.2",  "Sr": "1.95", "Y": "1.84",  "Zr": "1.73","Nb": "1.66","Mo": "1.57","Ru": "1.58", "Rh": "1.63",
          "Pd": "1.68", "Ag": "1.56", "Cd": "1.56", "In": "1.53","Sn": "1.64","Sb": "1.64","Te": "1.65", "I": "1.58",
          "Cs": "1.85", "Ba": "1.52", "La": "1.91", "Ce": "1.98","Pr": "1.75","Nd": "1.92","Sm": "1.89", "Eu": "1.83",
          "Gd": "1.79", "Tb": "1.82", "Dy": "1.79", "Ho": "1.63","Er": "1.80", "Tm": "1.84","Yb": "1.80","Lu": "1.86",
          "Hf": "1.73", "W": "1.33",  "Re": "1.29", "Ir": "1.50","Pt": "1.66", "Au": "1.68","Hg": "1.88","Pb": "2.2",
          "Bi": "1.72", "Th": "1.97", "U": "1.76",  "Np": "1.73","Pu": "1.71", "Am": "1.80","Cm": "1.69",    }

nonmetal_list = ["H", "C", "B", "N", "O", "F", "Si", "P", "S", "Se", "Cl", "Br", "I","As","Te","Se"] #the nonmetalic elemnents
job_type=1
charge_r=1

def coa_in_grid(coa_input,a,b,c):
    rcoa_f_input=coa_conv2_rcoa(coa_input,a,b,c)
    R0_f,R1_f,R2_f =0,0,0

    if rcoa_f_input[0] >= 0 and rcoa_f_input[0] <= 1:
        R0_f = rcoa_f_input[0]
    elif rcoa_f_input[0] < 0:
        R0_f = round(rcoa_f_input[0] - int(rcoa_f_input[0]) + 1, 6)
    elif rcoa_f_input[0] > 1:
        R0_f = round(rcoa_f_input[0] - int(rcoa_f_input[0]), 6)

    if rcoa_f_input[1] >= 0 and rcoa_f_input[1] <= 1:
        R1_f = rcoa_f_input[1]
    elif rcoa_f_input[1] < 0:
        R1_f = round(rcoa_f_input[1] - int(rcoa_f_input[1]) + 1, 6)
    elif rcoa_f_input[1] > 1:
        R1_f = round(rcoa_f_input[1] - int(rcoa_f_input[1]), 6)

    if rcoa_f_input[2] >= 0 and rcoa_f_input[2] <= 1:
        R2_f = rcoa_f_input[2]
    elif rcoa_f_input[2] < 0:
        R2_f = round(rcoa_f_input[2] - int(rcoa_f_input[2]) + 1, 6)
    elif rcoa_f_input[2] > 1:
        R2_f = round(rcoa_f_input[2] - int(rcoa_f_input[2]), 6)

    rcoa_f_out=[R0_f,R1_f,R2_f]
    coa_f_out = rcoa_conv2_coa(rcoa_f_out,a,b,c)


    return coa_f_out,rcoa_f_out

# a function to annalysis the bonded connectivity with the input like this：
# con=[(0, 4), (1, 6), (2, 3), (3, 2), (3, 12), (4, 0), (4, 10)...]
def con_find(con):

    linker_index = []
    while len(con) > 0:
        akk = []
        aj = 0
        akk.append(con[aj][0])
        akk.append(con[aj][1])
        while aj < len(con):
            # print("***con[aj]:",con[aj])
            if (con[aj][0] in akk or con[aj][1] in akk):
                if (con[aj][0] in akk and con[aj][1] in akk):
                    akk.append(con[aj][0])
                    akk.append(con[aj][1])
                    con.pop(aj)  # improve the running speed
                else:
                    akk.append(con[aj][0])
                    akk.append(con[aj][1])
                    con.pop(aj)
                    aj = 0
            else:
                aj = aj + 1

        akk_set = list(set(akk))
        linker_index.append(akk_set)

    # print(linker_index)
    return linker_index  # generate a result like this: [[18,1,12,3,4],[5,6],[7,8,9,10]]

# convert Cartesian coordinates to Relative coordinates
def coa_conv2_rcoa(coa_input,acod_input,bcod_input,ccod_input):
    rcoa_output = [0,0,0]
    b1 = coa_input[0]
    b2 = coa_input[1]
    b3 = coa_input[2]

    a11 = acod_input[0]
    a12 = bcod_input[0]
    a13 = ccod_input[0]

    a21 = acod_input[1]
    a22 = bcod_input[1]
    a23 = ccod_input[1]

    a31 = acod_input[2]
    a32 = bcod_input[2]
    a33 = ccod_input[2]

    Len = a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32 - (a11 * a23 * a32 + a12 * a21 * a33 + a13 * a22 * a31)

    #
    xin = round((b1 * a22 * a33 + b2 * a13 * a32 + b3 * a12 * a23 - (b1 * a23 * a32 + b2 * a12 * a33 + b3 * a13 * a22)) / Len, 6)
    yin = round((b1 * a23 * a31 + b2 * a11 * a33 + b3 * a13 * a21 - (b1 * a21 * a33 + b2 * a13 * a31 + b3 * a11 * a23)) / Len, 6)
    zin = round((b1 * a21 * a32 + b2 * a12 * a31 + b3 * a11 * a22 - (b1 * a22 * a31 + b2 * a11 * a32 + b3 * a12 * a21)) / Len, 6)
    rcoa_output[0] = xin
    rcoa_output[1] = yin
    rcoa_output[2] = zin
    return rcoa_output

# convert Relative coordinates to Cartesian coordinates
def rcoa_conv2_coa(cord,acod_input,bcod_input,ccod_input):
    f_coordinate = [0, 0, 0]
    #print(acod, bcod, ccod)
    f_coordinate[0] = cord[0] * acod_input[0] + cord[1] * bcod_input[0] + cord[2] * ccod_input[0]
    f_coordinate[1] = cord[0] * acod_input[1] + cord[1] * bcod_input[1] + cord[2] * ccod_input[1]
    f_coordinate[2] = cord[0] * acod_input[2] + cord[1] * bcod_input[2] + cord[2] * ccod_input[2]
    return f_coordinate

# calculate the distance between two points with Cartesian coordinates
def distance_comput(a1, a2):
    x = a1[0] - a2[0]
    y = a1[1] - a2[1]
    z = a1[2] - a2[2]
    distance = math.sqrt(x*x+y*y+z*z)
    return distance

# Calculate the distance in a periodic structure with Cartesian coordinates as inputs;
#a_input,b_input and c_input represent the axis of the periodic structure
#x_n,y_n and z_n control the direction that should consider the structure priod
def p_distance_comput(coordination1, coordination2,a_input,b_input,c_input,x_n,y_n,z_n):
    p_d0 = distance_comput(coordination1, coordination2)
    rc1 = coa_conv2_rcoa(coordination1, a_input, b_input, c_input)
    rc2 = coa_conv2_rcoa(coordination2, a_input, b_input, c_input)
    rc2_0 = 0
    rc2_1 = 0
    rc2_2 = 0

    if x_n==1:
        if abs(rc1[0]- rc2[0]) > 0.5:
            if rc2[0] >= 0.5:
                rc2_0 = rc2[0] - 1
            elif rc2[0] < 0.5:
                rc2_0 = rc2[0] + 1
        else:
            rc2_0 = rc2[0]
    elif x_n==0:
        rc2_0 = rc2[0]

    if y_n==1:
        if abs(rc1[1]- rc2[1]) > 0.5:
            if rc2[1] >= 0.5:
                rc2_1 = rc2[1] - 1
            elif rc2[1] < 0.5:
                rc2_1 = rc2[1] + 1
        else:
            rc2_1 = rc2[1]
    elif y_n==0:
        rc2_1 = rc2[1]

    if z_n==1:
        if abs(rc1[2]- rc2[2]) > 0.5:
            if rc2[2] >= 0.5:
                rc2_2 = rc2[2] - 1
            elif rc2[2] < 0.5:
                rc2_2 = rc2[2] + 1
        else:
            rc2_2 = rc2[2]
    elif z_n==0:
        rc2_2 = rc2[2]

    coa2 = rcoa_conv2_coa([rc2_0, rc2_1, rc2_2], a_input, b_input, c_input)

    p_d1 = distance_comput(coordination1, coa2)


    d = min(p_d0, p_d1)
    # print(d)
    return d

#repeat N_of_cell times of the coordinates of a unit cell to that od supercell along vector_input
def supercell_coa(coa_input,N_of_cell,vector_input):
    coa_output = [i for i in coa_input]
    #print(coa_output)

    for ii in range (1,N_of_cell):
        for i_coa in coa_input:
            i_coa_out = [0, 0, 0]
            i_coa_out[0] = i_coa[0] + ii*vector_input[0]
            i_coa_out[1] = i_coa[1] + ii*vector_input[1]
            i_coa_out[2] = i_coa[2] + ii*vector_input[2]

            #print(i_coa)
            #print(i_coa_out)
            coa_output.append(i_coa_out)
            #print(coa_output)
    return  coa_output

#repeat struture units of the unit cell to the supercell
#for example: L_index_input=[[1,2,3,4,5], [4,5,6], ...] ,
# N_of_cell=5,
# length_all=6 #the value means the number of all the nonmetalic atoms in unit cell
def supercell_index(L_index_input,N_of_cell,length_all):
    L_index_output = [i for i in L_index_input]
    #print(L_index_input)

    for i in range(1, N_of_cell):
        for L_i in L_index_input:  #repeat L_i=[0,2,3,4,5] to [6,8,9,10,11],[12,14,15,16,17]...

            #print("L_i:",L_i)
            L_i_out =[Lii+i*length_all for Lii in L_i]
            #print("L_i_out",L_i_out)
            L_index_output.append (L_i_out)
    return  L_index_output


#repeat the z-coordinate of atoms in unit cell to that in supercell (repeat N_of_cell times);
#z_input is the z-coordinate of atoms in unit cell;h_input is the height of the unit cell
def supercell_Z(z_input,N_of_cell,h_input):  #[1.1, 2.2, 3.0, 5.9] ,vector_input=ccod
    z_output = [i for i in z_input]
    #print(z_output)

    for ii in range (1,N_of_cell):
        for i_z in z_input:
            i_z_out = i_z + ii*h_input
            #print(i_z)
            #print(i_z_out)
            z_output.append(i_z_out)
            #print(z_output)
    return z_output

#repeat a list N_of_cell times by adding a number of length for each element
#for example: input :([1,2,3],2,5)--->output: [1,2,3,6,7,8]
def supercell_list(list_input, N_of_cell, length):
    list_output = list_input

    for i in range(1, N_of_cell):
        list_output=list_output+[Lii + i * length for Lii in list_input]
    return list_output

#the input record the index of linkers that connected to each connected metals;
#the variable length means the number of the linkers
#the output record the index of metals that connected to each linker
def linker_metal(input,length):  #input=[[0,2],[0,1]] length=3; output:[[0,1],[1],[0]]
    index_array=np.arange(0,length)
    list_out=[]
    for i in index_array:
        list_out_each = []
        for j in range(0,len(input )):
            if i in input[j]:
                list_out_each.append(j)
            else: pass
        list_out.append(list_out_each)
    return list_out

# a very important function, used to estimate the pore diameters
# coa_input: the coordinates of atoms;
# a_input,b_input and c_input represent the axis of the periodic structure
#pld type = 1 means using projection along a-,b-,c- axis to estimate the pore diameter
#pld type = 2 means using projection exactly onto the facet plane  to estimate the pore diameter

def pld_comput(coa_input,acod_input,bcod_input,ccod_input,pld_type):
    Nf = 2  # N_for grid supercell
    if pld_type==1:
        point_o = [0, 0, 0]
        point_a = [Nf * acod_input[0], Nf * acod_input[1], Nf * acod_input[2]]
        point_b = [Nf * bcod_input[0], Nf * bcod_input[1], Nf * bcod_input[2]]
        axis_rB = bcod_input
        axis_rC = ccod_input
        axis_rA = acod_input

        ra_len = math.sqrt(axis_rA[0] * axis_rA[0] + axis_rA[1] * axis_rA[1] + axis_rA[2] * axis_rA[2])
        rb_len = math.sqrt(axis_rB[0] * axis_rB[0] + axis_rB[1] * axis_rB[1] + axis_rB[2] * axis_rB[2])
        rc_len = math.sqrt(axis_rC[0] * axis_rC[0] + axis_rC[1] * axis_rC[1] + axis_rC[2] * axis_rC[2])

        cos_bc = (axis_rB[0] * axis_rC[0] + axis_rB[1] * axis_rC[1] + axis_rB[2] * axis_rC[2]) / (rb_len * rc_len)
        cos_ac = (axis_rA[0] * axis_rC[0] + axis_rA[1] * axis_rC[1] + axis_rA[2] * axis_rC[2]) / (ra_len * rc_len)

        alpha_f = round(180 * math.acos(cos_bc) / 3.14159265, 6)
        beta_f = round(180 * math.acos(cos_ac) / 3.14159265, 6)
        if alpha_f < 55 and beta_f < 55:
            Nfc = 3
            point_c = [Nfc * ccod_input[0], Nfc * ccod_input[1], Nfc * ccod_input[2]]

        else:
            point_c = [Nf * ccod_input[0], Nf * ccod_input[1], Nf * ccod_input[2]]
            Nfc = Nf
        #
        # point_c = [Nf *ccod_input[0], Nf *ccod_input[1], Nf *ccod_input[2]]

        pld_coa_all = supercell_coa(supercell_coa(supercell_coa(coa_input, Nf, acod_input), Nf, bcod_input), Nfc,
                                    ccod_input)

        point_d = [point_a[0] + point_b[0], point_a[1] + point_b[1], point_a[2] + point_b[2]]  # 2(a+b)
        point_e = [point_a[0] + point_c[0], point_a[1] + point_c[1], point_a[2] + point_c[2]]
        point_f = [point_b[0] + point_c[0], point_b[1] + point_c[1], point_b[2] + point_c[2]]
        point_g = [point_a[0] + point_b[0] + point_c[0], point_a[1] + point_b[1] + point_c[1],
                   point_a[2] + point_b[2] + point_c[2]]
        point_x_min = min(point_a[0], point_b[0], point_c[0], point_d[0], point_e[0], point_f[0], point_g[0], point_o[0])
        point_x_max = max(point_a[0], point_b[0], point_c[0], point_d[0], point_e[0], point_f[0], point_g[0], point_o[0])
        point_y_min = min(point_a[1], point_b[1], point_c[1], point_d[1], point_e[1], point_f[1], point_g[1], point_o[1])
        point_y_max = max(point_a[1], point_b[1], point_c[1], point_d[1], point_e[1], point_f[1], point_g[1], point_o[1])
        point1 = [0.5 * (point_o[0] + point_c[0]), 0.5 * (point_o[1] + point_c[1])]
        point2 = [0.5 * (point_a[0] + point_e[0]), 0.5 * (point_a[1] + point_e[1])]
        point3 = [0.5 * (point_b[0] + point_f[0]), 0.5 * (point_b[1] + point_f[1])]
        point4 = [0.5 * (point_g[0] + point_d[0]), 0.5 * (point_g[1] + point_d[1])]

        grid_size = float(grid_size_input) #0.8  # unit: A
        grid_hy = (point_y_max - point_y_min)
        grid_Ny = math.ceil(grid_hy / grid_size)
        grid_hx = (point_x_max - point_x_min)
        grid_Nx = math.ceil(grid_hx / grid_size)

        # grid_hx = (point_x_max - point_x_min)
        # grid_Nx = 50
        # grid_size = grid_hx/grid_Nx
        # grid_hy = (point_y_max - point_y_min)
        # grid_Ny = math.ceil(grid_hy / grid_size)

        grid_ij = []    # [[x1,y1,1,20,19..],[x2,y2,2,5,8...]....]


        for i in range(0, grid_Ny):
            line_ij = []
            for j in range(0, grid_Nx):
                line_ij.append([point_x_min + (j + 0.5) * grid_size, (i + 0.5) * grid_size + point_y_min])
            grid_ij.append(line_ij)

        for k in range(0, len(pld_coa_all)):
            pld_rcoa_k = coa_conv2_rcoa(pld_coa_all[k], point_a, point_b, point_c)
            if pld_rcoa_k[0] > 1:
                pld_coa_all[k] = coa_conv2_rcoa([pld_rcoa_k[0] - 1, pld_rcoa_k[1], pld_rcoa_k[2]], point_a, point_b,
                                                point_c)
            if pld_rcoa_k[0] < 0:
                pld_coa_all[k] = coa_conv2_rcoa([pld_rcoa_k[0] + 1, pld_rcoa_k[1], pld_rcoa_k[2]], point_a, point_b,
                                                point_c)
            if pld_rcoa_k[1] > 1:
                pld_coa_all[k] = coa_conv2_rcoa([pld_rcoa_k[0], pld_rcoa_k[1] - 1, pld_rcoa_k[2]], point_a, point_b,
                                                point_c)
            if pld_rcoa_k[1] < 0:
                pld_coa_all[k] = coa_conv2_rcoa([pld_rcoa_k[0], pld_rcoa_k[1] + 1, pld_rcoa_k[2]], point_a, point_b,
                                                point_c)

            grid_i = int((pld_coa_all[k][1] - point_y_min) / grid_size)
            grid_j = int((pld_coa_all[k][0] - point_x_min) / grid_size)

            grid_ij[grid_i][grid_j].append(k)

    if pld_type==2:
        pld_coa_all = supercell_coa(supercell_coa(coa_input, Nf, bcod_input), Nf, ccod_input)

        point_a = acod_input
        point_b = [Nf * bcod_input[0], Nf * bcod_input[1], Nf * bcod_input[2]]
        point_c = [Nf * ccod_input[0], Nf * ccod_input[1], Nf * ccod_input[2]]
        point1 = point_a[1:]
        point2 = point_b[1:]
        point3 = point_c[1:]
        point4= [point2[0]+point3[0],point2[1]+point3[1]]

        point_x_min = min(point1[0], point2[0], point3[0], point4[0])
        point_x_max = max(point1[0], point2[0], point3[0], point4[0])
        point_y_min = min(point1[1], point2[1], point3[1], point4[1])
        point_y_max = max(point1[1], point2[1], point3[1], point4[1])

        grid_size = float(grid_size_input) #0.8  # unit: A
        grid_hy = (point_y_max - point_y_min)
        grid_Ny = math.ceil(grid_hy / grid_size)
        grid_hx = (point_x_max - point_x_min)
        grid_Nx = math.ceil(grid_hx / grid_size)
        # grid_hx = (point_x_max - point_x_min)
        # grid_Nx = 50
        # grid_size = grid_hx/grid_Nx
        # grid_hy = (point_y_max - point_y_min)
        # grid_Ny = math.ceil(grid_hy / grid_size)

        grid_ij = []  # [[x1,y1,1,20,19..],[x2,y2,2,5,8...]....]

        for i in range(0, grid_Ny):
            line_ij = []
            for j in range(0, grid_Nx):
                line_ij.append([point_x_min + (j + 0.5) * grid_size, (i + 0.5) * grid_size + point_y_min])


            grid_ij.append(line_ij)

        for k in range(0, len(pld_coa_all)):
            grid_i = int((pld_coa_all[k][2] - point_y_min) / grid_size)
            grid_j = int((pld_coa_all[k][1] - point_x_min) / grid_size)

            grid_ij[grid_i][grid_j].append(k)

    notnull_ij = []  #the grid cell has atoms
    pore_ij = [] #the hollow zone
    edge_ij = []  #record the edge grid of hollow zone
    edge2_ij = [] #record the edge grid of the projection
    all_ij = []  #all the grid cell

    pts = [point1, point2, point4, point3]
    polygon = geometry.Polygon(pts)

    for i in range(0, len(grid_ij)):
        for j in range(0, len(grid_ij[i])):
            all_ij.append(grid_ij[i][j][0:2])
            point = geometry.Point(grid_ij[i][j][0:2])
            line1 = geometry.LineString([point1, point2])
            line2 = geometry.LineString([point2, point4])
            line3 = geometry.LineString([point4, point3])
            line4 = geometry.LineString([point3, point1])
            edge_mark = 0
            # if len(grid_ij[i][j]) > 0:
            if polygon.covers(point) == True:
                if len(grid_ij[i][j]) > 2:  # means the presence of atoms
                    notnull_ij.append(grid_ij[i][j][0:2])
                    if i < len(grid_ij) - 1:

                        if len(grid_ij[i + 1][j]) == 2:
                            edge_ij.append(grid_ij[i][j][0:2])
                        if j < len(grid_ij[i]) - 1:
                            if len(grid_ij[i + 1][j + 1]) == 2:
                                edge_ij.append(grid_ij[i][j][0:2])
                        if j > 0:
                            if len(grid_ij[i + 1][j - 1]) == 2:
                                edge_ij.append(grid_ij[i][j][0:2])

                    if j < len(grid_ij[i]) - 1:
                        if len(grid_ij[i][j + 1]) == 2:
                            edge_ij.append(grid_ij[i][j][0:2])

                    if i > 0:
                        if len(grid_ij[i - 1][j]) == 2:
                            edge_ij.append(grid_ij[i][j][0:2])
                        if j < len(grid_ij[i]) - 1:
                            if len(grid_ij[i - 1][j + 1]) == 2:
                                edge_ij.append(grid_ij[i][j][0:2])
                        if j > 0:
                            if len(grid_ij[i - 1][j - 1]) == 2:
                                edge_ij.append(grid_ij[i][j][0:2])

                    if j > 0:
                        if len(grid_ij[i][j - 1]) == 2:
                            edge_ij.append(grid_ij[i][j][0:2])

                if line1.distance(point) < 1.3 * grid_size:
                    edge_ij.append(grid_ij[i][j][0:2])
                    edge2_ij.append(grid_ij[i][j][0:2])
                    edge_mark = 1

                if line2.distance(point) < 1.3 * grid_size:
                    edge_ij.append(grid_ij[i][j][0:2])
                    edge2_ij.append(grid_ij[i][j][0:2])
                    edge_mark = 1

                if line3.distance(point) < 1.3 * grid_size:
                    edge_ij.append(grid_ij[i][j][0:2])
                    edge2_ij.append(grid_ij[i][j][0:2])
                    edge_mark = 1

                if line4.distance(point) < 1.3 * grid_size:
                    edge_ij.append(grid_ij[i][j][0:2])
                    edge2_ij.append(grid_ij[i][j][0:2])
                    edge_mark = 1

                if len(grid_ij[i][j]) == 2 and not edge_mark == 1:  # means empty cell and not close to edge2
                    pore_ij.append(grid_ij[i][j])


    ld_ij = []
    ij_ld_ij=[]
    ld_ld_ij = []


    for el in edge_ij[:]:
        if edge_ij.count(el) > 1:
            edge_ij.remove(el)
    for el in pore_ij[:]:
        if pore_ij.count(el) > 1:
            print("clear please")

    if len(pore_ij)>0:
        for ij in range(0, len(pore_ij)):
            ld_ij_each = []
            ij_ld_ij_each = []
            ld_ld_ij_each = []
            for kl in range(0, len(edge_ij)):
                p_e_x = pore_ij[ij][0] - edge_ij[kl][0]
                p_e_y = pore_ij[ij][1] - edge_ij[kl][1]
                ij_lk_d = math.sqrt(p_e_x*p_e_x+p_e_y*p_e_y)
                ld_ij_each.append(ij_lk_d)
                ij_ld_ij_each.append(pore_ij[ij])
                ld_ld_ij_each.append(edge_ij[kl])
            min_ld_ij_each = min(ld_ij_each)
            min_N = ld_ij_each.index(min_ld_ij_each)
            if edge_ij[min_N] not in edge2_ij:
                ld_ij.append(min(ld_ij_each))
    elif len(pore_ij) == 0:
        ld_ij=[0]

    pld_out = 2 * max(ld_ij)

    return pld_out

#this function is used to rotate the crystal, making any facets become {001}
def lattice_rotation(facet_input,coa_input,acod_input,bcod_input,ccod_input):
    R_rcoa_all=[]

    if facet_input == "1 1 1":
        axis_rA = [bcod_input[0] - acod_input[0], bcod_input[1] - acod_input[1], bcod_input[2] - acod_input[2]]  # a  AB_rA
        axis_rB = [ccod_input[0] - acod_input[0], ccod_input[1] - acod_input[1], ccod_input[2] - acod_input[2]]  # b  AC _rB
        axis_rC = [0 - acod_input[0], 0 - acod_input[1], 0 - acod_input[2]]  # c  AO  rC

    # facet ="011"
    if facet_input == "0 1 1":
        axis_rA = [acod_input[0], acod_input[1], acod_input[2]]
        axis_rB = [ccod_input[0] - bcod_input[0], ccod_input[1] - bcod_input[1], ccod_input[2] - bcod_input[2]]
        axis_rC = [bcod_input[0], bcod_input[1], bcod_input[2]]

    # facet = "101"
    if facet_input == "1 0 1":
        axis_rA = [bcod_input[0], bcod_input[1], bcod_input[2]]  # OB
        axis_rB = [ccod_input[0] - acod_input[0], ccod_input[1] - acod_input[1], ccod_input[2] - acod_input[2]]  # CA
        axis_rC = [acod_input[0], acod_input[1], acod_input[2]]  # OA

    # facet == "110"
    if facet_input == "1 1 0":  #
        axis_rA = [ccod_input[0], ccod_input[1], ccod_input[2]]  # OC
        axis_rB = [bcod_input[0] - acod_input[0], bcod_input[1] - acod_input[1], bcod_input[2] - acod_input[2]]  # AB
        axis_rC = [acod_input[0], acod_input[1], acod_input[2]]  # OA

    # facet == "001"
    if facet_input == "0 0 1":
        axis_rA = [acod_input[0], acod_input[1], acod_input[2]]
        axis_rB = [bcod_input[0], bcod_input[1], bcod_input[2]]  # AB
        axis_rC = [ccod_input[0], ccod_input[1], ccod_input[2]]  # OA

    # facet == "010"
    if facet_input == "0 1 0":
        axis_rA = [acod_input[0], acod_input[1], acod_input[2]]
        axis_rB = [ccod_input[0], ccod_input[1], ccod_input[2]]
        axis_rC = [bcod_input[0], bcod_input[1], bcod_input[2]]

    # facet == "100"
    if facet_input == "1 0 0":
        axis_rA = [ccod_input[0], ccod_input[1], ccod_input[2]]
        axis_rB = [bcod_input[0], bcod_input[1], bcod_input[2]]
        axis_rC = [acod_input[0], acod_input[1], acod_input[2]]

    ra_len = math.sqrt(axis_rA[0] * axis_rA[0] + axis_rA[1] * axis_rA[1] + axis_rA[2] * axis_rA[2])
    rb_len = math.sqrt(axis_rB[0] * axis_rB[0] + axis_rB[1] * axis_rB[1] + axis_rB[2] * axis_rB[2])
    rc_len = math.sqrt(axis_rC[0] * axis_rC[0] + axis_rC[1] * axis_rC[1] + axis_rC[2] * axis_rC[2])

    cos_bc = (axis_rB[0] * axis_rC[0] + axis_rB[1] * axis_rC[1] + axis_rB[2] * axis_rC[2]) / (rb_len * rc_len)
    cos_ac = (axis_rA[0] * axis_rC[0] + axis_rA[1] * axis_rC[1] + axis_rA[2] * axis_rC[2]) / (ra_len * rc_len)
    cos_ab = (axis_rA[0] * axis_rB[0] + axis_rA[1] * axis_rB[1] + axis_rA[2] * axis_rB[2]) / (ra_len * rb_len)

    alpha_out = round(180 * math.acos(cos_bc) / 3.14159265, 6)
    beta_out = round(180 * math.acos(cos_ac) / 3.14159265, 6)
    gama_out = round(180 * math.acos(cos_ab) / 3.14159265, 6)

    for i in range(0, len(coa_input)):
        R_rcoa_each = coa_conv2_rcoa(coa_input[i], axis_rA, axis_rB, axis_rC)
        if R_rcoa_each[0] >= 0 and R_rcoa_each[0] <= 1:
            R0 = R_rcoa_each[0]
        elif R_rcoa_each[0] < 0:
            R0 = round(R_rcoa_each[0] - int(R_rcoa_each[0]) + 1, 6)
        elif R_rcoa_each[0] > 1:
            R0 = round(R_rcoa_each[0] - int(R_rcoa_each[0]), 6)

        if R_rcoa_each[1] >= 0 and R_rcoa_each[1] <= 1:
            R1 = R_rcoa_each[1]
        elif R_rcoa_each[1] < 0:
            R1 = round(R_rcoa_each[1] - int(R_rcoa_each[1]) + 1, 6)
        elif R_rcoa_each[1] > 1:
            R1 = round(R_rcoa_each[1] - int(R_rcoa_each[1]), 6)

        if R_rcoa_each[2] >= 0 and R_rcoa_each[2] <= 1:
            R2 = R_rcoa_each[2]
        elif R_rcoa_each[2] < 0:
            R2 = round(R_rcoa_each[2] - int(R_rcoa_each[2]) + 1, 6)
        elif R_rcoa_each[2] > 1:
            R2 = round(R_rcoa_each[2] - int(R_rcoa_each[2]), 6)

        R_rcoa_all.append([R0, R1, R2])

    R_A = cos_bc
    R_B = cos_ac
    R_C = cos_ab

    R_ca = 3.14159265 * gama_out / 180
    R_h = rc_len * math.sqrt(1 - R_A * R_A - R_B * R_B - R_C * R_C + 2 * R_A * R_B * R_C) / math.sin(R_ca)

    R_acod = [ra_len, 0, 0]
    R_bcod = [rb_len * R_C, rb_len * math.sin(R_ca), 0]
    R_ccodx = rc_len * R_B
    R_ccodz = R_h
    R_ccody = (R_A * rb_len * rc_len - R_bcod[0] * R_ccodx - R_bcod[2] * R_ccodz) / R_bcod[1]
    R_ccod = [R_ccodx, R_ccody, R_ccodz]

    R_coa_all = []
    for i in range(0, len(R_rcoa_all)):
        R_coa_all.append(rcoa_conv2_coa(R_rcoa_all[i], R_acod, R_bcod, R_ccod))
    #print("R_coa_all",R_coa_all)
    return R_coa_all, R_acod, R_bcod, R_ccod,[ra_len,rb_len,rc_len, alpha_out, beta_out, gama_out]


####################Tkinker##################################

def button_click():
    global running  # create global
    global slab_h
    global vacuum
    global facet_index_input
    global grid_size_input
    global file_dir
    global job_type  # 0 means cleave surface; 1 means supercell expansion
    global min_X_set
    global min_Y_set
    global min_X_length
    global min_Y_length
    global start_t


    running = True
    slab_h = float(expression3.get())
    vacuum = float(expression4.get())
    facet_index_input = expression7.get()
    grid_size_input = expression7_1.get()
    min_X_length=float(expression5.get())
    min_Y_length = float(expression6.get())

    # Create new thread
    t = Thread(target=start)
    # Start new thread
    start_t = time.time()
    t.start()


##the function to cleave surface
def start():
    log_lines = ["MOF-Membrane Constructor-Job Log\n\n"]

    def slab_construction(elemt_all_input, coa_all_input, acod_input, bcod_input, ccod_input, lattice_input, yn):

        alen_input = lattice_input[0]
        blen_input = lattice_input[1]
        clen_input = lattice_input[2]
        aang_input = lattice_input[3]
        bang_input = lattice_input[4]
        cang_input = lattice_input[5]
        h = ccod_input[2]
        print("h",h)
        f_break_atm = []

        if yn == 0:  # AUTO mode
            nonmetal_elemt = []  # elements without metal
            metal_elemt = []
            nonmetal_coa = []  # coordinate of atom, no metals
            metal_coa = []
            nonmetal_charge = []
            metal_charge = []

            #print("start metal recognition:")
            for i in range(0, len(elemt_all_input)):  # recognize metals
                if elemt_all_input[i] not in nonmetal_list:
                    metal_elemt.append(elemt_all_input[i])
                    metal_coa.append(coa_all_input[i])
                    if charge_r == 1:
                        metal_charge.append(charge_all_N[i])
                else:
                    nonmetal_elemt.append(elemt_all_input[i])
                    nonmetal_coa.append(coa_all_input[i])
                    if charge_r == 1:
                        nonmetal_charge.append(charge_all_N[i])

            print("atom number of nonmetal:", len(nonmetal_elemt))
            print("atom number of metal:", len(metal_elemt))
            len_nonmetal_elemt_u = len(nonmetal_elemt)

            # it is expected to get a value like this:
            # f_con=[(0, 4), (1, 6), (2, 3), (3, 2), (3, 12), (4, 0), (4, 10)...]
            f_con = []

            for i in range(0, len(nonmetal_elemt)):

                kk = []  # record the index of atoms connected to atom i
                for j in range(i, len(nonmetal_elemt)):
                    if j != i:
                        d0 = distance_comput(nonmetal_coa[i], nonmetal_coa[j])

                        nonmetal_rcoai = coa_conv2_rcoa(nonmetal_coa[i], acod_input, bcod_input, ccod_input)
                        nonmetal_rcoaj = coa_conv2_rcoa(nonmetal_coa[j], acod_input, bcod_input, ccod_input)

                        z_break_mark = 0
                        jrcoa0 = 0
                        jrcoa1 = 0
                        jrcoa2 = 0
                        if abs(nonmetal_rcoai[0] - nonmetal_rcoaj[0]) > 0.5:
                            if nonmetal_rcoaj[0] >= 0.5:
                                jrcoa0 = nonmetal_rcoaj[0] - 1
                            elif nonmetal_rcoaj[0] < 0.5:
                                jrcoa0 = nonmetal_rcoaj[0] + 1
                        else:
                            jrcoa0 = nonmetal_rcoaj[0]

                        if abs(nonmetal_rcoai[1] - nonmetal_rcoaj[1]) > 0.5:
                            if nonmetal_rcoaj[1] >= 0.5:
                                jrcoa1 = nonmetal_rcoaj[1] - 1
                            elif nonmetal_rcoaj[1] < 0.5:
                                jrcoa1 = nonmetal_rcoaj[1] + 1
                        else:
                            jrcoa1 = nonmetal_rcoaj[1]

                        if abs(nonmetal_rcoai[2] - nonmetal_rcoaj[2]) > 0.5:
                            if nonmetal_rcoaj[2] >= 0.5:
                                jrcoa2 = nonmetal_rcoaj[2] - 1
                            elif nonmetal_rcoaj[2] < 0.5:
                                jrcoa2 = nonmetal_rcoaj[2] + 1
                            z_break_mark = 1
                        else:
                            jrcoa2 = nonmetal_rcoaj[2]

                        coaj = rcoa_conv2_coa([jrcoa0, jrcoa1, jrcoa2], acod_input, bcod_input, ccod_input)
                        d1 = distance_comput(nonmetal_coa[i], coaj)
                        d = min(d0, d1)

                        bond = float(atom_r[nonmetal_elemt[i]]) + float(atom_r[nonmetal_elemt[j]])

                        if d <= bond:
                            kk.append(j)
                            if d == d1 and z_break_mark == 1:
                                f_break_atm.append(i)

                if len(kk) > 0:
                    for ii in kk:
                        con_single = (i, ii)
                        f_con.append(con_single)
                elif len(kk) == 0:
                    pass

            L_index = con_find(f_con)  # The index here is the index of "nonmetal_elemt" list

            # the following is performed analyze the connectivity between metal and linker;
            # the linker does not include the single nonmetal atoms
            metal_connect_structure = []  # [[221, 232, 234, 243], [228, 230, 241, 251],....]
            L_index_list = []  # [1,2,9,11...]

            for i in range(0, len(metal_elemt)):
                metal_connect_structure_atoms = []
                for linker_i in range(0, len(L_index)):
                    for linker_i_atomindex in L_index[linker_i]:
                        L_index_list.append(linker_i_atomindex)
                        d = p_distance_comput(metal_coa[i], nonmetal_coa[linker_i_atomindex], acod_input, bcod_input,
                                              ccod_input, 1, 1, 1)
                        m_bond = float(atom_r[metal_elemt[i]]) + float(atom_r[nonmetal_elemt[linker_i_atomindex]])
                        if d <= m_bond:
                            metal_connect_structure_atoms.append(linker_i_atomindex)
                        else:
                            pass
                # print(metal_connect_structure_atoms)
                metal_connect_structure.append(metal_connect_structure_atoms)
        elif yn == 1:  # "SET" mode
            nonmetal_elemt = []  # elements without metal
            metal_elemt = []
            nonmetal_coa = []  # coordinate of atom, no metals
            metal_coa = []
            nonmetal_charge = []
            metal_charge = []
            #print("start metal recognition:")

            for i in range(0, len(elemt_all_input)):  # recognize metals
                if elemt_all_input[i] not in nonmetal_list:
                    metal_elemt.append(elemt_all_input[i])
                    metal_coa.append(coa_all_input[i])
                    if charge_r == 1:
                        metal_charge.append(charge_all_N[i])
                else:
                    nonmetal_elemt.append(elemt_all_input[i])
                    nonmetal_coa.append(coa_all_input[i])
                    if charge_r == 1:
                        nonmetal_charge.append(charge_all_N[i])
            f_con1 = []
            f_con2 = []
            len_nonmetal_elemt_u = len(nonmetal_elemt_u)

            for n in range(0, N2):
                for i in f_con0:
                    f_con1.append((i[0] + n * len(nonmetal_elemt_u), i[1] + n * len_nonmetal_elemt_u))
            # print("f_con1", f_con1)
            for i in f_con1:
                d = p_distance_comput(nonmetal_coa[i[0]], nonmetal_coa[i[1]], acod_input, bcod_input, ccod_input, 1, 1,
                                      1)
                bond_d = float(atom_r[nonmetal_elemt[i[0]]]) + float(atom_r[nonmetal_elemt[i[1]]])
                if d <= bond_d:
                    f_con2.append(i)
                else:
                    bx = [k for k in range(0, len(nonmetal_coa)) if (k - i[1]) % len_nonmetal_elemt_u == 0]
                    mark = 0
                    for j in bx:
                        d = p_distance_comput(nonmetal_coa[i[0]], nonmetal_coa[j], acod_input, bcod_input, ccod_input,
                                              1, 1, 1)
                        bond_d = float(atom_r[nonmetal_elemt[i[0]]]) + float(atom_r[nonmetal_elemt[j]])
                        if d <= bond_d:
                            f_con2.append((i[0], j))
                            mark = 1
                            break
                        else:
                            pass
                    if mark == 0:
                        print("error")
                        error_lines.append(Filename+':  '+"failed to build bond connectivity"+"\n")
            # print("f_con2", f_con2)
            L_index = con_find(f_con2)

            print("start get break atoms")
            ##mark the atom pairs that across the top/bottom periodic boundary;

            for k in L_index:
                for i in k:
                    for j in k:
                        if j != i:
                            d0 = distance_comput(nonmetal_coa[i], nonmetal_coa[j])
                            nonmetal_rcoai = coa_conv2_rcoa(nonmetal_coa[i], acod_input, bcod_input, ccod_input)
                            nonmetal_rcoaj = coa_conv2_rcoa(nonmetal_coa[j], acod_input, bcod_input, ccod_input)
                            z_break_mark = 0
                            jrcoa0 = 0
                            jrcoa1 = 0
                            jrcoa2 = 0
                            if abs(nonmetal_rcoai[0] - nonmetal_rcoaj[0]) > 0.5:
                                if nonmetal_rcoaj[0] >= 0.5:
                                    jrcoa0 = nonmetal_rcoaj[0] - 1
                                elif nonmetal_rcoaj[0] < 0.5:
                                    jrcoa0 = nonmetal_rcoaj[0] + 1
                            else:
                                jrcoa0 = nonmetal_rcoaj[0]

                            if abs(nonmetal_rcoai[1] - nonmetal_rcoaj[1]) > 0.5:
                                if nonmetal_rcoaj[1] >= 0.5:
                                    jrcoa1 = nonmetal_rcoaj[1] - 1
                                elif nonmetal_rcoaj[1] < 0.5:
                                    jrcoa1 = nonmetal_rcoaj[1] + 1
                            else:
                                jrcoa1 = nonmetal_rcoaj[1]

                            if abs(nonmetal_rcoai[2] - nonmetal_rcoaj[2]) > 0.5:
                                if nonmetal_rcoaj[2] >= 0.5:
                                    jrcoa2 = nonmetal_rcoaj[2] - 1
                                elif nonmetal_rcoaj[2] < 0.5:
                                    jrcoa2 = nonmetal_rcoaj[2] + 1
                                z_break_mark = 1
                            else:
                                jrcoa2 = nonmetal_rcoaj[2]

                            coaj = rcoa_conv2_coa([jrcoa0, jrcoa1, jrcoa2], acod_input, bcod_input, ccod_input)
                            d1 = distance_comput(nonmetal_coa[i], coaj)
                            d = min(d0, d1)
                            bond = float(atom_r[nonmetal_elemt[i]]) + float(atom_r[nonmetal_elemt[j]])
                            if d <= bond and d == d1 and z_break_mark == 1:
                                f_break_atm.append(i)
            f_break_atm = list(set(f_break_atm))

            # the following is performed analyze the connectivity between metal and linker;
            # the linker does not include the single nonmetal atoms

            metal_connect_structure = []  # [[221, 232, 234, 243], [228, 230, 241, 251],....]
            L_index_list = []  # [1,2,9,11...]
            metal_connect_structure = supercell_index(metal_connect_structure_u, N2, len_nonmetal_elemt_u)
            # The metal_connect_structure here is not exact, mainly used to narrow the range

            for linker_i in range(0, len(L_index)):
                for linker_i_atomindex in L_index[linker_i]:
                    L_index_list.append(linker_i_atomindex)

        # the following is performed analyze the connectivity between metal and single nonmetal atoms
        nonmetal_elemt_list = np.arange(0, len(nonmetal_elemt))
        single_atom_list = [i for i in nonmetal_elemt_list if
                            i not in L_index_list]  # 其中的数字对应的是nonmetal_elemt和nonmetal_coa 的index
        # print(single_atom_list)

        single_atom_elemt = []
        single_atom_connect = []
        for i in single_atom_list:
            single_atom_elemt.append(nonmetal_elemt[i])
            n_connect = 0
            for j in range(0,len(metal_coa)):
                d_single = p_distance_comput(metal_coa[j], nonmetal_coa[i], acod_input, bcod_input,
                                              ccod_input, 1, 1, 1)
                bond_single = float(atom_r[nonmetal_elemt[i]]) + float(atom_r[metal_elemt[j]])

                if d_single<= bond_single:
                    n_connect = n_connect+1
            single_atom_connect.append(n_connect)
        #print("single_atom_connect",single_atom_connect)

        #####
        # Regularize the structure of the linker that across the top/bottom periodic boundary;
        # connect it into a whole, reduce the computing cost in the later stage;
        # and get its lowest point and highest point at the same time
        L_Z_min = []
        L_Z_max = []
        f_break_atm = list(set(f_break_atm))
        break_linker = []  #
        atom_linker_atm = []
        atom_linker_L = []

        for i in range(0, len(L_index)):  # supercell_L_index[i]= [2,3,5,9]
            for j in L_index[i]:
                # supercell_atom_linker.append([j,i])  #[2,i],[3,i],[5,i],[9,i]
                atom_linker_atm.append(j)
                atom_linker_L.append(i)
        for i in f_break_atm:
            L_k = atom_linker_atm.index(i)
            break_linker.append(atom_linker_L[L_k])

        break_linker = list(set(break_linker))
        # print("len(L_index) ", len(L_index))
        # print("break_linker = ", break_linker)

        for linker_i in range(0, len(L_index)):
            denote = 0
            L_Z_each = []
            for linker_i_atomindex in L_index[linker_i]:
                L_Z_each.append(nonmetal_coa[linker_i_atomindex][2])

            L_Z_each.sort()

            z_mark = 0
            z_mark0 = 0
            z_break = 0
            if linker_i in break_linker:

                for i in range(0, len(L_Z_each) - 1):
                    d3 = L_Z_each[i + 1] - L_Z_each[i]
                    if d3 > 3.01:  #
                        z_mark = L_Z_each[i + 1]
                        z_mark0 = L_Z_each[i]
                        z_break = 1
                        break

                if z_break == 1 and z_mark <= 0.5 * h:
                    for linker_i_atomindex in L_index[linker_i]:
                        # print("***", nonmetal_coa[linker_i_atomindex], "****")
                        if nonmetal_coa[linker_i_atomindex][2] < z_mark:
                            # print("***", nonmetal_coa[linker_i_atomindex], "****")
                            nonmetal_coa[linker_i_atomindex][0] = nonmetal_coa[linker_i_atomindex][0] + ccod_input[0]
                            nonmetal_coa[linker_i_atomindex][1] = nonmetal_coa[linker_i_atomindex][1] + ccod_input[1]
                            nonmetal_coa[linker_i_atomindex][2] = nonmetal_coa[linker_i_atomindex][2] + ccod_input[2]

                            L_Z_each[0] = z_mark
                            L_Z_each[-1] = z_mark0 + ccod_input[2]

                        else:
                            pass
                        # print("####", nonmetal_coa[linker_i_atomindex], "###")
                if z_break == 1 and z_mark > 0.5 * h:
                    for linker_i_atomindex in L_index[linker_i]:
                        # print("***", nonmetal_coa[linker_i_atomindex], "****")
                        if nonmetal_coa[linker_i_atomindex][2] >= z_mark:
                            # print("***", nonmetal_coa[linker_i_atomindex], "****")
                            nonmetal_coa[linker_i_atomindex][0] = nonmetal_coa[linker_i_atomindex][0] - ccod_input[0]
                            nonmetal_coa[linker_i_atomindex][1] = nonmetal_coa[linker_i_atomindex][1] - ccod_input[1]
                            nonmetal_coa[linker_i_atomindex][2] = nonmetal_coa[linker_i_atomindex][2] - ccod_input[2]

                            L_Z_each[0] = z_mark - ccod_input[2]
                            L_Z_each[-1] = z_mark0
                        else:
                            pass
                        # print("####", nonmetal_coa[linker_i_atomindex], "###")

                elif z_break == 0:
                    #
                    ###########################
                    linker_f_con = []
                    linker_single = []
                    linker_part_index = []

                    break_atoms_each = []

                    fx = 0  #
                    for l in L_index[linker_i]:  # i-->l
                        kkk = []  # kk-->kkk
                        for m in L_index[linker_i]:  # j-->m  range(i-->rang(0
                            if m != l:
                                d00 = distance_comput(nonmetal_coa[l], nonmetal_coa[m])  # d0-->d00
                                # print(d)
                                bond00 = float(atom_r[nonmetal_elemt[m]]) + float(
                                    atom_r[nonmetal_elemt[l]])  # bond-->bond00

                                nonmetal_rcoal = coa_conv2_rcoa(nonmetal_coa[l], acod_input, bcod_input, ccod_input)
                                nonmetal_rcoam = coa_conv2_rcoa(nonmetal_coa[m], acod_input, bcod_input, ccod_input)

                                mrcoa0, mrcoa1, mrcoa2, mrcoa22 = 0, 0, 0, 0
                                if abs(nonmetal_rcoal[0] - nonmetal_rcoam[0]) > 0.5:
                                    # print("4444444444")
                                    if nonmetal_rcoam[0] >= 0.5:
                                        mrcoa0 = nonmetal_rcoam[0] - 1
                                    elif nonmetal_rcoam[0] < 0.5:
                                        mrcoa0 = nonmetal_rcoam[0] + 1
                                else:
                                    mrcoa0 = nonmetal_rcoam[0]

                                if abs(nonmetal_rcoal[1] - nonmetal_rcoam[1]) > 0.5:
                                    if nonmetal_rcoam[1] >= 0.5:
                                        mrcoa1 = nonmetal_rcoam[1] - 1
                                    elif nonmetal_rcoam[1] < 0.5:
                                        mrcoa1 = nonmetal_rcoam[1] + 1
                                    # print("55555555")
                                else:
                                    mrcoa1 = nonmetal_rcoam[1]

                                if abs(nonmetal_rcoal[2] - nonmetal_rcoam[2]) > 0.5:
                                    if nonmetal_rcoam[2] >= 0.5:
                                        mrcoa22 = nonmetal_rcoam[2] - 1
                                    elif nonmetal_rcoam[2] < 0.5:
                                        mrcoa22 = nonmetal_rcoam[2] + 1

                                else:
                                    mrcoa22 = nonmetal_rcoam[2]

                                mrcoa2 = nonmetal_rcoam[2]
                                coam = rcoa_conv2_coa([mrcoa0, mrcoa1, mrcoa2], acod_input, bcod_input, ccod_input)
                                # print(rcoa[j],"----",jrcoa0,jrcoa1,jrcoa2)
                                # print(rcoa[i])
                                d11 = distance_comput(nonmetal_coa[l], coam)

                                coam_break = rcoa_conv2_coa([mrcoa0, mrcoa1, mrcoa22], acod_input, bcod_input,
                                                            ccod_input)
                                d11_break = distance_comput(nonmetal_coa[l], coam_break)

                                dml = min(d00, d11)
                                dml_break = min(d00, d11_break)

                                if dml <= bond00:
                                    kkk.append(m)
                                    # print("no break")
                                elif dml > bond00 and dml_break <= bond00:
                                    break_atoms_each.append(m)
                                    break_atoms_each.append(l)

                        if len(kkk) > 0:
                            for ll in kkk:
                                con_single = (l, ll)
                                linker_f_con.append(con_single)
                        elif len(kkk) == 0:
                            linker_single.append(l)

                    break_atoms_each = list(set(break_atoms_each))

                    linker_part_index = con_find(linker_f_con)
                    # print("linker_part_index:", linker_part_index)

                    for s in linker_single:
                        linker_part_index.append([s])
                    len_part = []
                    for eachpart in linker_part_index:
                        len_part.append(len(eachpart))
                    len_max_index = len_part.index(max(len_part))

                    ll_z_reserve = []
                    ll_z_linker = []

                    for ll_i in linker_part_index[len_max_index]:
                        ll_z_reserve.append(nonmetal_coa[ll_i][2])
                        ll_z_linker.append(nonmetal_coa[ll_i][2])

                    ll_z_reserve.sort()
                    # ll_middle = (ll_z_reserve[0] + ll_z_reserve[-1]) / 2

                    break_part_index = []
                    # break_max = []
                    for eachpart_i in range(0, len(linker_part_index)):
                        break_part_index.append([i for i in linker_part_index[eachpart_i] if i in break_atoms_each])

                    len_break_part_index = 0
                    for i in break_part_index:
                        len_break_part_index = len_break_part_index + len(i)

                    F_break = [i for i in break_part_index[len_max_index]]
                    F_linker = [i for i in linker_part_index[len_max_index]]
                    break_part_index.pop(len_max_index)
                    linker_part_index.pop(len_max_index)

                    while len(break_part_index) > 0:
                        break_mark = 0
                        for i in range(0, len(break_part_index)):
                            if break_mark == 1:
                                break
                            if break_mark == 0:
                                for p in break_part_index[i]:
                                    # if r in F_break can be found that is connected to p in  break_part_index[i]?
                                    # yes: break_part_index[i] move from  break_part_index to F_break,linker_part_index[i] move from linker_part_index to F_linker
                                    # yes: modify the nonmetal_coa of  linker_part_index[i]
                                    # yes: "break" to search the next, i.e., break_part_index[i+1]
                                    # no: pass
                                    yn_mark = 0
                                    fx = 0
                                    for r in F_break:
                                        d_break = distance_comput(nonmetal_coa[r], nonmetal_coa[p])
                                        bond_rp = float(atom_r[nonmetal_elemt[r]]) + float(
                                            atom_r[nonmetal_elemt[p]])  # bond-->bond00

                                        nonmetal_rcoa_r = coa_conv2_rcoa(nonmetal_coa[r], acod_input, bcod_input,
                                                                         ccod_input)
                                        nonmetal_rcoa_p = coa_conv2_rcoa(nonmetal_coa[p], acod_input, bcod_input,
                                                                         ccod_input)

                                        prcoa0, prcoa1, prcoa2, prcoa3, prcoa4 = 0, 0, 0, 0, 0
                                        if abs(nonmetal_rcoa_r[0] - nonmetal_rcoa_p[0]) > 0.5:
                                            if nonmetal_rcoa_p[0] >= 0.5:
                                                prcoa0 = nonmetal_rcoa_p[0] - 1
                                            elif nonmetal_rcoa_p[0] < 0.5:
                                                prcoa0 = nonmetal_rcoa_p[0] + 1
                                        else:
                                            prcoa0 = nonmetal_rcoa_p[0]

                                        if abs(nonmetal_rcoa_r[1] - nonmetal_rcoa_p[1]) > 0.5:
                                            if nonmetal_rcoa_p[1] >= 0.5:
                                                prcoa1 = nonmetal_rcoa_p[1] - 1
                                            elif nonmetal_rcoa_p[1] < 0.5:
                                                prcoa1 = nonmetal_rcoa_p[1] + 1
                                        else:
                                            prcoa1 = nonmetal_rcoa_p[1]

                                        prcoa2_0 = nonmetal_rcoa_p[2]
                                        prcoa2_1 = nonmetal_rcoa_p[2] + 1
                                        prcoa2_2 = nonmetal_rcoa_p[2] + 2
                                        prcoa2_m1 = nonmetal_rcoa_p[2] - 1
                                        prcoa2_m2 = nonmetal_rcoa_p[2] - 2

                                        coa_p0 = rcoa_conv2_coa([prcoa0, prcoa1, prcoa2_0], acod_input, bcod_input,
                                                                ccod_input)
                                        d_break_0 = distance_comput(nonmetal_coa[r], coa_p0)

                                        coa_p1 = rcoa_conv2_coa([prcoa0, prcoa1, prcoa2_1], acod_input, bcod_input,
                                                                ccod_input)
                                        d_break_1 = distance_comput(nonmetal_coa[r], coa_p1)

                                        coa_p2 = rcoa_conv2_coa([prcoa0, prcoa1, prcoa2_2], acod_input, bcod_input,
                                                                ccod_input)
                                        d_break_2 = distance_comput(nonmetal_coa[r], coa_p2)

                                        coa_pm1 = rcoa_conv2_coa([prcoa0, prcoa1, prcoa2_m1], acod_input, bcod_input,
                                                                 ccod_input)
                                        d_break_m1 = distance_comput(nonmetal_coa[r], coa_pm1)

                                        coa_pm2 = rcoa_conv2_coa([prcoa0, prcoa1, prcoa2_m2], acod_input, bcod_input,
                                                                 ccod_input)
                                        d_break_m2 = distance_comput(nonmetal_coa[r], coa_pm2)

                                        d_break_final = min(d_break, d_break_0, d_break_1, d_break_2, d_break_m1,
                                                            d_break_m2)

                                        if d_break_final > bond_rp:
                                            pass
                                        elif d_break_final < bond_rp and d_break_final == d_break:
                                            fx = 0
                                            yn_mark = 1
                                            break
                                        elif d_break_final < bond_rp and d_break_final == d_break_0:
                                            fx = 0
                                            yn_mark = 1
                                            break

                                        elif d_break_final < bond_rp and d_break_final == d_break_1:
                                            fx = 1
                                            yn_mark = 1
                                            break

                                        elif d_break_final < bond_rp and d_break_final == d_break_2:
                                            fx = 2
                                            yn_mark = 1
                                            break

                                        elif d_break_final < bond_rp and d_break_final == d_break_m1:
                                            fx = -1
                                            yn_mark = 1
                                            break

                                        elif d_break_final < bond_rp and d_break_final == d_break_m2:
                                            fx = -2
                                            yn_mark = 1
                                            break
                                    if yn_mark == 0:  # means no
                                        pass

                                    elif yn_mark == 1:
                                        F_break = F_break + break_part_index[i]
                                        F_linker = F_linker + linker_part_index[i]
                                        # print("p,break_part_index[i] before pop", p, break_part_index[i])
                                        # modify the nonmetal_coa of  linker_part_index[i]
                                        if fx == 0:
                                            break_part_index.pop(i)
                                            linker_part_index.pop(i)
                                        elif fx == -1:
                                            for j in linker_part_index[i]:
                                                nonmetal_coa[j][0] = nonmetal_coa[j][0] - ccod_input[0]
                                                nonmetal_coa[j][1] = nonmetal_coa[j][1] - ccod_input[1]
                                                nonmetal_coa[j][2] = nonmetal_coa[j][2] - ccod_input[2]

                                                ll_z_linker.append(nonmetal_coa[j][2])
                                            break_part_index.pop(i)
                                            linker_part_index.pop(i)

                                        elif fx == -2:
                                            for j in linker_part_index[i]:
                                                nonmetal_coa[j][0] = nonmetal_coa[j][0] - 2 * ccod_input[0]
                                                nonmetal_coa[j][1] = nonmetal_coa[j][1] - 2 * ccod_input[1]
                                                nonmetal_coa[j][2] = nonmetal_coa[j][2] - 2 * ccod_input[2]

                                                ll_z_linker.append(nonmetal_coa[j][2])
                                            break_part_index.pop(i)
                                            linker_part_index.pop(i)

                                        elif fx == 1:
                                            for j in linker_part_index[i]:
                                                nonmetal_coa[j][0] = nonmetal_coa[j][0] + ccod_input[0]
                                                nonmetal_coa[j][1] = nonmetal_coa[j][1] + ccod_input[1]
                                                nonmetal_coa[j][2] = nonmetal_coa[j][2] + ccod_input[2]

                                                ll_z_linker.append(nonmetal_coa[j][2])
                                            break_part_index.pop(i)
                                            linker_part_index.pop(i)

                                        elif fx == 2:
                                            for j in linker_part_index[i]:
                                                nonmetal_coa[j][0] = nonmetal_coa[j][0] + 2 * ccod_input[0]
                                                nonmetal_coa[j][1] = nonmetal_coa[j][1] + 2 * ccod_input[1]
                                                nonmetal_coa[j][2] = nonmetal_coa[j][2] + 2 * ccod_input[2]
                                                ll_z_linker.append(nonmetal_coa[j][2])
                                            break_part_index.pop(i)
                                            linker_part_index.pop(i)
                                        break_mark = 1
                                        break
                    # here get break_max #print("break_max", break_max)
                    # here get F_break #print("F_break", F_break)
                    # here get linker_single #print("linker_single",linker_single)
                    ll_z_linker.sort()
                    L_Z_each[0] = ll_z_linker[0]
                    L_Z_each[-1] = ll_z_linker[-1]

            L_Z_min.append(round(L_Z_each[0], 4))
            L_Z_max.append(round(L_Z_each[-1], 4))

        #
        nonmetal_length = len(nonmetal_elemt)
        NOC = math.ceil((slab_h + vacuum) / h + 2)  # number of unit-cell in the super-cell;
        # NOC =1

        # print("h:", h, "  NOC:", NOC)
        # print("length:", len(nonmetal_elemt))

        # the followings is performed to ananlysis the connectivities in supercell

        supercell_acod = acod_input
        supercell_bcod = bcod_input
        supercell_ccod = [NOC * ccod_input[0], NOC * ccod_input[1], NOC * ccod_input[2]]

        supercell_nonmetal_elemt = NOC * nonmetal_elemt  # [C,H,C,H,O]-->[C,H,C,H,O,  C,H,C,H,O,...]
        supercell_metal_elemt = NOC * metal_elemt  # [Cu,Zn]-->[Cu,Zn,  Cu,Zn...]

        supercell_nonmetal_coa = supercell_coa(nonmetal_coa, NOC,
                                               ccod_input)  # [[1.0,2.0,3.0],[9.0,8.0,6.0]]-->[[1.0,2.0,3.0],[9.0,8.0,6.0], [1.0+x,2.0+y,3.0+z],[9.0+x,8.0+y,6.0+z]...]
        supercell_metal_coa = supercell_coa(metal_coa, NOC, ccod_input)
        supercell_L_Z_min = supercell_Z(L_Z_min, NOC, round(h, 4))  # [-0.9,2.0]-->[-0.9,2.0, 9.1,12, ...]
        supercell_L_Z_max = supercell_Z(L_Z_max, NOC, round(h, 4))  #

        supercell_single_atom_connect = NOC * single_atom_connect

        if charge_r == 1:
            supercell_nonmetal_charge = NOC * nonmetal_charge
            supercell_metal_charge = NOC * metal_charge

        supercell_L_index = supercell_index(L_index, NOC,
                                            nonmetal_length)  # [[1,2,3],[5,7,9]]-->[[1,2,3],[5,7,9], [11,12,13],[15,17,19],...]
        supercell_single_atom_list = supercell_list(single_atom_list, NOC,
                                                    nonmetal_length)  # [1,2,3,4]-->[1,2,3,4,21,22,23,24,....]

        supercell_metal_connect_structure = []
        metal_connect_atoms = []

        for i in metal_connect_structure:
            metal_connect_atoms = metal_connect_atoms + i

        metal_connect_atoms = list(set(metal_connect_atoms))
        supercell_metal_connect_atoms = supercell_list(metal_connect_atoms, NOC, nonmetal_length)
        supercell_metal_connect_structure0 = supercell_index(metal_connect_structure, NOC, nonmetal_length)

        for i in range(0, len(supercell_metal_elemt)):
            each = []
            for k in supercell_metal_connect_structure0[i]:

                sd_k = p_distance_comput(supercell_metal_coa[i], supercell_nonmetal_coa[k], supercell_acod,
                                         supercell_bcod,
                                         supercell_ccod, 1, 1, 1)
                bond_ik = float(atom_r[supercell_metal_elemt[i]]) + float(
                    atom_r[supercell_nonmetal_elemt[k]])
                if sd_k <= bond_ik:
                    each.append(k)

            if len(each) == len(supercell_metal_connect_structure0[i]):
                supercell_metal_connect_structure.append(each)
            else:
                # print("each", each)
                # print("enter else")
                miss = [m for m in supercell_metal_connect_structure0[i] if m not in each]
                # print("miss", miss)
                for n in miss:
                    possible = [m for m in supercell_metal_connect_atoms if (m - n) % len_nonmetal_elemt_u == 0]
                    # print("possible", possible)
                    for j in possible:
                        sd = p_distance_comput(supercell_metal_coa[i], supercell_nonmetal_coa[j], supercell_acod,
                                               supercell_bcod, supercell_ccod, 1, 1, 1)
                        bond_ij = float(atom_r[supercell_metal_elemt[i]]) + float(
                            atom_r[supercell_nonmetal_elemt[j]])  # bond-->bond00
                        if sd <= bond_ij:
                            each.append(j)
                        each = list(set(each))
                        if len(each) == len(supercell_metal_connect_structure0[i]):
                            break

                    # print("final each", each)
                supercell_metal_connect_structure.append(each)

        supercell_metal_connect_code = []
        supercell_atom_linker_atm = []
        supercell_atom_linker_L = []
        for i in range(0, len(supercell_L_index)):  # supercell_L_index[i]= [2,3,5,9]
            for j in supercell_L_index[i]:
                supercell_atom_linker_atm.append(j)
                supercell_atom_linker_L.append(i)

        for i in supercell_metal_connect_structure:  # i=[2,3,5,9]
            supercell_metal_connect_code_each = []
            for j in i:
                k = supercell_atom_linker_atm.index(j)
                # print(k)
                supercell_metal_connect_code_each.append(supercell_atom_linker_L[k])

            supercell_metal_connect_code_each = list(set(supercell_metal_connect_code_each))
            supercell_metal_connect_code.append(supercell_metal_connect_code_each)

        linker_metal_code = linker_metal(supercell_metal_connect_code, len(supercell_L_index))

        ############################
        # define slab_region
        sr_top = 0.5 * (supercell_ccod[2] - 2 * ccod_input[2]) + 0.5 * slab_h  # slab_region_to
        sr_bottom = 0.5 * (supercell_ccod[2] - 2 * ccod_input[2]) - 0.5 * slab_h
        print("slab_region_top", sr_top)
        print("slab_region_bottom", sr_bottom)

        # identify the linkers that preserved in the output model
        linker_out_index = []
        for i in range(0, len(supercell_L_index)):
            min_Z = supercell_L_Z_min[i]
            max_Z = supercell_L_Z_max[i]
            if (min_Z >= sr_bottom and min_Z < sr_top):
                linker_out_index.append(i)  # i is index of linker in the list of supercell_L_index
            elif (max_Z <= sr_top and max_Z > sr_bottom):
                linker_out_index.append(i)
            elif (0.5 * (min_Z + max_Z) <= sr_top and 0.5 * (min_Z + max_Z) > sr_bottom):
                linker_out_index.append(i)
            else:
                pass

        # print("linker_out_index-step1:",linker_out_index)
        #### preserved metals######
        linker_out_atoms = []
        for i in linker_out_index:
            linker_out_atoms = linker_out_atoms + [j for j in supercell_L_index[i]]
        metal_out_index = []
        add_groups_elemt = []
        add_groups_coa = []
        add_groups_charge = []
        H_linker_index = []
        add_index = []
        for i in range(0, len(supercell_metal_elemt)):
            metal_z = supercell_metal_coa[i][2]

            if (metal_z >= sr_bottom and metal_z <= sr_top):
                metal_out_index.append(i)

            else:
                if len(supercell_metal_connect_structure[i]) > 0:
                    metal_connect = [j for j in linker_out_atoms if j in supercell_metal_connect_structure[i]]

                    if (len(metal_connect) / len(supercell_metal_connect_structure[i])) > 0.5 and len(supercell_metal_connect_structure[i])>1:  # 接近配位饱和的金属保留

                        metal_out_index.append(i)

                        print("metal_connect",metal_connect)
                        # print("supercell_metal_connect_structure[i]", supercell_metal_connect_structure[i])

                    elif (len(metal_connect) / len(supercell_metal_connect_structure[i])) <= 0.5 and len(
                            metal_connect) >= 3:  #####3 coordinated metals are also preserved in the output

                        metal_out_index.append(i)
                        # print("metal_connect",metal_connect)
                        # print("supercell_metal_connect_structure[i]", supercell_metal_connect_structure[i])

                    elif (len(metal_connect) / len(supercell_metal_connect_structure[i])) <= 0.5 and len(
                            metal_connect) < 3:  # less coordinated metals are not preserved

                        add_index = add_index + [k for k in metal_connect]
                        # print("add_index",add_index)
                    elif (len(metal_connect) / len(supercell_metal_connect_structure[i])) > 0.5 and len(supercell_metal_connect_structure[i]) == 1:
                        add_index = add_index + [k for k in metal_connect]

        if protonation_yn == 1:
            add_index = list(set(add_index))
            # print("add_index",add_index)
            add_H_n = []
            i01 = [0, 1]
            for l in add_index:
                add_yn = 1
                add_H_n.append(l)
                #######the following defines several function groups for the surface saturation
                if supercell_nonmetal_elemt[l] == "O":
                    # print("l", l)
                    k = supercell_atom_linker_atm.index(l)

                    for m in supercell_L_index[supercell_atom_linker_L[k]]:
                        if m != l and supercell_nonmetal_elemt[m] == "C":
                            dm = p_distance_comput(supercell_nonmetal_coa[m], supercell_nonmetal_coa[l], supercell_acod,
                                                   supercell_bcod, supercell_ccod, 1, 1, 1)
                            if dm <= float(atom_r[supercell_nonmetal_elemt[m]]) + float(
                                    atom_r[supercell_nonmetal_elemt[l]]):
                                # means C-O bond
                                C_connect = []
                                for n in supercell_L_index[supercell_atom_linker_L[k]]:
                                    if n != l and n != m:
                                        dn = p_distance_comput(supercell_nonmetal_coa[m],
                                                               supercell_nonmetal_coa[n], supercell_acod,
                                                               supercell_bcod, supercell_ccod, 1, 1, 1)
                                        if dn <= float(atom_r[supercell_nonmetal_elemt[m]]) + float(
                                                atom_r[supercell_nonmetal_elemt[n]]):
                                            C_connect.append(n)
                                if len(C_connect) == 3:
                                    add_yn = 1
                                elif len(C_connect) == 1:
                                    add_yn = 0
                                elif len(C_connect) == 2:
                                    i0, i1 = i01[0], i01[1]
                                    if supercell_nonmetal_elemt[C_connect[i0]] == "O":  ######-COO detected for -COOH
                                        if C_connect[i0] in add_index:
                                            add_index.remove(C_connect[i0])
                                            if i01 == [0, 1]:
                                                i01 = [1, 0]
                                            elif i01 == [1, 0]:
                                                i01 = [0, 1]

                                        add_yn = 1

                                    elif supercell_nonmetal_elemt[C_connect[i1]] == "O":
                                        if C_connect[i1] in add_index:
                                            add_index.remove(C_connect[i1])
                                            if i01 == [0, 1]:
                                                i01 = [1, 0]
                                            elif i01 == [1, 0]:
                                                i01 = [0, 1]
                                        add_yn = 1
                                    else:
                                        add_yn = 0  #

                        if m != l and supercell_nonmetal_elemt[m] == "N":  # -NO2,-NO
                            dm = p_distance_comput(supercell_nonmetal_coa[m], supercell_nonmetal_coa[l], supercell_acod,
                                                   supercell_bcod, supercell_ccod, 1, 1, 1)
                            if dm <= float(atom_r[supercell_nonmetal_elemt[m]]) + float(
                                    atom_r[supercell_nonmetal_elemt[l]]):
                                # means N-O bond,
                                add_yn = 0
                                break

                        if m != l and supercell_nonmetal_elemt[m] == "S":
                            dm = p_distance_comput(supercell_nonmetal_coa[m], supercell_nonmetal_coa[l], supercell_acod,
                                                   supercell_bcod, supercell_ccod, 1, 1, 1)
                            if dm <= float(atom_r[supercell_nonmetal_elemt[m]]) + float(
                                    atom_r[supercell_nonmetal_elemt[l]]):
                                # means C-O bond
                                S_connect = []
                                for n in supercell_L_index[supercell_atom_linker_L[k]]:
                                    if n != l and n != m:
                                        dn = p_distance_comput(supercell_nonmetal_coa[m],
                                                               supercell_nonmetal_coa[n], supercell_acod,
                                                               supercell_bcod, supercell_ccod, 1, 1, 1)
                                        if dn <= float(atom_r[supercell_nonmetal_elemt[m]]) + float(
                                                atom_r[supercell_nonmetal_elemt[n]]):
                                            S_connect.append(n)

                                if len(S_connect) == 2:  # R-SO-R'
                                    add_yn = 0

                                elif len(S_connect) == 3:  # -SO3H -SO2-
                                    S_O = [x for x in S_connect if supercell_nonmetal_elemt[x] == "O"]
                                    if len(S_O) == 1:
                                        add_yn = 0
                                    elif len(S_O) == 2:
                                        add_yn = 1
                                        i0, i1 = i01[0], i01[1]
                                        if S_O[i0] in add_index:
                                            add_index.remove(S_O[i0])
                                            if i01 == [0, 1]:
                                                i01 = [1, 0]
                                            elif i01 == [1, 0]:
                                                i01 = [0, 1]
                                        elif S_O[i1] in add_index:
                                            add_index.remove(S_O[i1])
                                            if i01 == [0, 1]:
                                                i01 = [1, 0]
                                            elif i01 == [1, 0]:
                                                i01 = [0, 1]

                        if m != l and supercell_nonmetal_elemt[m] == "P":
                            dm = p_distance_comput(supercell_nonmetal_coa[m], supercell_nonmetal_coa[l], supercell_acod,
                                                   supercell_bcod, supercell_ccod, 1, 1, 1)
                            if dm <= float(atom_r[supercell_nonmetal_elemt[m]]) + float(
                                    atom_r[supercell_nonmetal_elemt[l]]):

                                P_connect = []
                                for n in supercell_L_index[supercell_atom_linker_L[k]]:
                                    if n != l and n != m:
                                        dn = p_distance_comput(supercell_nonmetal_coa[m],
                                                               supercell_nonmetal_coa[n], supercell_acod,
                                                               supercell_bcod, supercell_ccod, 1, 1, 1)
                                        if dn <= float(atom_r[supercell_nonmetal_elemt[m]]) + float(
                                                atom_r[supercell_nonmetal_elemt[n]]):
                                            P_connect.append(n)
                                P_O = [x for x in P_connect if supercell_nonmetal_elemt[x] == "O"]
                                if len(P_O) == 2:
                                    add_yn = 1
                                    if P_O[0] in add_index and P_O[1] in add_index:
                                        add_index.remove(P_O[0])
                                    else:
                                        pass
                                if len(P_O) == 3:
                                    add_yn = 1
                                    if P_O[0] in add_index and P_O[1] in add_index:
                                        add_index.remove(P_O[0])
                                    elif P_O[0] in add_index and P_O[2] in add_index:
                                        add_index.remove(P_O[0])
                                    elif P_O[1] in add_index and P_O[2] in add_index:
                                        add_index.remove(P_O[1])
                                    else:
                                        pass


                elif supercell_nonmetal_elemt[l] == "C":
                    k = supercell_atom_linker_atm.index(l)
                    C_connect = []
                    for m in supercell_L_index[supercell_atom_linker_L[k]]:
                        if m != l:
                            dm = p_distance_comput(supercell_nonmetal_coa[m], supercell_nonmetal_coa[l],
                                                   supercell_acod, supercell_bcod, supercell_ccod, 1, 1, 1)
                            if dm <= float(atom_r[supercell_nonmetal_elemt[m]]) + float(
                                    atom_r[supercell_nonmetal_elemt[l]]):
                                # means N-X bond
                                C_connect.append(m)
                    if len(C_connect) > 1:
                        add_yn = 0
                    else:
                        pass


                elif supercell_nonmetal_elemt[l] == "N":
                    k = supercell_atom_linker_atm.index(l)
                    N_connect = []
                    for m in supercell_L_index[supercell_atom_linker_L[k]]:
                        if m != l:
                            dm = p_distance_comput(supercell_nonmetal_coa[m], supercell_nonmetal_coa[l],
                                                   supercell_acod, supercell_bcod, supercell_ccod, 1, 1, 1)
                            if dm <= float(atom_r[supercell_nonmetal_elemt[m]]) + float(
                                    atom_r[supercell_nonmetal_elemt[l]]):
                                # means N-X bond
                                N_connect.append(m)
                    if len(N_connect) == 3:
                        add_yn = 0

                    elif len(N_connect) == 1:
                        if supercell_nonmetal_elemt[N_connect[0]] == "C":
                            add_yn = 0

                    elif len(N_connect) == 2:
                        if supercell_nonmetal_elemt[N_connect[0]] == "C" and supercell_nonmetal_elemt[
                            N_connect[1]] == "C":
                            NC0_0 = supercell_nonmetal_coa[N_connect[0]][0] - supercell_nonmetal_coa[l][0]
                            NC0_1 = supercell_nonmetal_coa[N_connect[0]][1] - supercell_nonmetal_coa[l][1]
                            NC0_2 = supercell_nonmetal_coa[N_connect[0]][2] - supercell_nonmetal_coa[l][2]
                            NC1_0 = supercell_nonmetal_coa[N_connect[1]][0] - supercell_nonmetal_coa[l][0]
                            NC1_1 = supercell_nonmetal_coa[N_connect[1]][1] - supercell_nonmetal_coa[l][1]
                            NC1_2 = supercell_nonmetal_coa[N_connect[1]][2] - supercell_nonmetal_coa[l][2]
                            NC0 = [NC0_0, NC0_1, NC0_2]
                            NC1 = [NC1_0, NC1_1, NC1_2]

                            NC0_len = math.sqrt(NC0[0] * NC0[0] + NC0[1] * NC0[1] + NC0[2] * NC0[2])
                            NC1_len = math.sqrt(NC1[0] * NC1[0] + NC1[1] * NC1[1] + NC1[2] * NC1[2])

                            cos_CNC = (NC0[0] * NC1[0] + NC0[1] * NC1[1] + NC0[2] * NC1[2]) / (
                                    NC0_len * NC1_len)

                            CNC = round(180 * math.acos(cos_CNC) / 3.14159265, 6)
                            if abs(CNC - 120) <= 6:
                                add_yn = 0
                            elif abs(CNC - 108) < 6:
                                for n in supercell_L_index[supercell_atom_linker_L[k]]:
                                    if n != l and n != N_connect[0] and n != N_connect[1] and supercell_nonmetal_elemt[
                                        n] == "N":
                                        dn0 = p_distance_comput(supercell_nonmetal_coa[N_connect[0]],
                                                                supercell_nonmetal_coa[n], supercell_acod,
                                                                supercell_bcod, supercell_ccod, 1, 1, 1)
                                        dn1 = p_distance_comput(supercell_nonmetal_coa[N_connect[1]],
                                                                supercell_nonmetal_coa[n], supercell_acod,
                                                                supercell_bcod, supercell_ccod, 1, 1, 1)
                                        C_N_d = float(atom_r[supercell_nonmetal_elemt[N_connect[0]]]) + float(
                                            atom_r[supercell_nonmetal_elemt[n]])
                                        if dn0 <= C_N_d or dn1 <= C_N_d:  # 表示发现C-N-C-N结构
                                            if n in add_index:
                                                add_yn = 1
                                                add_index.remove(n)
                                                break
                                            else:
                                                break_yn = 0
                                                for o in supercell_L_index[supercell_atom_linker_L[k]]:
                                                    if o != l and o != N_connect[0] and o != N_connect[1] and o != n and  supercell_nonmetal_elemt[o] == "H":
                                                        dh0 = p_distance_comput(supercell_nonmetal_coa[o],
                                                                                supercell_nonmetal_coa[n],
                                                                                supercell_acod,
                                                                                supercell_bcod, supercell_ccod, 1, 1, 1)
                                                        if dh0 <= float(atom_r[supercell_nonmetal_elemt[n]]) + float(
                                                                atom_r[supercell_nonmetal_elemt[o]]):
                                                            add_yn = 0  # C-N-NH-C in 5-membe-ring
                                                            break_yn = 1
                                                            break
                                                        else:
                                                            pass

                                                if break_yn == 1:
                                                    break
                                                elif break_yn == 0:
                                                    add_yn = 1
                            else:
                                add_yn = 1
                        elif (supercell_nonmetal_elemt[N_connect[0]] == "C" and supercell_nonmetal_elemt[
                            N_connect[1]] == "N"):
                            if N_connect[1] in add_index:
                                add_yn = 1
                                add_index.remove(N_connect[1])
                                print("findN-N1")

                            else:
                                break_yn = 0
                                for o in supercell_L_index[supercell_atom_linker_L[k]]:
                                    if o != l and o != N_connect[0] and o != N_connect[1] and supercell_nonmetal_elemt[o] == "H":
                                        dh0 = p_distance_comput(supercell_nonmetal_coa[o],
                                                                supercell_nonmetal_coa[N_connect[1]],
                                                                supercell_acod,
                                                                supercell_bcod, supercell_ccod, 1, 1, 1)
                                        if dh0 <= float(atom_r[supercell_nonmetal_elemt[N_connect[1]]]) + float(
                                                atom_r[supercell_nonmetal_elemt[o]]):
                                            add_yn = 0
                                            break_yn = 1
                                            break
                                        else:
                                            pass

                                if break_yn == 1:
                                    break
                                elif break_yn == 0:  #
                                    add_yn = 1
                        elif (supercell_nonmetal_elemt[N_connect[1]] == "C" and supercell_nonmetal_elemt[
                            N_connect[0]] == "N"):
                            if N_connect[0] in add_index:
                                add_yn = 1
                                add_index.remove(N_connect[0])
                                #print("findN-N2")

                            else:
                                break_yn = 0
                                for o in supercell_L_index[supercell_atom_linker_L[k]]:
                                    if o != l and o != N_connect[0] and o != N_connect[1] and \
                                            supercell_nonmetal_elemt[o] == "H":
                                        dh0 = p_distance_comput(supercell_nonmetal_coa[o],
                                                                supercell_nonmetal_coa[N_connect[0]],
                                                                supercell_acod,
                                                                supercell_bcod, supercell_ccod, 1, 1, 1)
                                        if dh0 <= float(atom_r[supercell_nonmetal_elemt[N_connect[0]]]) + float(
                                                atom_r[supercell_nonmetal_elemt[o]]):
                                            add_yn = 0
                                            break_yn = 1
                                            break
                                        else:
                                            pass

                                if break_yn == 1:
                                    break
                                elif break_yn == 0:
                                    add_yn = 1
                if add_yn == 1:
                    add_groups_elemt.append("H")
                    k = supercell_atom_linker_atm.index(l)
                    H_linker_index.append(supercell_atom_linker_L[k])
                    i = 0
                    for mi in range(0, len(supercell_metal_connect_structure)):
                        if l in supercell_metal_connect_structure[mi] and mi not in metal_out_index:
                            i = mi
                            break
                    metal_Linker_d = distance_comput(supercell_metal_coa[i], supercell_nonmetal_coa[
                        l])
                    bond_value = float(atom_r[supercell_nonmetal_elemt[l]]) + float(
                        atom_r[supercell_metal_elemt[i]])
                    # print(bond_value)
                    if metal_Linker_d <= bond_value:
                        H0 = 0.96 * (supercell_metal_coa[i][0] - supercell_nonmetal_coa[l][0]) / metal_Linker_d + \
                             supercell_nonmetal_coa[l][0]  # the length is 0.96 A
                        H1 = 0.96 * (supercell_metal_coa[i][1] - supercell_nonmetal_coa[l][1]) / metal_Linker_d + \
                             supercell_nonmetal_coa[l][1]
                        H2 = 0.96 * (supercell_metal_coa[i][2] - supercell_nonmetal_coa[l][2]) / metal_Linker_d + \
                             supercell_nonmetal_coa[l][2]
                        add_groups_coa.append([H0, H1, H2])
                        if abs(H2 - 67.869376) < 0.001 and abs(H1 - 96.9915) < 0.01 and abs(H0 - 275.252) < 0.01:
                            print("l1,mi1", l, mi)
                            print(supercell_nonmetal_coa[l], supercell_metal_coa[i])

                        if abs(H2 - 67.0883) < 0.001 and abs(H1 - 97.7266) < 0.01 and abs(H0 - 275.342) < 0.01:
                            print("l2,mi2", l, mi)
                            print(supercell_nonmetal_coa[l], supercell_metal_coa[i])

                        add_H_n.append(mi)

                        if charge_r == 1:
                            if supercell_metal_charge[i]!=0:
                                H_charge = H_charge_set  #H_charge = supercell_metal_charge[i] / len(supercell_metal_connect_structure[i])
                            elif supercell_metal_charge[i]==0:
                                H_charge=0
                            add_groups_charge.append(str(round(H_charge, 2)))

                    elif metal_Linker_d > bond_value:
                        m_rcoai = coa_conv2_rcoa(supercell_metal_coa[i], supercell_acod, supercell_bcod,
                                                 supercell_ccod)
                        L_rcoal = coa_conv2_rcoa(supercell_nonmetal_coa[l], supercell_acod, supercell_bcod,
                                                 supercell_ccod)
                        m_rcoai0, m_rcoai1, m_rcoai2 = 0, 0, 0
                        if abs(m_rcoai[0] - L_rcoal[0]) > 0.5:
                            if L_rcoal[0] >= 0.5:
                                m_rcoai0 = m_rcoai[0] + 1
                            elif L_rcoal[0] < 0.5:
                                m_rcoai0 = m_rcoai[0] - 1
                        else:
                            m_rcoai0 = m_rcoai[0]

                        if abs(m_rcoai[1] - L_rcoal[1]) > 0.5:
                            if L_rcoal[1] >= 0.5:
                                m_rcoai1 = m_rcoai[1] + 1
                            elif L_rcoal[1] < 0.5:
                                m_rcoai1 = m_rcoai[1] - 1
                        else:
                            m_rcoai1 = m_rcoai[1]

                        m_rcoai2 = m_rcoai[2]

                        coa_m = rcoa_conv2_coa([m_rcoai0, m_rcoai1, m_rcoai2], supercell_acod, supercell_bcod,
                                               supercell_ccod)

                        metal_Linker_d2 = distance_comput(supercell_nonmetal_coa[l], coa_m)
                        H0 = 0.96 * (coa_m[0] - supercell_nonmetal_coa[l][
                            0]) / metal_Linker_d2 + supercell_nonmetal_coa[l][0]
                        H1 = 0.96 * (coa_m[1] - supercell_nonmetal_coa[l][
                            1]) / metal_Linker_d2 + supercell_nonmetal_coa[l][1]
                        H2 = 0.96 * (coa_m[2] - supercell_nonmetal_coa[l][
                            2]) / metal_Linker_d2 + supercell_nonmetal_coa[l][2]
                        if abs(H2 - 67.869376) < 0.001 and abs(H1 - 96.9915) < 0.01 and abs(H0 - 275.252) < 0.01:
                            print("2,l1,mi1", l, mi)

                        if abs(H2 - 67.0883) < 0.001 and abs(H1 - 97.7266) < 0.01 and abs(H0 - 275.342) < 0.01:
                            print("2,l2,mi2", l, mi)
                        add_H_n.append(mi)
                        add_groups_coa.append([H0, H1, H2])
                        if charge_r == 1:
                            if supercell_metal_charge[i]!=0:
                                H_charge = H_charge_set  #H_charge = supercell_metal_charge[i] / len(supercell_metal_connect_structure[i])
                            elif supercell_metal_charge[i]==0:
                                H_charge=0
                            add_groups_charge.append(str(round(H_charge, 3)))

        # print("add_H_n",add_H_n)
        # the linkers connected to the preserved metals are also preserved
        linker_out_final = []
        for i in linker_out_index:  # i in [2,4,6,8,10]
            if len(linker_metal_code[i]) == 0:
                linker_out_final.append(i)
            if len(linker_metal_code[i]) > 0:
                for j in linker_metal_code[i]:
                    # linker_metal_code[i]=[1, 11, 15] means linker i connected to metal 1，11，15
                    if j in metal_out_index:  # metal_out_index=[2,5,11,7,8,...]
                        linker_out_final.append(i)
                        break
                    else:
                        pass

        # the preserved single nonmetalic atoms
        single_atom_out = []
        for single_index in supercell_single_atom_list:
            ndx=supercell_single_atom_list.index(single_index)
            connect_1 = supercell_single_atom_connect[ndx]
            connect_2 = 0
            for ii in range(0, len(metal_out_index)):
                i = metal_out_index[ii]

                d2 = p_distance_comput(supercell_metal_coa[i], supercell_nonmetal_coa[single_index], supercell_acod, supercell_bcod,
                                                   supercell_ccod, 1, 1, 1,)
                m_bond2 = float(atom_r[supercell_metal_elemt[i]]) + float(
                    atom_r[supercell_nonmetal_elemt[single_index]])

                if d2 <= m_bond2:
                    connect_2 = connect_2+1

            if connect_2==1 and connect_1==1:
                single_atom_out.append(single_index)
            elif connect_2>=2:
                single_atom_out.append(single_index)
                    # print("supercell_metal_coa[i]", supercell_metal_coa[i])
                    # print("supercell_nonmetal_coa[single_index]", supercell_nonmetal_coa[single_index])
                    # print("m_d02, m_d12", m_d02, m_d12)

        alpha = 90.0000
        beta = 90.0000



        linker_out_atoms2 = []
        for i in linker_out_final:
            linker_out_atoms2 = linker_out_atoms2 + [j for j in supercell_L_index[i]]

        ###coa_in_grid####
        out_coa=[]
        out_rcoa=[]
        out_elemt=[]
        out_charge=[]
        out_z =[]

        for i in linker_out_atoms2:
            c_final = coa_in_grid(supercell_nonmetal_coa[i],supercell_acod,supercell_bcod,[0,0,slab_h + vacuum])
            coa_final=c_final[0]
            rcoa_final=c_final[1]
            out_coa.append(coa_final)
            out_z.append(coa_final[2])
            out_rcoa.append(rcoa_final)
            out_elemt.append(supercell_nonmetal_elemt[i])
            if charge_r == 1:
                out_charge.append(supercell_nonmetal_charge[i])
            elif charge_r == 0:
                out_charge.append("0.000")
        for i in single_atom_out:
            c_final = coa_in_grid(supercell_nonmetal_coa[i],supercell_acod,supercell_bcod,[0,0,slab_h + vacuum])
            coa_final=c_final[0]
            rcoa_final=c_final[1]
            out_coa.append(coa_final)
            out_z.append(coa_final[2])
            out_rcoa.append(rcoa_final)
            out_elemt.append(supercell_nonmetal_elemt[i])
            if charge_r == 1:
                out_charge.append(supercell_nonmetal_charge[i])
            elif charge_r == 0:
                out_charge.append("0.000")

        for i in metal_out_index:
            c_final = coa_in_grid(supercell_metal_coa[i],supercell_acod,supercell_bcod,[0,0,slab_h + vacuum])
            coa_final=c_final[0]
            rcoa_final=c_final[1]

            out_coa.append(coa_final)
            out_z.append(coa_final[2])
            out_rcoa.append(rcoa_final)

            out_elemt.append(supercell_metal_elemt[i])
            if charge_r == 1:
                out_charge.append(supercell_metal_charge[i])
            elif charge_r == 0:
                out_charge.append("0.000")


        for i in range(0, len(add_groups_elemt)):
            if H_linker_index[i] in linker_out_final:
                c_final = coa_in_grid(add_groups_coa[i], supercell_acod, supercell_bcod,[0, 0, slab_h + vacuum])
                coa_final = c_final[0]
                rcoa_final = c_final[1]
                out_coa.append(coa_final)
                out_z.append(coa_final[2])
                out_rcoa.append(rcoa_final)

                out_elemt.append(add_groups_elemt[i])
                if charge_r == 1:
                    out_charge.append(add_groups_charge[i])
                elif charge_r == 0:
                    out_charge.append("0.000")


        if min_X_set==1:
            Nx_out=math.ceil(min_X_length/supercell_acod[0])
            out_coa=supercell_coa(out_coa,Nx_out,supercell_acod)
            out_elemt=Nx_out*out_elemt
            out_charge=Nx_out*out_charge
            alen_input=Nx_out*alen_input
            supercell_acod2=[Nx_out*supercell_acod[0],Nx_out*supercell_acod[1],Nx_out*supercell_acod[2]]
        else:
            supercell_acod2=supercell_acod

        if min_Y_set==1:
            Ny_out=math.ceil(min_Y_length/supercell_bcod[1])
            out_coa=supercell_coa(out_coa,Ny_out,supercell_bcod)
            out_elemt=Ny_out*out_elemt
            out_charge=Ny_out*out_charge
            blen_input=Ny_out*blen_input
            supercell_bcod2=[Ny_out*supercell_bcod[0],Ny_out*supercell_bcod[1],Ny_out*supercell_bcod[2]]
        else:
            supercell_bcod2=supercell_bcod

#############add graphene layer############
        sl_coa_output = []
        if separate_layer == 1:
            with open(".\\Barrier_layer\\barrier_layer.car",encoding="utf-8") as sl_f:  #users can choose other nanosheets
                sl_lines = sl_f.readlines()  # car_lines-->sl_lines
                #              ####
                sl_z = []
                for l_number in range(0, len(sl_lines)):
                    if len(re.findall(r"PBC", sl_lines[l_number])) > 0 and len(
                            re.findall(r"[0-9]", sl_lines[l_number])) > 0:
                        sl_lines_line_t = re.sub("\t", " ", sl_lines[l_number])
                        sl_lines_line_t = sl_lines_line_t.strip('\n') + "  "
                        sl_lattice = re.findall(
                            r'([+-]?(?:[0-9]\d*)(?:\.\d*)? |[+-]?(?:[0-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+))',
                            sl_lines_line_t)
                        #print(lattice)
                        sl_alen = float(sl_lattice[0])
                        sl_blen = float(sl_lattice[1])
                        sl_clen = float(sl_lattice[2])
                        sl_aang = float(sl_lattice[3])
                        sl_bang = float(sl_lattice[4])
                        sl_cang = float(sl_lattice[5])
                        print("lattice:", sl_alen, sl_blen, sl_clen, sl_aang, sl_bang, sl_cang)

                sl_aa = 3.14159265 * sl_aang / 180
                sl_ba = 3.14159265 * sl_bang / 180
                # print(pbc[1])
                sl_ca = 3.14159265 * sl_cang / 180
                sl_A = math.cos(sl_aa)
                sl_B = math.cos(sl_ba)
                sl_C = math.cos(sl_ca)

                sl_h0 = sl_clen * math.sqrt(
                    1 - sl_A * sl_A - sl_B * sl_B - sl_C * sl_C + 2 * sl_A * sl_B * sl_C) / math.sin(sl_ca)

                sl_acod = [sl_alen, 0, 0]
                sl_bcod = [sl_blen * sl_C, sl_blen * math.sin(sl_ca), 0]
                sl_ccodx = sl_clen * sl_B
                sl_ccodz = sl_h0
                sl_ccody = (sl_A * sl_blen * sl_clen - sl_bcod[0] * sl_ccodx - sl_bcod[2] * sl_ccodz) / sl_bcod[1]
                sl_ccod = [sl_ccodx, sl_ccody, sl_ccodz]
                print("sl_ccod=", sl_ccod)

                sl_structure = sl_lines[5:]  #
                N_of_end = sl_structure.count("end\n")
                for jN in range(0, N_of_end):
                    sl_structure.remove("end\n")
                # print(N_of_end)

                N_of_n = sl_structure.count("\n")
                for iN in range(0, N_of_n):
                    sl_structure.remove("\n")

                # print("*****")
                # print(structure)
                sl_elemt_all = []  #[H,C,O,H,C,N,N]
                sl_coa_all1 = []  # coordinate of atom
                sl_coa_all = []
                sl_charge_all = []

                for line in sl_structure:

                    # print(L[0:9])
                    # print(elemt_single)
                    sl_line_t = re.sub("\t", " ", line)
                    sl_line_t = sl_line_t.strip('\n') + "  "
                    # print("car lines:",line_t)
                    sl_coa_charge = re.findall(
                        r'( [+-]?(?:[0-9]\d*)(?:\.\d*)? |[+-]?(?:[0-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+))', sl_line_t)
                    L = sl_line_t.split()
                    elemt_single = L[0:9][7]

                    sl_coordinate = sl_coa_charge[0:3]
                    if charge_r == 1 and len(sl_coa_charge)>4:
                        sl_charge=sl_coa_charge[4]

                    else:
                        sl_charge="0.0"

                    # print("coa_charge_car:",coa_charge)

                    coa_single = [0, 0, 0]
                    coa_single[0] = float(sl_coordinate[0])
                    coa_single[1] = float(sl_coordinate[1])
                    coa_single[2] = float(sl_coordinate[2])
                    sl_z.append(coa_single[2])
                    # print("coa_single:",coa_single)

                    sl_elemt_all.append(elemt_single)
                    sl_coa_all1.append(coa_single)
                    sl_charge_all.append(sl_charge)

            print(sl_elemt_all)
            print(sl_coa_all1)
            total_z = (min(sl_z) + max(sl_z)) / 2
            print(total_z)

            for i in sl_coa_all1:
                sl_coa_all.append([i[0], i[1], i[2] - total_z])

            print(sl_coa_all)

            pt1 = [0, 0, 0]
            pt2 = supercell_acod2
            print("pt2,acod_input",pt2,acod_input)
            pt3 = supercell_bcod2
            print("pt3,bcod_input", pt3, bcod_input)
            pt4 = [pt2[0] + pt3[0], pt2[1] + pt3[1], pt2[2] + pt3[2]]

            if cang_input <= sl_cang:
                pt4_rcoa = coa_conv2_rcoa(pt4, sl_acod, sl_bcod, [0, 0, 1])
                Nx = math.ceil(pt4_rcoa[0])
                print("Nx",Nx)
                Ny = math.ceil(pt3[1] / sl_bcod[1])
                supercell_sl_coa = supercell_coa(supercell_coa(sl_coa_all, Nx, sl_acod), Ny, sl_bcod)
                supercell_sl_elemt_all = Nx * Ny * sl_elemt_all
                supercell_sl_charge_all = Nx * Ny * sl_charge_all

            elif cang_input > sl_cang:
                Nx1 = math.ceil(pt2[0] / sl_alen)
                print("Nx1",Nx1)
                Ny = math.ceil(pt3[1] / sl_bcod[1])
                supercell_sl_coa = supercell_coa(supercell_coa(sl_coa_all, Nx1, sl_acod), Ny, sl_bcod)

                pt3_rcoa = coa_conv2_rcoa(pt3, sl_acod, sl_bcod, [0, 0, 1])
                Nx2 = abs(math.ceil(pt3_rcoa[0]))
                print("Nx2",Nx2)
                sl_acod2 = [-1 * i for i in sl_acod]
                supercell_sl_coa = supercell_sl_coa + supercell_coa(supercell_coa(sl_coa_all, Nx2, sl_acod2), Ny,
                                                                    sl_bcod)
                supercell_sl_elemt_all = Nx1 * Ny * sl_elemt_all + Nx2 * Ny * sl_elemt_all
                supercell_sl_charge_all = Nx1 * Ny * sl_charge_all + Nx2 * Ny * sl_charge_all

            sl_elemt_output = []
            sl_charge_output = []

            pt_supercell = [pt1, pt2, pt4, pt3]
            out_polygon = geometry.Polygon(pt_supercell)
            for i in range(0, len(supercell_sl_coa)):
                point = geometry.Point(supercell_sl_coa[i][0:2])
                out_line1 = geometry.LineString([pt1, pt2])
                out_line2 = geometry.LineString([pt2, pt4])
                out_line3 = geometry.LineString([pt4, pt3])
                out_line4 = geometry.LineString([pt3, pt1])
                # if len(grid_ij[i][j]) > 0:
                line_distance = min(out_line2.distance(point), out_line3.distance(point))
                edge1_distance=out_line1.distance(point)
                edge4_distance = out_line4.distance(point)

                if out_polygon.covers(point) == True and line_distance > 0.83: #usually 0.83 is enough to prevent leakage of the guest
                    sl_coa_output.append(supercell_sl_coa[i])
                    sl_elemt_output.append(supercell_sl_elemt_all[i])
                    sl_charge_output.append(supercell_sl_charge_all[i])

            for i in range(0, len(sl_coa_output)):
                c_final = coa_in_grid(sl_coa_output[i], supercell_acod2, supercell_bcod2,[0, 0, slab_h + vacuum])
                coa_final = c_final[0]
                rcoa_final = c_final[1]

                out_coa.append(coa_final)
                out_rcoa.append(rcoa_final)
                out_elemt.append(sl_elemt_output[i])

                if charge_r == 1:
                    out_charge.append(sl_charge_output[i])

                elif charge_r == 0:
                    out_charge.append(0.00)


        charge_sum=0
        for i in out_charge:
            charge_sum = charge_sum+float(i)
        log_lines.append("X Length (along a-axis): " + str(round(alen_input, 4)) + "\n")
        #print("blen_input,cang_input",blen_input,cang_input)
        yl=round(blen_input*math.sin(3.14159265 * cang_input / 180), 4)
        #print("blen_input,cang_input,yl", blen_input, cang_input,yl)
        log_lines.append("Y Length (perpendicular to a-axis): " + str(yl) + "\n")
        log_lines.append("Summary of Charge: " + str(charge_sum)+"\n################\n\n")


        #print("len(sl_coa_output)",len(sl_coa_output))
        if "car" in output_files_type:
            out_lines = []
            out_lines.append("!BIOSYM archive 3")
            out_lines.append("\nPBC=ON")
            out_lines.append("\nGenerated by MOF-Membrane Constructor")
            out_lines.append("\n!DATE")

            out_lines.append("\nPBC   " + str(round(alen_input, 4)) + "   " + str(round(blen_input, 4)) + "   " + str(
                round(slab_h + vacuum, 4)) + "   " + str(round(alpha, 4)) + "   " + str(round(beta, 4)) + "   " + str(
                round(cang_input, 4)) + " (P1)")

            for i in range(0,len(out_coa)-len(sl_coa_output)):
                out_lines.append(
                    "\n" + '{0:<5}'.format(out_elemt[i] + "1"))
                out_lines.append('{0:>15}'.format(str(round(out_coa[i][0], 9))) +
                                 '{0:>15}'.format(str(round(out_coa[i][1], 9))) +
                                 '{0:>15}'.format(str(round(out_coa[i][2], 9))) + " MOF  1      ")
                out_lines.append(
                    '{0:<8}'.format(out_elemt[i]) + '{0:<4}'.format(out_elemt[i]))
                if charge_r == 1:
                    out_lines.append("  " + '{0:>6}'.format(str(out_charge[i])))
                if charge_r == 0:
                    out_lines.append("  0.000")
                # out_lines.append("   " + str(i))
            for i in range(len(out_coa)-len(sl_coa_output),len(out_coa)):
                out_lines.append(
                    "\n" + '{0:<5}'.format(out_elemt[i] + "1"))
                out_lines.append('{0:>15}'.format(str(round(out_coa[i][0], 9))) +
                                 '{0:>15}'.format(str(round(out_coa[i][1], 9))) +
                                 '{0:>15}'.format(str(round(out_coa[i][2], 9))) + " GRA  1      ")
                out_lines.append(
                    '{0:<8}'.format(out_elemt[i]) + '{0:<4}'.format(out_elemt[i]))
                if charge_r == 1:
                    out_lines.append("  " + '{0:>6}'.format(str(out_charge[i])))
                if charge_r == 0:
                    out_lines.append("  0.000")

            out_lines.append("\nend\nend\n ")
            file = open(outpath + Filename + "_" + facet_index + ".car", 'w')
            for ij in out_lines:
                file.write(str(ij))
            file.close()
##############pdb##########################
        if "pdb" in output_files_type:
            pdb_lines = []
            pdb_lines.append("REMARK   Generated by MOF-Membrane Constructor")
            pdb_lines.append("\nCRYST1   " + str(round(alen_input, 3)) + "   " + str(round(blen_input, 3)) + "   " + str(
                round(slab_h + vacuum, 3)) + "   " + str(round(alpha, 2)) + "   " + str(round(beta, 2)) + "   " + str(
                round(cang_input, 2)) + " P 1           1  ")

            atom_n=1

            for i in range(0,len(out_coa)-len(sl_coa_output)):
                pdb_lines.append("\n" + '{0:<5}'.format("ATOM"))
                pdb_lines.append('{0:>6}'.format(str(atom_n)))
                atom_n=atom_n+1
                pdb_lines.append('{0:>3}'.format(out_elemt[i]))
                pdb_lines.append("   MOF X   1    " + '{0:>8}'.format(str(round(out_coa[i][0], 3))) +
                                 '{0:>8}'.format(str(round(out_coa[i][1], 3))) +
                                 '{0:>8}'.format(str(round(out_coa[i][2], 3))) + "  1.00  0.00          ")
                pdb_lines.append('{0:>2}'.format(out_elemt[i]))
                if charge_r == 1:
                    pdb_lines.append(" ")
                    pdb_lines.append(str(out_charge[i]))

                elif charge_r == 0:
                    pdb_lines.append(" ")
                # out_lines.append("   " + str(i))

            for i in range(len(out_coa)-len(sl_coa_output),len(out_coa)):
                pdb_lines.append("\n" + '{0:<5}'.format("ATOM"))
                pdb_lines.append('{0:>6}'.format(str(atom_n)))
                atom_n=atom_n+1
                pdb_lines.append('{0:>3}'.format(out_elemt[i]))
                pdb_lines.append("   GRA X   1    " + '{0:>8}'.format(str(round(out_coa[i][0], 3))) +
                                 '{0:>8}'.format(str(round(out_coa[i][1], 3))) +
                                 '{0:>8}'.format(str(round(out_coa[i][2], 3))) + "  1.00  0.00          ")
                pdb_lines.append('{0:>2}'.format(out_elemt[i]))
                if charge_r == 1:
                    pdb_lines.append(" ")
                    pdb_lines.append(str(out_charge[i]))

                elif charge_r == 0:
                    pdb_lines.append(" ")
                # out_lines.append("   " + str(i))

            pdb_lines.append("\nEND\n ")
            file = open(outpath + Filename + "_" + facet_index + ".pdb", 'w')
            for ij in pdb_lines:
                file.write(str(ij))
            file.close()
#####cif################
        if "cif" in output_files_type:
            if min_X_set == 1 or min_Y_set == 1:
                out_rcoa = [coa_conv2_rcoa(i, supercell_acod2, supercell_bcod2, [0, 0, slab_h + vacuum]) for i in
                            out_coa]
            # print("cif")
            ciflines = []
            #head lines
            ciflines.append("data_"+Filename)
            ciflines.append("\n\n#############################################")
            ciflines.append("\n# Generated by MOF-Membrane Constructor")
            ciflines.append("\n#############################################")

            ciflines.append("\n\n_cell_length_a            " + str(round(alen_input, 3)))
            ciflines.append("\n_cell_length_b            " + str(round(blen_input, 3)))
            ciflines.append("\n_cell_length_c            " + str(round(slab_h + vacuum, 3)))
            ciflines.append("\n_cell_angle_alpha		" + str(round(alpha, 2)))
            ciflines.append("\n_cell_angle_beta		" + str(round(beta, 2)))
            ciflines.append("\n_cell_angle_gamma		" + str(round(cang_input, 2)))
            ciflines.append("\n\n_symmetry_space_group_name_H-M		'P1'" )
            ciflines.append("\n_symmetry_Int_Tables_number		1")
            ciflines.append("\n_symmetry_cell_setting		Triclinic\n")
            ciflines.append("\nloop_\n_symmetry_equiv_pos_as_xyz\n'+x,+y,+z'\n")
            ciflines.append("\nloop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_charge")

            for i in range(0, len(out_rcoa)):
                if charge_r == 1:
                    incoord = '{0:<9}'.format(str(round(out_rcoa[i][0], 5))) + "   " + '{0:<9}'.format(str(round(out_rcoa[i][1], 5))) + "   " + '{0:<9}'.format(str(round(out_rcoa[i][2], 5))) + "      " + str(out_charge[i])
                elif charge_r == 0:
                    incoord = '{0:<9}'.format(str(round(out_rcoa[i][0], 5))) + "   " + '{0:<9}'.format(str(round(out_rcoa[i][1], 5))) + "   " + '{0:<9}'.format(str(round(out_rcoa[i][2], 5))) + "      " + "0.00 "
                ciflines.append("\n"+'{0:<5}'.format(out_elemt[i]) + '{0:<5}'.format(out_elemt[i]) + incoord)

            ciflines.append("\n")
            file = open(outpath + Filename + "_" + facet_index + ".cif", 'w')
            for ij in ciflines:
                file.write(str(ij))
            file.close()

#########main function start##############
    global charge_r
    charge_r = 1
    for dirpath, dirnames, filenames in os.walk(file_dir):

        print("dirpath,dirnames,filenames", dirpath, dirnames, filenames)
        total_len = len(filenames)
        p_i=0

        log_lines.append("########Job Parameters########\n\n")
        log_lines.append("Input Path: \n" + file_dir + "\n\n")
        log_lines.append("Membrane Thickness: " + str(slab_h) + "\n")
        log_lines.append("Vacuum Thickness: " + str(vacuum) + "\n")
        if min_X_set == 1:
            log_lines.append("Min X Length: " + str(min_X_length) + "\n")
        elif min_X_set == 0:
            log_lines.append("Min X Length: Not Set \n")
        if min_Y_set == 1:
            log_lines.append("Min Y Length: " + str(min_Y_length) + "\n")
        elif min_Y_set == 0:
            log_lines.append("Min Y Length: Not Set \n")

        if facet_mode == "SET":
            log_lines.append("Exposed Facet: Set; Facet index:" + facet_index_input + "\n")
        elif facet_mode == "AUTO":

            log_lines.append("Exposed Facet: Auto; Cell size:"+grid_size_input  + "\n")

        if protonation_yn == 0:
            log_lines.append("Surface Protonation: No\n")
        elif protonation_yn == 1:
            log_lines.append("Surface Protonation: Yes\n")
        if separate_layer == 0:
            log_lines.append("Barrier Layer: No\n")
        elif separate_layer == 1:
            log_lines.append("Barrier Layer: Yes\n")
        log_lines.append("\n##############################\n\n")


        for file_name in filenames:
            label_11_12.config(text=str(p_i) + "/" + str(total_len) + " (" + str(round(100 * p_i / total_len, 1)) + "%)")
            if running == False:
                label_11_22.config(text="Cancelled")
                break
            if running != False:
                try:
                    # xlsdata_each=[]
                    if file_name.find(".cif") > 0:
                        Filename = file_name[:-4]
                        label_11_22.config(text="Running")


                        log_lines.append("File Name:" + Filename+"\n")
                        dir_file_name = file_dir + file_name  # file_name = x.cif
                        print("path", dir_file_name)
                        # turn(dir_file_name)
                        with open(dir_file_name, encoding="utf-8") as cif_f:
                            cif_lines = cif_f.readlines()
                            list1 = []
                            structure_index = 28
                            for l_number in range(0, len(cif_lines)):
                                each_line_cif = cif_lines[l_number].split()

                                if len(each_line_cif) > 0:
                                    # print("%%%%%%%%%", re_zst(each_line_cif[0], "atom_site_charge"))
                                    if len(re.findall(r"_cell_length_a", each_line_cif[0])) > 0:
                                        pbc = cif_lines[l_number].split()  # a
                                        alen = float(pbc[1])

                                    if len(re.findall(r"_cell_length_b", each_line_cif[0])) > 0:
                                        pbc = cif_lines[l_number].split()  # b
                                        blen = float(pbc[1])

                                    if len(re.findall(r"_cell_length_c", each_line_cif[0])) > 0:
                                        pbc = cif_lines[l_number].split()  # c
                                        clen = float(pbc[1])

                                    if len(re.findall(r"_cell_angle_alpha", each_line_cif[0])) > 0:
                                        pbc = cif_lines[l_number].split()  # A
                                        aang = float(pbc[1])

                                    if len(re.findall(r"_cell_angle_beta", each_line_cif[0])) > 0:
                                        pbc = cif_lines[l_number].split()  # B
                                        bang = float(pbc[1])

                                    if len(re.findall(r"_cell_angle_gamma", each_line_cif[0])) > 0:
                                        pbc = cif_lines[l_number].split()  # C
                                        cang = float(pbc[1])

                                    # if len(re.findall(r"symmetry", each_line_cif[0])) > 0:
                                    #     print("##############", l_number)
                                    if len(re.findall(r"atom_site_fract_z", each_line_cif[0])) > 0:
                                        # print("enter...enter....")
                                        # print(cif_lines[l_number+1].split())
                                        if l_number + 1 < len(cif_lines) and len(re.findall(
                                                r'( [+-]?(?:[0-9]\d*)(?:\.\d*)? |[+-]?(?:[0-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+))',
                                                cif_lines[l_number + 1].split()[0])) > 0:
                                            structure_index = l_number + 1
                                            print("structure_index,l_number", l_number)
                                        elif l_number + 1 < len(cif_lines) and len(
                                                re.findall(r'charge', cif_lines[l_number + 1].split()[0])) > 0:
                                            structure_index = l_number + 2
                                            print("charge detected, structure_index,l_number+2", l_number + 2)
                                            print("atom_site_charges detected. They will be reserved in the output")
                                            charge_r = 1  ###raspa format
                                        elif l_number + 1 < len(cif_lines) and len(
                                                re.findall(r'U_iso_or_equiv', cif_lines[l_number + 1].split()[0])) > 0:
                                            structure_index = l_number + 4
                                            print("iso detected, structure_index,l_number+2", l_number + 4)
                                            ###MS format

                            #
                            aa = 3.14159265 * aang / 180
                            ba = 3.14159265 * bang / 180
                            # print(pbc[1])
                            ca = 3.14159265 * cang / 180
                            A = math.cos(aa)
                            B = math.cos(ba)
                            C = math.cos(ca)
                            h0 = clen * math.sqrt(1 - A * A - B * B - C * C + 2 * A * B * C) / math.sin(ca)

                            acod = [alen, 0, 0]
                            bcod = [blen * C, blen * math.sin(ca), 0]
                            ccodx = clen * B
                            ccodz = h0
                            ccody = (A * blen * clen - bcod[0] * ccodx - bcod[2] * ccodz) / bcod[1]
                            ccod = [ccodx, ccody, ccodz]
                            # print("ccod=", ccod)

                            structure = cif_lines[structure_index:]  #
                            N_of_end = structure.count("_end\n")
                            for jN in range(0, N_of_end):
                                structure.remove("_end\n")

                            N_of_n = structure.count("\n")
                            for iN in range(0, N_of_n):
                                structure.remove("\n")

                            elemt_all = []  # get the elements，[H,C,O,H,C,N,N]
                            coa_all = []  # coordinate of atom
                            charge_all = []


                            for line in structure:
                                # print("line:",line)
                                line_0 = re.findall(r'([A-Z][a-z]|[A-Z])', line)
                                # print("line_0",line_0)
                                if "loop_" in line_0:
                                    break
                                elif len(line_0)>1:
                                    elemt_single=line_0[0]
                                    # print("elemt_single",elemt_single)

                                line_t = re.sub("\t", " ", line)
                                line_t = re.sub(chr(9), "  ", line_t)
                                line_t = line_t.strip('\n') + "  "  # for the re.findall perform
                                # print("car lines:",line_t)
                                coa_charge = re.findall(
                                    r'( [+-]?(?:[0-9]\d*)(?:\.\d*)? |[+-]?(?:[0-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+))', line_t)
                                if len(coa_charge)>0:
                                    coordinate = coa_charge[0:3]

                                    coa_single = [0, 0, 0]
                                    # print("coordinate",coordinate)

                                    cood = [0, 0, 0]
                                    cood[0] = float(coordinate[0])
                                    cood[1] = float(coordinate[1])
                                    cood[2] = float(coordinate[2])

                                    # print(cood)
                                    # print(alen,blen,clen,aang,bang,cang)
                                    coa_single = rcoa_conv2_coa(cood, acod, bcod, ccod)

                                    elemt_all.append(elemt_single)  # elements
                                    coa_all.append(coa_single)
                                    if charge_r == 1:
                                        charge_each = coa_charge[3]
                                        charge_all.append(float(charge_each))
                                if len(coa_charge) == 0:
                                    break

                            print("len(coa_all), len(elemt_all)", len(coa_all), len(elemt_all))
                            #############determine the facet##################
                            if job_type==1:
                                if facet_mode == "SET":
                                    nonmetal_elemt_u = []  # elements without metal
                                    metal_elemt_u = []
                                    nonmetal_coa_u = []  # coordinate of atom, no metals
                                    metal_coa_u = []
                                    nonmetal_charge_u = []
                                    metal_charge_u = []
                                    #print("start metal recognition:")
                                    # xlsdata_each.append("null")
                                    for i in range(0, len(elemt_all)):
                                        if elemt_all[i] not in nonmetal_list:
                                            metal_elemt_u.append(elemt_all[i])
                                            metal_coa_u.append(coa_all[i])
                                            if charge_r == 1:
                                                metal_charge_u.append(charge_all[i])
                                        else:
                                            nonmetal_elemt_u.append(elemt_all[i])
                                            nonmetal_coa_u.append(coa_all[i])
                                            if charge_r == 1:
                                                nonmetal_charge_u.append(charge_all[i])
                                    print("atom number of nonmetal:", len(nonmetal_elemt_u))
                                    print("atom number of metal:", len(metal_elemt_u))

                                    # it is expected to get a value like this:
                                    # f_con=[(0, 4), (1, 6), (2, 3), (3, 2), (3, 12), (4, 0), (4, 10)...]
                                    f_con = []
                                    for i in range(0, len(nonmetal_elemt_u)):
                                        kk = []
                                        for j in range(i, len(nonmetal_elemt_u)):
                                            if j != i:
                                                d = p_distance_comput(nonmetal_coa_u[i], nonmetal_coa_u[j], acod, bcod, ccod, 1, 1,
                                                                      1)
                                                bond = float(atom_r[nonmetal_elemt_u[i]]) + float(atom_r[nonmetal_elemt_u[j]])
                                                if d <= bond:
                                                    kk.append(j)

                                        if len(kk) > 0:
                                            for ii in kk:
                                                con_single = (i, ii)
                                                f_con.append(con_single)
                                        elif len(kk) == 0:
                                            pass

                                    f_con0 = [i for i in f_con]
                                    L_index_u = con_find(f_con)
                                    metal_connect_structure_u = []
                                    L_index_list_u = []
                                    for i in range(0, len(metal_elemt_u)):
                                        metal_connect_structure_atoms = []
                                        for linker_i in range(0, len(
                                                L_index_u)):
                                            for linker_i_atomindex in L_index_u[linker_i]:
                                                L_index_list_u.append(linker_i_atomindex)
                                                d = p_distance_comput(metal_coa_u[i], nonmetal_coa_u[linker_i_atomindex], acod,
                                                                      bcod, ccod, 1, 1, 1)
                                                m_bond = float(atom_r[metal_elemt_u[i]]) + float(
                                                    atom_r[nonmetal_elemt_u[linker_i_atomindex]])
                                                if d <= m_bond:
                                                    metal_connect_structure_atoms.append(linker_i_atomindex)
                                                else:
                                                    pass
                                        metal_connect_structure_u.append(metal_connect_structure_atoms)

                                    # facet_index_input = "2 1 1"
                                    facet_str = re.findall(r'\d+', facet_index_input)
                                    facet_123 = [int(i) for i in facet_str]
                                    print("facet_123", facet_123)
                                    N = 1
                                    N2 = 1

                                    for i in facet_123:
                                        if i != 0:
                                            N = N * i

                                    if facet_123[0] >= 1:
                                        n0 = int(N / facet_123[0])
                                        coa_all_a = supercell_coa(coa_all, n0, acod)
                                        nonmetal_coa_a = supercell_coa(nonmetal_coa_u, n0,
                                                                       acod)  # coordinate of atom, no metals
                                        metal_coa_a = supercell_coa(metal_coa_u, n0, acod)
                                        f0 = "1"
                                        acod_N = [acod[0] * n0, acod[1] * n0, acod[2] * n0]
                                        N2 = n0 * N2
                                    elif facet_123[0] == 0:
                                        coa_all_a = coa_all
                                        nonmetal_coa_a = nonmetal_coa_u  # coordinate of atom, no metals
                                        metal_coa_a = metal_coa_u
                                        f0 = "0"
                                        acod_N = acod
                                    else:
                                        print("facet_input_error")

                                    if facet_123[1] >= 1:
                                        n1 = int(N / facet_123[1])
                                        coa_all_b = supercell_coa(coa_all_a, n1, bcod)
                                        nonmetal_coa_b = supercell_coa(nonmetal_coa_a, n1, bcod)
                                        metal_coa_b = supercell_coa(metal_coa_a, n1, bcod)
                                        f1 = "1"
                                        bcod_N = [bcod[0] * n1, bcod[1] * n1, bcod[2] * n1]
                                        N2 = n1 * N2

                                    elif facet_123[1] == 0:
                                        coa_all_b = coa_all_a
                                        nonmetal_coa_b = nonmetal_coa_a
                                        metal_coa_b = metal_coa_a
                                        f1 = "0"
                                        bcod_N = bcod
                                    else:
                                        print("facet_input_error")

                                    if facet_123[2] >= 1:
                                        n2 = int(N / facet_123[2])
                                        print("n2", n2)
                                        coa_all_N = supercell_coa(coa_all_b, n2, ccod)
                                        #nonmetal_coa_N = supercell_coa(nonmetal_coa_b, n2, ccod)
                                        #metal_coa_N = supercell_coa(metal_coa_b, n2, ccod)
                                        f2 = "1"
                                        ccod_N = [ccod[0] * n2, ccod[1] * n2, ccod[2] * n2]
                                        N2 = n2 * N2

                                    elif facet_123[2] == 0:
                                        coa_all_N = coa_all_b
                                        #nonmetal_coa_N = nonmetal_coa_b
                                        #metal_coa_N = metal_coa_b
                                        f2 = "0"
                                        ccod_N = ccod
                                    else:
                                        print("facet_input_error")

                                    facet_1 = f0 + " " + f1 + " " + f2
                                    print("facet_1", facet_1)

                                    elemt_all_N = N2 * elemt_all
                                    charge_all_N = N2 * charge_all
                                    print("start roration")
                                    coa_all_2, acod_2, bcod_2, ccod_2, lattice_2 = lattice_rotation(facet_1, coa_all_N, acod_N,
                                                                                                    bcod_N, ccod_N)
                                    print("rotation done")

                                    print(len(coa_all_2), len(elemt_all_N))

                                    print("start construction")
                                    facet_index = facet_str[0] + facet_str[1] + facet_str[2]

                                    ####test point2####


                                    slab_construction(elemt_all_N, coa_all_2, acod_2, bcod_2, ccod_2, lattice_2, 1)
                                    print("construction done")
                                elif facet_mode == "AUTO":
                                    facet = "0 0 1"  # 0
                                    coa_all_R001, acod_R001, bcod_R001, ccod_R001, lattice_R001 = lattice_rotation(facet, coa_all,
                                                                                                                   acod,
                                                                                                                   bcod, ccod)
                                    pld001_1 = round(pld_comput(coa_all_R001, acod_R001, bcod_R001, ccod_R001,1), 3)
                                    pld100_2 = round(pld_comput(coa_all_R001, acod_R001, bcod_R001, ccod_R001,2), 3)


                                    facet = "0 1 0"  # 5
                                    coa_all_R010, acod_R010, bcod_R010, ccod_R010, lattice_R010 = lattice_rotation(
                                        facet, coa_all,
                                        acod,
                                        bcod, ccod)
                                    pld010_1 = round(pld_comput(coa_all_R010, acod_R010, bcod_R010, ccod_R010,1), 3)

                                    facet = "1 0 0"  # 6
                                    coa_all_R100, acod_R100, bcod_R100, ccod_R100, lattice_R100 = lattice_rotation(
                                        facet, coa_all,
                                        acod,
                                        bcod, ccod)

                                    pld100_1 = round(pld_comput(coa_all_R100, acod_R100, bcod_R100, ccod_R100,1), 3)
                                    pld001_2 = round(pld_comput(coa_all_R100, acod_R100, bcod_R100, ccod_R100,2), 3)

                                    facet = "0 1 1"  # 3
                                    coa_all_R011, acod_R011, bcod_R011, ccod_R011, lattice_R011 = lattice_rotation(facet, coa_all,
                                                                                                                   acod,
                                                                                                                   bcod, ccod)
                                    pld011 = round(pld_comput(coa_all_R011, acod_R011, bcod_R011, ccod_R011,1), 3)

                                    facet = "1 0 1"  # 2
                                    coa_all_R101, acod_R101, bcod_R101, ccod_R101, lattice_R101 = lattice_rotation(
                                        facet, coa_all,
                                        acod,
                                        bcod, ccod)
                                    pld101 = round(pld_comput(coa_all_R101, acod_R101, bcod_R101, ccod_R101,1), 3)
                                    pld010_2 = round(pld_comput(coa_all_R101, acod_R101, bcod_R101, ccod_R101, 2), 3)


                                    facet = "1 1 0"  # 1
                                    coa_all_R110, acod_R110, bcod_R110, ccod_R110, lattice_R110 = lattice_rotation(
                                        facet, coa_all,
                                        acod,
                                        bcod, ccod)
                                    pld110 = round(pld_comput(coa_all_R110, acod_R110, bcod_R110, ccod_R110,1), 3)

                                    facet = "1 1 1"  # 0
                                    coa_all_R111, acod_R111, bcod_R111, ccod_R111, lattice_R111 = lattice_rotation(facet, coa_all,
                                                                                                                   acod,
                                                                                                                   bcod, ccod)
                                    pld111 = round(pld_comput(coa_all_R111, acod_R111, bcod_R111, ccod_R111,1), 3)


                                    pld001 = max(pld001_1,pld001_2)
                                    pld010 = max(pld010_1, pld010_2)
                                    pld100 = max(pld100_1, pld100_2)
                                    # print("pld001,pld010,pld100,pld011,pld101,pld110,pld111:", pld001, pld010, pld100, pld011,
                                    #       pld101,
                                    #       pld110, pld111)

                                    facet_list = ["001", "010", "100", "011", "101", "110", "111"]


                                    pld_list = [pld001, pld010, pld100, pld011, pld101, pld110, pld111]
                                    pld_max = max(pld_list)
                                    # log_lines.append("Estimated pore_001: " + str(pld001) + "\n")
                                    # log_lines.append("Estimated pore_010: " + str(pld010) + "\n")
                                    # log_lines.append("Estimated pore_100: " + str(pld100) + "\n")
                                    # log_lines.append("Estimated pore_110: " + str(pld110) + "\n")
                                    # log_lines.append("Estimated pore_101: " + str(pld101) + "\n")
                                    # log_lines.append("Estimated pore_011: " + str(pld011) + "\n")
                                    # log_lines.append("Estimated pore_111: " + str(pld111) + "\n")

                                    facet_list_index = pld_list.index(pld_max)
                                    facet_index = facet_list[facet_list_index]
                                    print("facet_index = ", facet_index)
                                    charge_all_N = [i for i in charge_all]

                                    ################重新定义晶体参数##################
                                    lattice_2 = []
                                    if facet_list_index==0:
                                        log_lines.append("Selected facet: 001 \n" )
                                        coa_all_2 = coa_all_R001
                                        acod_2 = acod_R001
                                        bcod_2 = bcod_R001
                                        ccod_2 = ccod_R001
                                        lattice_2 = lattice_R001
                                    elif facet_list_index==1:
                                        log_lines.append("Selected facet: 010 \n")
                                        coa_all_2 = coa_all_R010
                                        acod_2 = acod_R010
                                        bcod_2 = bcod_R010
                                        ccod_2 = ccod_R010
                                        lattice_2 = lattice_R010
                                    elif facet_list_index == 2:
                                        log_lines.append("Selected facet: 100 \n")
                                        coa_all_2 = coa_all_R100
                                        acod_2 = acod_R100
                                        bcod_2 = bcod_R100
                                        ccod_2 = ccod_R100
                                        lattice_2 = lattice_R100
                                    elif facet_list_index == 3:
                                        log_lines.append("Selected facet: 011 \n")
                                        coa_all_2 = coa_all_R011
                                        acod_2 = acod_R011
                                        bcod_2 = bcod_R011
                                        ccod_2 = ccod_R011
                                        lattice_2 = lattice_R011
                                    elif facet_list_index == 4:
                                        log_lines.append("Selected facet: 101 \n")
                                        coa_all_2 = coa_all_R101
                                        acod_2 = acod_R101
                                        bcod_2 = bcod_R101
                                        ccod_2 = ccod_R101
                                        lattice_2 = lattice_R101
                                    elif facet_list_index == 5:
                                        log_lines.append("Selected facet: 110 \n")
                                        coa_all_2 = coa_all_R110
                                        acod_2 = acod_R110
                                        bcod_2 = bcod_R110
                                        ccod_2 = ccod_R110
                                        lattice_2 = lattice_R110
                                    elif facet_list_index == 6:
                                        log_lines.append("Selected facet: 111 \n")
                                        coa_all_2 = coa_all_R111
                                        acod_2 = acod_R111
                                        bcod_2 = bcod_R111
                                        ccod_2 = ccod_R111
                                        lattice_2 = lattice_R111
                                    slab_construction(elemt_all, coa_all_2, acod_2, bcod_2, ccod_2, lattice_2, 0)

                    if file_name.find(".car") > 0:

                        Filename = file_name[:-4]
                        label_11_22.config(text="Running")
                        log_lines.append("File Name:" + Filename + "\n")

                        dir_file_name = file_dir + file_name  # file_name = dasda.cif
                        print("path", dir_file_name)
                        charge_r = 1
                        acod = [0, 0, 0]
                        bcod = [0, 0, 0]
                        ccod = [0, 0, 0]

                        with open(dir_file_name, encoding="utf-8") as car_f:
                            car_lines = car_f.readlines()
                            #              ####
                            for l_number in range(0, len(car_lines)):
                                # print(cif_lines[l_number])
                                if len(re.findall(r"PBC", car_lines[l_number])) > 0 and len(
                                        re.findall(r"[0-9]", car_lines[l_number])) > 0:
                                    # print(car_lines[l_number])
                                    car_lines_line_t = re.sub("\t", " ", car_lines[l_number])
                                    car_lines_line_t = car_lines_line_t.strip('\n') + "  "
                                    lattice = re.findall(
                                        r'([+-]?(?:[0-9]\d*)(?:\.\d*)? |[+-]?(?:[0-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+))',
                                        car_lines_line_t)
                                    # print(lattice)
                                    alen = float(lattice[0])
                                    blen = float(lattice[1])
                                    clen = float(lattice[2])
                                    aang = float(lattice[3])
                                    bang = float(lattice[4])
                                    cang = float(lattice[5])
                                    print("lattice:", alen, blen, clen, aang, bang, cang)

                            aa = 3.14159265 * aang / 180
                            ba = 3.14159265 * bang / 180
                            # print(pbc[1])
                            ca = 3.14159265 * cang / 180
                            A = math.cos(aa)
                            B = math.cos(ba)
                            C = math.cos(ca)

                            h0 = clen * math.sqrt(1 - A * A - B * B - C * C + 2 * A * B * C) / math.sin(ca)

                            acod = [alen, 0, 0]
                            bcod = [blen * C, blen * math.sin(ca), 0]
                            ccodx = clen * B
                            ccodz = h0
                            ccody = (A * blen * clen - bcod[0] * ccodx - bcod[2] * ccodz) / bcod[1]
                            ccod = [ccodx, ccody, ccodz]
                            #print("ccod=", ccod)

                            structure = car_lines[5:]  #
                            N_of_end = structure.count("end\n")
                            for jN in range(0, N_of_end):
                                structure.remove("end\n")
                            # print(N_of_end)

                            N_of_n = structure.count("\n")
                            for iN in range(0, N_of_n):
                                structure.remove("\n")

                            # print("*****")
                            # print(structure)
                            elemt_all = []  # get the elements
                            coa_all = []  # coordinate of atom
                            charge_all = []

                            for line in structure:

                                # print(L[0:9])
                                # print(elemt_single)

                                line_t = re.sub("\t", " ", line)
                                line_t = line_t.strip('\n') + "  "  #modify for the re.findall
                                # print("car lines:",line_t)
                                L = line_t.split()
                                # print("L[0:9]",L[0:9])
                                elemt_single = L[0:9][7]
                                coa_charge = re.findall(
                                    r'( [+-]?(?:[0-9]\d*)(?:\.\d*)? |[+-]?(?:[0-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+))', line_t)
                                coordinate = coa_charge[0:3]
                                # print("coa_charge_car:",coa_charge)
                                if len(coa_charge) >= 5:
                                    charge_r = 1
                                else: charge_r = 0

                                coa_single = [0, 0, 0]
                                coa_single[0] = float(coordinate[0])
                                coa_single[1] = float(coordinate[1])
                                coa_single[2] = float(coordinate[2])
                                # print("coa_single:",coa_single)

                                elemt_all.append(elemt_single)  # elements
                                coa_all.append(coa_single)  # coordinations
                                if charge_r == 1:
                                    charge_each = coa_charge[4]  # car 4
                                    charge_all.append(float(charge_each))

                            #############判断晶面##################
                            if job_type==1:
                                if facet_mode == "SET":

                                    nonmetal_elemt_u = []  # elements without metal
                                    metal_elemt_u = []
                                    nonmetal_coa_u = []  # coordinate of atom, no metals
                                    metal_coa_u = []
                                    nonmetal_charge_u = []
                                    metal_charge_u = []
                                    #print("start metal recognition:")
                                    for i in range(0, len(elemt_all)):
                                        if elemt_all[i] not in nonmetal_list:
                                            metal_elemt_u.append(elemt_all[i])
                                            metal_coa_u.append(coa_all[i])
                                            if charge_r == 1:
                                                metal_charge_u.append(charge_all[i])
                                        else:
                                            nonmetal_elemt_u.append(elemt_all[i])
                                            nonmetal_coa_u.append(coa_all[i])
                                            if charge_r == 1:
                                                nonmetal_charge_u.append(charge_all[i])
                                    print("atom number of nonmetal:", len(nonmetal_elemt_u))
                                    print("atom number of metal:", len(metal_elemt_u))

                                    # it is expected to get a value like this:
                                    # f_con=[(0, 4), (1, 6), (2, 3), (3, 2), (3, 12), (4, 0), (4, 10)...]
                                    f_con = []
                                    for i in range(0, len(nonmetal_elemt_u)):
                                        kk = []
                                        for j in range(i, len(nonmetal_elemt_u)):
                                            if j != i:
                                                d = p_distance_comput(nonmetal_coa_u[i], nonmetal_coa_u[j], acod, bcod, ccod, 1, 1,
                                                                      1)
                                                bond = float(atom_r[nonmetal_elemt_u[i]]) + float(atom_r[nonmetal_elemt_u[j]])
                                                if d <= bond:
                                                    kk.append(j)

                                        if len(kk) > 0:
                                            for ii in kk:
                                                con_single = (i, ii)
                                                f_con.append(con_single)
                                        elif len(kk) == 0:
                                            pass

                                    f_con0 = [i for i in f_con]
                                    L_index_u = con_find(f_con)
                                    metal_connect_structure_u = []
                                    L_index_list_u = []
                                    for i in range(0, len(metal_elemt_u)):
                                        metal_connect_structure_atoms = []
                                        for linker_i in range(0, len(
                                                L_index_u)):
                                            for linker_i_atomindex in L_index_u[linker_i]:
                                                L_index_list_u.append(linker_i_atomindex)
                                                d = p_distance_comput(metal_coa_u[i], nonmetal_coa_u[linker_i_atomindex], acod,
                                                                      bcod, ccod, 1, 1, 1)
                                                m_bond = float(atom_r[metal_elemt_u[i]]) + float(
                                                    atom_r[nonmetal_elemt_u[linker_i_atomindex]])
                                                if d <= m_bond:
                                                    metal_connect_structure_atoms.append(linker_i_atomindex)
                                                else:
                                                    pass
                                        metal_connect_structure_u.append(metal_connect_structure_atoms)

                                    # facet_index_input = "2 1 1"
                                    facet_str = re.findall(r'\d+', facet_index_input)
                                    facet_123 = [int(i) for i in facet_str]
                                    print("facet_123", facet_123)
                                    N = 1
                                    N2 = 1

                                    for i in facet_123:
                                        if i != 0:
                                            N = N * i

                                    if facet_123[0] >= 1:
                                        n0 = int(N / facet_123[0])
                                        coa_all_a = supercell_coa(coa_all, n0, acod)
                                        nonmetal_coa_a = supercell_coa(nonmetal_coa_u, n0,
                                                                       acod)  # coordinate of atom, no metals
                                        metal_coa_a = supercell_coa(metal_coa_u, n0, acod)
                                        f0 = "1"
                                        acod_N = [acod[0] * n0, acod[1] * n0, acod[2] * n0]
                                        N2 = n0 * N2
                                    elif facet_123[0] == 0:
                                        coa_all_a = coa_all
                                        nonmetal_coa_a = nonmetal_coa_u  # coordinate of atom, no metals
                                        metal_coa_a = metal_coa_u
                                        f0 = "0"
                                        acod_N = acod
                                    else:
                                        print("facet_input_error")

                                    if facet_123[1] >= 1:
                                        n1 = int(N / facet_123[1])
                                        coa_all_b = supercell_coa(coa_all_a, n1, bcod)
                                        nonmetal_coa_b = supercell_coa(nonmetal_coa_a, n1, bcod)
                                        metal_coa_b = supercell_coa(metal_coa_a, n1, bcod)
                                        f1 = "1"
                                        bcod_N = [bcod[0] * n1, bcod[1] * n1, bcod[2] * n1]
                                        N2 = n1 * N2

                                    elif facet_123[1] == 0:
                                        coa_all_b = coa_all_a
                                        nonmetal_coa_b = nonmetal_coa_a
                                        metal_coa_b = metal_coa_a
                                        f1 = "0"
                                        bcod_N = bcod
                                    else:
                                        print("facet_input_error")

                                    if facet_123[2] >= 1:
                                        n2 = int(N / facet_123[2])
                                        print("n2", n2)
                                        coa_all_N = supercell_coa(coa_all_b, n2, ccod)
                                        #nonmetal_coa_N = supercell_coa(nonmetal_coa_b, n2, ccod)
                                        #metal_coa_N = supercell_coa(metal_coa_b, n2, ccod)
                                        f2 = "1"
                                        ccod_N = [ccod[0] * n2, ccod[1] * n2, ccod[2] * n2]
                                        N2 = n2 * N2

                                    elif facet_123[2] == 0:
                                        coa_all_N = coa_all_b
                                        #nonmetal_coa_N = nonmetal_coa_b
                                        #metal_coa_N = metal_coa_b
                                        f2 = "0"
                                        ccod_N = ccod
                                    else:
                                        print("facet_input_error")

                                    facet_1 = f0 + " " + f1 + " " + f2
                                    print("facet_1", facet_1)

                                    elemt_all_N = N2 * elemt_all
                                    charge_all_N = N2 * charge_all
                                    print("start roration")
                                    coa_all_2, acod_2, bcod_2, ccod_2, lattice_2 = lattice_rotation(facet_1, coa_all_N, acod_N,
                                                                                                    bcod_N, ccod_N)
                                    print("rotation done")

                                    print(len(coa_all_2), len(elemt_all_N))

                                    print("start construction")
                                    facet_index = facet_str[0] + facet_str[1] + facet_str[2]

                                    print("facet_index ",facet_index )

                                    slab_construction(elemt_all_N, coa_all_2, acod_2, bcod_2, ccod_2, lattice_2, 1)
                                    print("construction done")
                                elif facet_mode == "AUTO":
                                    facet = "0 0 1"  # 0
                                    coa_all_R001, acod_R001, bcod_R001, ccod_R001, lattice_R001 = lattice_rotation(
                                        facet, coa_all,
                                        acod,
                                        bcod, ccod)
                                    pld001_1 = round(pld_comput(coa_all_R001, acod_R001, bcod_R001, ccod_R001, 1), 3)
                                    pld100_2 = round(pld_comput(coa_all_R001, acod_R001, bcod_R001, ccod_R001, 2), 3)

                                    facet = "0 1 0"  # 5
                                    coa_all_R010, acod_R010, bcod_R010, ccod_R010, lattice_R010 = lattice_rotation(
                                        facet, coa_all,
                                        acod,
                                        bcod, ccod)
                                    pld010_1 = round(pld_comput(coa_all_R010, acod_R010, bcod_R010, ccod_R010, 1), 3)

                                    facet = "1 0 0"  # 6
                                    coa_all_R100, acod_R100, bcod_R100, ccod_R100, lattice_R100 = lattice_rotation(
                                        facet, coa_all,
                                        acod,
                                        bcod, ccod)

                                    pld100_1 = round(pld_comput(coa_all_R100, acod_R100, bcod_R100, ccod_R100, 1), 3)
                                    pld001_2 = round(pld_comput(coa_all_R100, acod_R100, bcod_R100, ccod_R100, 2), 3)

                                    facet = "0 1 1"  # 3
                                    coa_all_R011, acod_R011, bcod_R011, ccod_R011, lattice_R011 = lattice_rotation(
                                        facet, coa_all,
                                        acod,
                                        bcod, ccod)
                                    pld011 = round(pld_comput(coa_all_R011, acod_R011, bcod_R011, ccod_R011, 1), 3)

                                    facet = "1 0 1"  # 2
                                    coa_all_R101, acod_R101, bcod_R101, ccod_R101, lattice_R101 = lattice_rotation(
                                        facet, coa_all,
                                        acod,
                                        bcod, ccod)
                                    pld101 = round(pld_comput(coa_all_R101, acod_R101, bcod_R101, ccod_R101, 1), 3)
                                    pld010_2 = round(pld_comput(coa_all_R101, acod_R101, bcod_R101, ccod_R101, 2), 3)

                                    facet = "1 1 0"  # 1
                                    coa_all_R110, acod_R110, bcod_R110, ccod_R110, lattice_R110 = lattice_rotation(
                                        facet, coa_all,
                                        acod,
                                        bcod, ccod)
                                    pld110 = round(pld_comput(coa_all_R110, acod_R110, bcod_R110, ccod_R110, 1), 3)

                                    facet = "1 1 1"  # 0
                                    coa_all_R111, acod_R111, bcod_R111, ccod_R111, lattice_R111 = lattice_rotation(
                                        facet, coa_all,
                                        acod,
                                        bcod, ccod)
                                    pld111 = round(pld_comput(coa_all_R111, acod_R111, bcod_R111, ccod_R111, 1), 3)

                                    pld001 = max(pld001_1, pld001_2)
                                    pld010 = max(pld010_1, pld010_2)
                                    pld100 = max(pld100_1, pld100_2)

                                    facet_list = ["001", "010", "100", "011", "101", "110", "111"]
                                    pld_list = [pld001, pld010, pld100, pld011, pld101, pld110, pld111]
                                    pld_max = max(pld_list)
                                    # log_lines.append("Estimated pore_001: "+str(pld001)+"\n")
                                    # log_lines.append("Estimated pore_010: " + str(pld010) + "\n")
                                    # log_lines.append("Estimated pore_100: " + str(pld100) + "\n")
                                    # log_lines.append("Estimated pore_110: " + str(pld110) + "\n")
                                    # log_lines.append("Estimated pore_101: " + str(pld101) + "\n")
                                    # log_lines.append("Estimated pore_011: " + str(pld011) + "\n")
                                    # log_lines.append("Estimated pore_111: " + str(pld111) + "\n")
                                    facet_list_index = pld_list.index(pld_max)
                                    facet_index = facet_list[facet_list_index]
                                    print("facet_index = ", facet_index)
                                    charge_all_N = [i for i in charge_all]

                                    ################redefine the lattice##################
                                    lattice_2 = []
                                    if facet_list_index==0:
                                        log_lines.append("Selected facet: 001\n")
                                        coa_all_2 = coa_all_R001
                                        acod_2 = acod_R001
                                        bcod_2 = bcod_R001
                                        ccod_2 = ccod_R001
                                        lattice_2 = lattice_R001
                                    elif facet_list_index==1:
                                        log_lines.append("Selected facet: 010\n")
                                        coa_all_2 = coa_all_R010
                                        acod_2 = acod_R010
                                        bcod_2 = bcod_R010
                                        ccod_2 = ccod_R010
                                        lattice_2 = lattice_R010
                                    elif facet_list_index == 2:
                                        log_lines.append("Selected facet: 100\n")
                                        coa_all_2 = coa_all_R100
                                        acod_2 = acod_R100
                                        bcod_2 = bcod_R100
                                        ccod_2 = ccod_R100
                                        lattice_2 = lattice_R100
                                    elif facet_list_index == 3:
                                        log_lines.append("Selected facet: 011\n")
                                        coa_all_2 = coa_all_R011
                                        acod_2 = acod_R011
                                        bcod_2 = bcod_R011
                                        ccod_2 = ccod_R011
                                        lattice_2 = lattice_R011
                                    elif facet_list_index == 4:
                                        log_lines.append("Selected facet: 101\n")
                                        coa_all_2 = coa_all_R101
                                        acod_2 = acod_R101
                                        bcod_2 = bcod_R101
                                        ccod_2 = ccod_R101
                                        lattice_2 = lattice_R101
                                    elif facet_list_index == 5:
                                        log_lines.append("Selected facet: 110\n")
                                        coa_all_2 = coa_all_R110
                                        acod_2 = acod_R110
                                        bcod_2 = bcod_R110
                                        ccod_2 = ccod_R110
                                        lattice_2 = lattice_R110
                                    elif facet_list_index == 6:
                                        log_lines.append("Selected facet: 111\n")
                                        coa_all_2 = coa_all_R111
                                        acod_2 = acod_R111
                                        bcod_2 = bcod_R111
                                        ccod_2 = ccod_R111
                                        lattice_2 = lattice_R111

                                    slab_construction(elemt_all, coa_all_2, acod_2, bcod_2, ccod_2, lattice_2, 0)


                except Exception as e:
                    error_lines.append(Filename+'  '+str(e)+"\n")
            p_i = p_i + 1
            label_11_12.config(text=str(p_i) + "/" + str(total_len)+" (" + str(round(100 * p_i / total_len, 1)) + "%)")

        e_file = open(outpath + "error.log", 'w')
        for ij in error_lines:
            e_file.write(str(ij))
        e_file.close()



        log_file = open(outpath + "job.log", 'w')
        for ij in log_lines:
            log_file.write(str(ij))
        log_file.close()
        end = time.time()

        print('Running time: %s Seconds' % (end - start_t))
        label_11_22.config(text='Done ('+ str(round(end - start_t,1))+" s)")

def stop():
    global running  # create global
    running = False
    label_11_22.config(text="Cancelling")

def selectPath():
    global file_dir
    path_ = askdirectory(initialdir=file_dir)
    path.set(path_)
    file_dir = path_ +"/"
    file_dir2 = path_ + "/"
    mark_i=0
    global out_ini
    if file_dir2.count("/")>=2:
        for i in range(0,len(file_dir2)):
            if file_dir2[-1]=="/" and mark_i==0:
                file_dir2=file_dir2.strip("/")
                mark_i=1
            elif file_dir2[-1]!="/":
                file_dir2=file_dir2.strip(file_dir2[-1])
            elif file_dir2[-1]=="/" and mark_i==1:
                break

        out_ini=file_dir2
    else:
        out_ini=file_dir2


def outPath():
    global outpath
    path_2 = askdirectory(initialdir=out_ini)
    path2.set(path_2)
    outpath=path_2+"/"
    print("output_path:", outpath)
    global out_ini2


def button8_select():
    global output_files_type
    if var8_1.get() == 1:
        output_files_type.append("pdb")
        output_files_type = list(set(output_files_type))

    elif var8_1.get() == 0:
        N = output_files_type.count("pdb")
        for i in range(0,N):
            output_files_type.remove("pdb")
    if var8_2.get() == 1:
        output_files_type.append("cif")
        output_files_type = list(set(output_files_type))

    elif var8_2.get() == 0:
        N = output_files_type.count("cif")
        for i in range(0,N):
            output_files_type.remove("cif")

    if var8_3.get() == 1:
        output_files_type.append("car")
        output_files_type = list(set(output_files_type))

    elif var8_3.get() == 0:
        N = output_files_type.count("car")
        for i in range(0,N):
            output_files_type.remove("car")

win = tkinter.Tk()
win.title(string="MOF-Membrane Constructor")

frame1 = tkinter.Frame(win, relief=tkinter.RAISED,borderwidth=2,width=900,height=100)
frame1.pack(side=tkinter.TOP,fill=tkinter.BOTH,ipadx=13,ipady=13,expand=1)

frame2 = tkinter.Frame(win, relief=tkinter.RAISED, borderwidth=2,width=800,height=260)
frame2.pack(side=tkinter.TOP,fill=tkinter.BOTH,ipadx=13,ipady=13,expand=1)

frame3 = tkinter.Frame(win, relief=tkinter.RAISED, borderwidth=2,width=800,height=100)
frame3.pack(side=tkinter.BOTTOM, fill=tkinter.BOTH, ipadx=13, ipady=13, expand=1)

startButton = tkinter.Button(frame3, height=2, width=15, text="Start",command=button_click,font=('Arial', '12'))  # Change to call button_click instead start
stopButton = tkinter.Button(frame3, height=2, width=15, text="Stop", command=stop,font=('Arial', '12'))
startButton.place(x=300,y=85,anchor=tkinter.CENTER,height=40)
stopButton.place(x=600,y=85,anchor=tkinter.CENTER,height=40)
###################################
path = tkinter.StringVar()

label_1 = tkinter.Label(frame1,text = "Input Path",font=('Arial', '12'))
entry = tkinter.Entry(frame1, textvariable = path)
button1 = tkinter.Button(frame1,text="Select",command = selectPath,font=('Arial', '12'))

label_1.place(x=150,y=30,anchor=tkinter.CENTER,width=300,height=40)
entry.place(x=30,y=70,anchor=tkinter.W,width=150,height=40)
button1.place(x=185,y=70,anchor=tkinter.W,width=100,height=40)

############

path2 = tkinter.StringVar()

label_2 = tkinter.Label(frame1,text = "Output Path",font=('Arial', '12'))
entry_2 = tkinter.Entry(frame1, textvariable = path2)
button2 = tkinter.Button(frame1,text="Select",command = outPath,font=('Arial', '12'))

label_2.place(x=450,y=30,anchor=tkinter.CENTER,width=300,height=40)
entry_2.place(x=330,y=70,anchor=tkinter.W,width=150,height=40)
button2.place(x=485,y=70,anchor=tkinter.W,width=100,height=40)
###########

label_3 = tkinter.Label(frame1,text = "Output Type",font=('Arial', '12'))
label_3.place(x=770,y=30,anchor=tkinter.CENTER,width=300,height=40)

var8_1=tkinter.IntVar()
button8_1=tkinter.Checkbutton(frame1,text="*.pdb",variable=var8_1,onvalue=1, offvalue=0,command=button8_select,font=('Arial', '12'))
button8_1.place(x=650,y=70,anchor=tkinter.W,width=90,height=40)

var8_2=tkinter.IntVar(value=1)
button8_2=tkinter.Checkbutton(frame1,text="*.cif",variable=var8_2,onvalue=1, offvalue=0,command=button8_select,font=('Arial', '12'))
button8_2.place(x=730,y=70,anchor=tkinter.W,width=90,height=40)

var8_3=tkinter.IntVar()
button8_3=tkinter.Checkbutton(frame1,text="*.car",variable=var8_3,onvalue=1, offvalue=0,command=button8_select,font=('Arial', '12'))
button8_3.place(x=810,y=70,anchor=tkinter.W,width=90,height=40)
############################

#########
expression3 = tkinter.StringVar(value="45")
entry_3=tkinter.Entry(frame2,textvariable=expression3,font=('Arial', '12'))
label_3 = tkinter.Label(frame2,text="Membrane Thickness",font=('Arial', '12'))
label_3_2 = tkinter.Label(frame2,text="(Å)",font=('Arial', '12'))

label_3.place(x=40,y=40,anchor=tkinter.W,width=180,height=40)
entry_3.place(x=220,y=40,anchor=tkinter.W,width=90,height=40)
label_3_2.place(x=315,y=40,anchor=tkinter.W,width=50,height=40)
#########

expression4 = tkinter.StringVar(value="100")

entry_4 = tkinter.Entry(frame2,textvariable=expression4,font=('Arial', '12'))
label_4 = tkinter.Label(frame2,text="Vacuum Thickness",font=('Arial', '12'))
label_4_2 = tkinter.Label(frame2,text="(Å)",font=('Arial', '12'))

label_4.place(x=40,y=100,anchor=tkinter.W,width=180,height=40)
entry_4.place(x=220,y=100,anchor=tkinter.W,width=90,height=40)
label_4_2.place(x=315,y=100,anchor=tkinter.W,width=50,height=40)

########################

expression5 = tkinter.StringVar(value="25")

def button5_select():
    global min_X_set
    if var5.get()==1:
        print("set X-size")
        # label_5.config(foreground="#000000")
        entry_5.config(foreground="#000000")
        label_5_2.config(foreground="#000000")
        min_X_set=1

    elif var5.get()==0:
        print("not set")
        # label_5.config(foreground="#969696")
        entry_5.config(foreground="#969696")
        label_5_2.config(foreground="#969696")
        min_X_set=0

var5=tkinter.IntVar()
button5 = tkinter.Checkbutton(frame2,variable=var5,onvalue=1, offvalue=0,command=button5_select )
entry_5 = tkinter.Entry(frame2,textvariable=expression5,foreground="#969696",font=('Arial', '12'))
label_5 = tkinter.Label(frame2,text="Min X Length",font=('Arial', '12'))
label_5_2 = tkinter.Label(frame2,text="(Å)",foreground="#969696",font=('Arial', '12'))


button5.place(x=30,y=160,anchor=tkinter.W,width=20,height=40)

label_5.place(x=40,y=160,anchor=tkinter.W,width=180,height=40)
entry_5.place(x=220,y=160,anchor=tkinter.W,width=90,height=40)
label_5_2.place(x=315,y=160,anchor=tkinter.W,width=50,height=40)

###################

########################

expression6 = tkinter.StringVar(value="25")

def button6_select():
    global min_Y_set
    if var6.get()==1:
        print("set Y-size")
        # label_6.config(foreground="#000000")
        entry_6.config(foreground="#000000")
        label_6_2.config(foreground="#000000")
        min_Y_set=1
    elif var6.get()==0:
        print("not set")
        # label_6.config(foreground="#969696")
        entry_6.config(foreground="#969696")
        label_6_2.config(foreground="#969696")
        min_Y_set=0

var6=tkinter.IntVar(value=0)
button6 = tkinter.Checkbutton(frame2,variable=var6,onvalue=1, offvalue=0,command=button6_select )
entry_6 = tkinter.Entry(frame2,textvariable=expression6,foreground="#969696",font=('Arial', '12'))
label_6 = tkinter.Label(frame2,text="Min Y Length",font=('Arial', '12'))
label_6_2 = tkinter.Label(frame2,text="(Å)",foreground="#969696",font=('Arial', '12'))

button6.place(x=30,y=220,anchor=tkinter.W,width=20,height=40)
label_6.place(x=40,y=220,anchor=tkinter.W,width=180,height=40)
entry_6.place(x=220,y=220,anchor=tkinter.W,width=90,height=40)
label_6_2.place(x=315,y=220,anchor=tkinter.W,width=50,height=40)
#####################
expression7 = tkinter.StringVar(value="0 0 1")
expression7_1 = tkinter.StringVar(value="0.8")

def button7_select():
    global facet_mode
    if var7.get()==1:
        print("auto")
        # label_6.config(foreground="#000000")
        entry_7_2.config(foreground="#969696")
        label_7_2.config(foreground="#969696")
        label_7_2_1.config(foreground="#969696")
        button7_2.config(foreground="#969696")
        button7_1.config(foreground="#000000")
        entry_7_1.config(foreground="#000000")
        label_7_1.config(foreground="#000000")
        button7_1.config(foreground="#000000")
        label_7_1_1.config(foreground="#000000")

        facet_mode="AUTO"

    elif var7.get()==0:
        print("set")
        # label_6.config(foreground="#969696")
        entry_7_2.config(foreground="#000000")
        button7_1.config(foreground="#969696")
        button7_2.config(foreground="#000000")
        label_7_2.config(foreground="#000000")
        entry_7_1.config(foreground="#969696")
        label_7_1.config(foreground="#969696")
        button7_1.config(foreground="#969696")
        label_7_2_1.config(foreground="#000000")
        label_7_1_1.config(foreground="#969696")

        facet_mode="SET"



var7=tkinter.IntVar(value=1)
label_7 = tkinter.Label(frame2,text="Exposed Facet",font=('Arial', '12'))
button7_1 = tkinter.Checkbutton(frame2,text="Auto",variable=var7,onvalue=1, offvalue=0,command=button7_select,font=('Arial', '12') )
button7_2 = tkinter.Checkbutton(frame2,text="Set ",variable=var7,onvalue=0, offvalue=1,command=button7_select,foreground="#969696",font=('Arial', '12'))
entry_7_1 = tkinter.Entry(frame2,text="0.8",textvariable=expression7_1,foreground="#000000",font=('Arial', '12'))
label_7_1 = tkinter.Label(frame2,text="(Å)",foreground="#000000",font=('Arial', '12'))
label_7_1_1 = tkinter.Label(frame2,text="cell size",foreground="#000000",font=('Arial', '12'))

entry_7_2 = tkinter.Entry(frame2,text="0 0 1",textvariable=expression7,foreground="#969696",font=('Arial', '12'))
label_7_2 = tkinter.Label(frame2,text="(h k l)",foreground="#969696",font=('Arial', '12'))
label_7_2_1 = tkinter.Label(frame2,text="facet index",foreground="#969696",font=('Arial', '12'))

label_7.place(x=430,y=40,anchor=tkinter.W,width=160,height=40)
button7_1.place(x=470,y=80,anchor=tkinter.W,width=70,height=40)
entry_7_1.place(x=680,y=80,anchor=tkinter.W,width=100,height=25)
label_7_1.place(x=760,y=80,anchor=tkinter.W,width=70,height=35)
label_7_1_1.place(x=570,y=80,anchor=tkinter.W,width=70,height=35)

button7_2.place(x=470,y=110,anchor=tkinter.W,width=68,height=40)
entry_7_2.place(x=680,y=110,anchor=tkinter.W,width=100,height=25)
label_7_2.place(x=760,y=110,anchor=tkinter.W,width=70,height=35)
label_7_2_1.place(x=570,y=110,anchor=tkinter.W,width=85,height=35)
##########################
def button12_select():
    global protonation_yn
    if var12.get()==1:
        # label_8.config(foreground="#000000")
        protonation_yn = 1

    elif var12.get()==0:
        protonation_yn = 0
        # label_8.config(foreground="#969696")

var12=tkinter.IntVar(value=1)
label_12 = tkinter.Label(frame2,text="Surface Protonation",font=('Arial', '12'))
button12_1 = tkinter.Checkbutton(frame2,variable=var12,onvalue=1, offvalue=0,command=button12_select )

label_12.place(x=479,y=160,anchor=tkinter.W,width=180,height=40)
button12_1.place(x=445,y=160,anchor=tkinter.W,width=40,height=40)

##########################
def button9_select():
    global separate_layer
    if var9.get()==1:
        separate_layer =1
        print("set ")
        # label_9.config(foreground="#000000")

    elif var9.get()==0:
        separate_layer = 0
        print("not set")
        # label_9.config(foreground="#969696")

var9=tkinter.IntVar(value=0)
label_9 = tkinter.Label(frame2,text="Barrier Layer",font=('Arial', '12'))
button9_1 = tkinter.Checkbutton(frame2,variable=var9,onvalue=1, offvalue=0,command=button9_select )

label_9.place(x=449,y=220,anchor=tkinter.W,width=200,height=40)
button9_1.place(x=445,y=220,anchor=tkinter.W,width=40,height=40)


label_11_1 = tkinter.Label(frame3, font=('Arial', '12'),text="Progress:")
label_11_12 = tkinter.Label(frame3, font=('Arial', '12'),bg="white",text="0/0 (0%)")
label_11_2 = tkinter.Label(frame3, font=('Arial', '12'),text="Status:")
label_11_22 = tkinter.Label(frame3, font=('Arial', '12'),bg="white",text="Setting Job")
label_11_1.place(x=80, y=35,anchor=tkinter.W, width=100, height=35)
label_11_2.place(x=430, y=35,anchor=tkinter.W, width=100, height=35)
label_11_12.place(x=170, y=35,anchor=tkinter.W, width=300, height=35)
label_11_22.place(x=510, y=35,anchor=tkinter.W, width=300, height=35)

###################################
win.mainloop()
