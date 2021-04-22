#!/usr/bin/env python3
# -*- coding: utf-8 -*-

print('   -----------------------------------------------------------------------------------   ')
print('                                         Inputs.py                                       ')
print(' ')
print('   Importing data from Input files...')

from numpy import seterr, sin, cos, sqrt, pi, array, linspace, meshgrid, arctan2, deg2rad, savetxt
import sys
from os import getcwd
import xlrd

seterr( divide = 'ignore')
seterr( invalid= 'ignore')

#%%
# IsFloat will check that the data is floating point or not
def IsFloat(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def MatAll(row_number, Nmbr_inCAE):

    ''' MatAll will replace data(x,:) operator
        this funtion is only for CAE_1.xlsx not for CAE_2.xlsx
        row_number is the row address of cell in CAE_1.xlsx file
        start is the starting point which is 4 '''

    stop = Nmbr_inCAE + 3
    start = 4
    AA = []
    for i in range(start - 1, stop):
        AA.append( data.cell_value( row_number - 1, i))

    return AA

def Cells(q):

    ''' Cells will creat the vector of empty cells
        e.g. [[],[],[], ... ] '''

    A = []
    for i in range(q):
        A.append([])

    return A

# Cart2Pol will convert cartesian points into the polar points
def Cart2Pol(x, y, z):
    ρ = sqrt(x**2 + y**2)
    φ = arctan2(y, x)
    return ρ, φ, z

# Pol2Cart will convert polar points into the cartesian points
def Pol2Cart(ρ, φ, z):
    x = ρ* cos(φ)
    y = ρ* sin(φ)
    return x, y, z

def Euler(θ, φ, ψ):
    return array([[ cos(θ)* cos(φ), sin(ψ)* sin(θ)* cos(φ) - cos(ψ)* sin(φ), cos(ψ)* sin(θ)* cos(φ) + sin(ψ)* sin(φ)],
                  [ cos(θ)* sin(φ), sin(ψ)* sin(θ)* sin(φ) + cos(ψ)* cos(φ), cos(ψ)* sin(θ)* sin(φ) - sin(ψ)* cos(φ)],
                  [        -sin(θ),                          sin(ψ)* cos(θ),                          cos(ψ)* cos(θ)] ])

def Mr_Exit(*args):
    print('\033[1;31;40m ERROR!!! \033[0m : ' + ''.join(args) )
    input('Press ENTER to exit')
    sys.exit(1)

#%%
Path = getcwd()
Path = Path[:-14]                   # 14 is the length of folder name "Magnet_modules"

FILE_CAE_1 = xlrd.open_workbook( Path + 'Inputs/CAE_1.xlsx')
FILE_CAE_2 = xlrd.open_workbook( Path + 'Inputs/CAE_2.xlsx')
data = FILE_CAE_1.sheet_by_index(0)
data_2 = FILE_CAE_2.sheet_by_index(0)

µk = data_2.cell_value(19,6)
if (type(µk) == str):
    Mr_Exit('Please check the permeability, 20 row, 7 column in CAE_2.xlsx')
µ0 = 4* pi* 1e-7
µ  = µk* µ0

need_plot = data_2.cell_value(21,6)
if (type(need_plot) == str):
    Mr_Exit('Please check the plotting option, 22 row, 7 column in CAE_2.xlsx')
elif (need_plot != 1) and (need_plot != 0):
    Mr_Exit('Please check the plotting option, 22 row, 7 column in CAE_2.xlsx')

need_curr_dir = data_2.cell_value(25,6)
if (type(need_curr_dir) == str):
    Mr_Exit('Please check the current directon option, 26 row, 7 column in CAE_2.xlsx')
elif (need_curr_dir != 1) and (need_curr_dir != 0):
    Mr_Exit('Please check the current directon option, 26 row, 7 column in CAE_2.xlsx')

fied_length_quiver = data_2.cell_value(27,6)
if (type(fied_length_quiver) == str):
    Mr_Exit('Please check the field length option, 28 row, 7 column in CAE_2.xlsx')


#%% filament

if (type(data.cell_value(2, 3)) == str):
    Mr_Exit('Please check the filament number option, 3 row, 4 column in CAE_1.xlsx')

f_Nmbr = int(data.cell_value(2, 3))

if f_Nmbr > 0:
    print('   Importing Filaments...')

    f_cnt = array([MatAll(9, f_Nmbr), MatAll(10, f_Nmbr), MatAll(11, f_Nmbr)])
    f_l = MatAll(13, f_Nmbr)
    # euler transformation angles
    f_ppye = MatAll(17, f_Nmbr)
    f_ttta = MatAll(18, f_Nmbr)
    f_pphi = MatAll(19, f_Nmbr)
    # current in filament
    f_I = MatAll(21, f_Nmbr)

    if not(isinstance(f_cnt[0, 0], float)):
        for i in range(0, f_Nmbr):
            if not( IsFloat( f_cnt[0, i])):
                Mr_Exit('Please check the center of FILAMENT, 9 row, ', str(i + 4),' column in CAE_1.xlsx')

            if not( IsFloat( f_cnt[1, i])):
                Mr_Exit('Please check the center of FILAMENT, 10 row, ', str(i + 4),' column in CAE_1.xlsx')

            if not( IsFloat( f_cnt[2, i])):
                Mr_Exit('Please check the center of FILAMENT, 11 row, ', str(i + 4),' column in CAE_1.xlsx')

    #
    f_operator = Cells(f_Nmbr)
    f_center   = Cells(f_Nmbr)

    for i in range(0, f_Nmbr):
        if not( IsFloat( f_l[i])) or ( f_l[i] < 0):
            Mr_Exit('Please check the length of FILAMENT, 13 row, ', str(i + 4),' column in CAE_1.xlsx')

        if not( IsFloat( f_pphi[i])):
            Mr_Exit('Please check the Euler angles of FILAMENT, 19 row, ', str(i + 4),' column in CAE_1.xlsx')

        if not( IsFloat( f_ppye[i])):
            Mr_Exit('Please check the Euler angles of FILAMENT, 17 row, ', str(i + 4),' column in CAE_1.xlsx')

        if not( IsFloat( f_ttta[i])):
            Mr_Exit('Please check the Euler angles of FILAMENT, 18 row, ', str(i + 4),' column in CAE_1.xlsx')

        if not( IsFloat( f_I[i])):
            Mr_Exit('Please check the current of FILAMENT, 21 row, ', str(i + 4),' column in CAE_1.xlsx')

        # degree to radian converter
        f_ttta[i] = deg2rad( f_ttta[i] )
        f_pphi[i] = deg2rad( f_pphi[i] )
        f_ppye[i] = deg2rad( f_ppye[i] )

        # separation
        f_operator[i] = Euler(f_ttta[i], f_pphi[i], f_ppye[i])
        f_center[i]   = f_cnt[:,i]

    print('   Done.')

#%% slab
s_Nmbr = int(data.cell_value(25, 3))

if s_Nmbr > 0:
    print('   Importing Slabs...')

    s_cnt = array([MatAll(32, s_Nmbr), MatAll(33, s_Nmbr), MatAll(34, s_Nmbr)])
    s_l = MatAll(36, s_Nmbr)
    s_b = MatAll(37, s_Nmbr)
    s_h = MatAll(38, s_Nmbr)
    # euler transformation angels
    s_ppye = MatAll(42, s_Nmbr)
    s_ttta = MatAll(43, s_Nmbr)
    s_pphi = MatAll(44, s_Nmbr)
    # current density in slab
    s_J = MatAll(46, s_Nmbr)

    if not(isinstance(s_cnt[0, 0], float)):
        for i in range(0, s_Nmbr):
            if not( IsFloat( s_cnt[0, i])):
                Mr_Exit('Please check the center of RECTANGULAR SLAB, 32 row, ', str(i + 4), ' column in CAE_1.xlsx')

            if not( IsFloat( s_cnt[1, i])):
                Mr_Exit('Please check the center of RECTANGULAR SLAB, 33 row, ', str(i + 4), ' column in CAE_1.xlsx')

            if not( IsFloat( s_cnt[2, i])):
                Mr_Exit('Please check the center of RECTANGULAR SLAB, 34 row, ', str(i + 4), ' column in CAE_1.xlsx')

    #
    s_ver_org = Cells(s_Nmbr)
    s_operator = Cells(s_Nmbr)
    s_center   = Cells(s_Nmbr)

    for i in range(0, s_Nmbr):
        if not( IsFloat( s_l[i])) or ( s_l[i] < 0):
            Mr_Exit('Please check the length of RECTANGULAR SLAB, 36 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( s_b[i])) or ( s_b[i] < 0):
            Mr_Exit('Please check the breadth of RECTANGULAR SLAB, 37 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( s_h[i])) or ( s_h[i] < 0):
            Mr_Exit('Please check the height of RECTANGULAR SLAB, 38 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( s_pphi[i])):
            Mr_Exit('Please check the Euler angles of RECTANGULAR SLAB, 44 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( s_ppye[i])):
            Mr_Exit('Please check the Euler angles of RECTANGULAR SLAB, 42 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( s_ttta[i])):
            Mr_Exit('Please check the Euler angles of RECTANGULAR SLAB, 43 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( s_J[i])):
            Mr_Exit('Please check the current density of RECTANGULAR SLAB, 46 row, ', str(i + 4), ' column in CAE_1.xlsx')

        # degree to radian converter
        s_ttta[i] = deg2rad( s_ttta[i])
        s_pphi[i] = deg2rad( s_pphi[i])
        s_ppye[i] = deg2rad( s_ppye[i])

        # separation
        s_operator[i] = Euler(s_ttta[i], s_pphi[i], s_ppye[i])
        s_center[i] = s_cnt[:,i]

        ## position of ideal slab :-
        # means one corner is [0 0 0] and other is [b l h]
        # vertex [0] and [6] are dimetrically opposite points (marked by '#%' sign)
        # '-[b/2;l/2;h/2]' will shift the slab to origin
        s_ver_org[i] = Cells(8)
        s_ver_org[i][0] = array([0,      0,      0]) - array([s_b[i]/2, s_l[i]/2, s_h[i]/2]) #%
        s_ver_org[i][1] = array([s_b[i], 0,      0]) - array([s_b[i]/2, s_l[i]/2, s_h[i]/2])
        s_ver_org[i][2] = array([s_b[i], s_l[i], 0]) - array([s_b[i]/2, s_l[i]/2, s_h[i]/2])
        s_ver_org[i][3] = array([0,      s_l[i], 0]) - array([s_b[i]/2, s_l[i]/2, s_h[i]/2])
        s_ver_org[i][4] = array([0,      0,      s_h[i]]) - array([s_b[i]/2, s_l[i]/2, s_h[i]/2])
        s_ver_org[i][5] = array([s_b[i], 0,      s_h[i]]) - array([s_b[i]/2, s_l[i]/2, s_h[i]/2])
        s_ver_org[i][6] = array([s_b[i], s_l[i], s_h[i]]) - array([s_b[i]/2, s_l[i]/2, s_h[i]/2]) #%
        s_ver_org[i][7] = array([0,      s_l[i], s_h[i]]) - array([s_b[i]/2, s_l[i]/2, s_h[i]/2])

    print('   Done.')

#%% circular arc
a_Nmbr = int(data.cell_value(50, 3))

if a_Nmbr > 0:
    print('   Importing Circular arcs...')

    a_cnt = array([MatAll(57, a_Nmbr), MatAll(58, a_Nmbr), MatAll(59, a_Nmbr)])
    a_a0 = MatAll(61, a_Nmbr)
    a_phy_1 = MatAll(62, a_Nmbr)
    a_phy_2 = MatAll(63, a_Nmbr)
    # euler transformation angels
    a_ppye = MatAll(67, a_Nmbr)
    a_ttta = MatAll(68, a_Nmbr)
    a_pphi = MatAll(69, a_Nmbr)
    # current in circular arc
    a_I = MatAll(71, a_Nmbr)

    if not(isinstance(a_cnt[0, 0], float)):
        for i in range(0, a_Nmbr):
            if not( IsFloat( a_cnt[0, i])):
                Mr_Exit('Please check the center of CIRCULAR ARC, 57 row, ', str(i + 4), ' column in CAE_1.xlsx')

            if not( IsFloat( a_cnt[1, i])):
                Mr_Exit('Please check the center of CIRCULAR ARC, 58 row, ', str(i + 4), ' column in CAE_1.xlsx')

            if not( IsFloat( a_cnt[2, i])):
                Mr_Exit('Please check the center of CIRCULAR ARC, 59 row, ', str(i + 4), ' column in CAE_1.xlsx')

    #
    a_operator = Cells(a_Nmbr)
    a_cnt_cart = Cells(a_Nmbr)

    for i in range(0, a_Nmbr):
        if not( IsFloat( a_a0[i])) or ( a_a0[i] < 0):
            Mr_Exit('Please check the radius of CIRCULAR ARC, 61 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( a_phy_1[i])):
            Mr_Exit('Please check the initial angle of CIRCULAR ARC, 62 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( a_phy_2[i])):
            Mr_Exit('Please check the final angle of CIRCULAR ARC, 63 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if ( a_phy_1[i] >= a_phy_2[i]):
            Mr_Exit('Please check the initial and final angles of CIRCULAR ARC, 62 or 63 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( a_pphi[i])):
            Mr_Exit('Please check the Euler angles of CIRCULAR ARC, 69 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( a_ppye[i])):
            Mr_Exit('Please check the Euler angles of CIRCULAR ARC, 67 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( a_ttta[i])):
            Mr_Exit('Please check the Euler angles of CIRCULAR ARC, 68 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( a_I[i])):
            Mr_Exit('Please check the current of CIRCULAR ARC, 71 row, ', str(i + 4), ' column in CAE_1.xlsx')

        # degree to radian converter
        a_phy_1[i] = deg2rad( a_phy_1[i])
        a_phy_2[i] = deg2rad( a_phy_2[i])
        #
        a_ttta[i] = deg2rad( a_ttta[i])
        a_pphi[i] = deg2rad( a_pphi[i])
        a_ppye[i] = deg2rad( a_ppye[i])

        # separation
        a_operator[i] = Euler(a_ttta[i], a_pphi[i], a_ppye[i])
        a_cnt_cart[i] = a_cnt[:,i]

    print('   Done.')

#%% ractangular circular bar
ra_Nmbr = int(data.cell_value(74, 3))

if ra_Nmbr > 0:
    print('   Importing Circular arcs of rectangular cross-section...')

    ra_cnt = array([MatAll(82, ra_Nmbr), MatAll(83, ra_Nmbr), MatAll(84, ra_Nmbr)])
    ra_r_inner = MatAll(86, ra_Nmbr)
    ra_r_outer = MatAll(87, ra_Nmbr)
    ra_phy_1 = MatAll(88, ra_Nmbr)
    ra_phy_2 = MatAll(89, ra_Nmbr)
    ra_tckns = MatAll(90, ra_Nmbr)
    # euler transformation angels
    ra_ppye = MatAll(94, ra_Nmbr)
    ra_ttta = MatAll(95, ra_Nmbr)
    ra_pphi = MatAll(96, ra_Nmbr)
    # current density in slab
    ra_J = MatAll(98, ra_Nmbr)

    if not(isinstance(ra_cnt[0, 0], float)):
        for i in range(0, ra_Nmbr):
            if not( IsFloat( ra_cnt[0, i])):
                Mr_Exit('Please check the center of RACTANGULAR ARC, 82 row , ', str(i + 4), ' column in CAE_1.xlsx')

            if not( IsFloat( ra_cnt[1, i])):
                Mr_Exit('Please check the center of RACTANGULAR ARC, 83 row, ', str(i + 4), ' column in CAE_1.xlsx')

            if not( IsFloat( ra_cnt[2, i])):
                Mr_Exit('Please check the center of RACTANGULAR ARC, 84 row, ', str(i + 4), ' column in CAE_1.xlsx')

    #
    ra_operator = Cells(ra_Nmbr)
    ra_cnt_cart = Cells(ra_Nmbr)

    for i in range(0, ra_Nmbr):
        if not( IsFloat( ra_r_inner[i])) or ( ra_r_inner[i] < 0):
            Mr_Exit('Please check the inner radius of RACTANGULAR ARC, 86 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( ra_r_outer[i])) or ( ra_r_outer[i] < 0):
            Mr_Exit('Please check the outer radius of RACTANGULAR ARC, 87 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if ra_r_inner[i] >= ra_r_outer[i]:
            Mr_Exit('Please check the inner and outer radius of RACTANGULAR ARC, 86 or 87 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( ra_tckns[i])) or ( ra_tckns[i] <= 0):
            Mr_Exit('Please check the thickness of RACTANGULAR ARC, 90 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( ra_J[i])):
            Mr_Exit('Please check the current density of RACTANGULAR ARC, 98 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( ra_phy_1[i])):
            Mr_Exit('Please check the initial angle of RACTANGULAR ARC, 88 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( ra_phy_2[i])):
            Mr_Exit('Please check the final angle of RACTANGULAR ARC, 89 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if (ra_phy_1[i] >= ra_phy_2[i]):
            Mr_Exit('Please check the initial and final angle of RACTANGULAR ARC, 88 or 89 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( ra_pphi[i])):
            Mr_Exit('Please check the Euler angles of RACTANGULAR ARC, 96 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( ra_ppye[i])):
            Mr_Exit('Please check the Euler angles of RACTANGULAR ARC, 94 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( ra_ttta[i])):
            Mr_Exit('Please check the Euler angles of RACTANGULAR ARC, 95 row, ', str(i + 4), ' column in CAE_1.xlsx')

        # degree to radian converter
        ra_phy_1[i] = deg2rad( ra_phy_1[i])
        ra_phy_2[i] = deg2rad( ra_phy_2[i])
        #
        ra_ttta[i] = deg2rad( ra_ttta[i])
        ra_pphi[i] = deg2rad( ra_pphi[i])
        ra_ppye[i] = deg2rad( ra_ppye[i])

        # separation
        ra_operator[i] = Euler(ra_ttta[i], ra_pphi[i], ra_ppye[i])
        ra_cnt_cart[i] = ra_cnt[:,i]

    print('   Done.')

#%% cylindrical wire
cy_Nmbr = int(data.cell_value(102, 3))

if cy_Nmbr > 0:

    print('   Importing Cylindrical wires...')

    cy_cnt = array([MatAll(109, cy_Nmbr), MatAll(110, cy_Nmbr), MatAll(111, cy_Nmbr)])
    cy_ro = MatAll(113, cy_Nmbr)
    cy_l = MatAll(114, cy_Nmbr)
    # euler transformation angels
    cy_ppye = MatAll(118, cy_Nmbr)
    cy_ttta = MatAll(119, cy_Nmbr)
    cy_pphi = MatAll(120, cy_Nmbr)
    # current in wire
    cy_J = MatAll(122, cy_Nmbr)

    if not(isinstance(cy_cnt[0, 0], float)):
        for i in range(0, cy_Nmbr):
            if not( IsFloat( cy_cnt[0, i])):
                Mr_Exit('Please check the center of CYLINDRICAL WIRE, 109 row, ', str(i + 4), ' column in CAE_1.xlsx')

            if not( IsFloat( cy_cnt[1, i])):
                Mr_Exit('Please check the center of CYLINDRICAL WIRE, 110 row, ', str(i + 4), ' column in CAE_1.xlsx')

            if not( IsFloat( cy_cnt[2, i])):
                Mr_Exit('Please check the center of CYLINDRICAL WIRE, 111 row, ', str(i + 4), ' column in CAE_1.xlsx')
    #
    cy_operator = Cells(cy_Nmbr)
    cy_center   = Cells(cy_Nmbr)

    for i in range(0, cy_Nmbr):
        if not( IsFloat( cy_ro[i])) or ( cy_ro[i] < 0):
            Mr_Exit('Please check the radius of CYLINDRICAL WIRE, 113 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( cy_l[i])) or ( cy_l[i] < 0):
            Mr_Exit('Please check the length of CYLINDRICAL WIRE, 114 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( cy_pphi[i])):
            Mr_Exit('Please check the Euler angles of CYLINDRICAL WIRE, 120 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( cy_ppye[i])):
            Mr_Exit('Please check the Euler angles of CYLINDRICAL WIRE, 118 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( cy_ttta[i])):
            Mr_Exit('Please check the Euler angles of CYLINDRICAL WIRE, 119 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( cy_J[i])):
            Mr_Exit('Please check the current density of CYLINDRICAL WIRE, 122 row, ', str(i + 4), ' column in CAE_1.xlsx')

        # degree to radian converter
        cy_ttta[i] = deg2rad( cy_ttta[i])
        cy_pphi[i] = deg2rad( cy_pphi[i])
        cy_ppye[i] = deg2rad( cy_ppye[i])

        # separation
        cy_operator[i] = Euler(cy_ttta[i], cy_pphi[i], cy_ppye[i])
        cy_center[i] = cy_cnt[:,i]

    print('   Done.')

#%%
Hl_Nmbr = int(data.cell_value(126, 3))

if Hl_Nmbr > 0:

    print('   Importing Helical wires...')

    Hl_cnt = array([MatAll(131, Hl_Nmbr), MatAll(132, Hl_Nmbr), MatAll(133, Hl_Nmbr)])
    Hl_a   = MatAll(135, Hl_Nmbr)
    Hl_L   = MatAll(136, Hl_Nmbr)
    Hl_p   = MatAll(137, Hl_Nmbr)
    # euler transformation angels
    Hl_ppye = MatAll(141, Hl_Nmbr)
    Hl_ttta = MatAll(142, Hl_Nmbr)
    Hl_pphi = MatAll(143, Hl_Nmbr)
    # current in wire
    Hl_I    = MatAll(145, Hl_Nmbr)

    if not(isinstance(Hl_cnt[0, 0], float)):
        for i in range(0, Hl_Nmbr):
            if not( IsFloat( Hl_cnt[0, i])):
                Mr_Exit('Please check the center of HELICAL WIRE, 131 row, ', str(i + 4), ' column in CAE_1.xlsx')

            if not( IsFloat( Hl_cnt[1, i])):
                Mr_Exit('Please check the center of HELICAL WIRE, 132 row, ', str(i + 4), ' column in CAE_1.xlsx')

            if not( IsFloat( Hl_cnt[2, i])):
                Mr_Exit('Please check the center of HELICAL WIRE, 133 row, ', str(i + 4), ' column in CAE_1.xlsx')
    #
    Hl_operator = Cells(Hl_Nmbr)
    Hl_center   = Cells(Hl_Nmbr)

    for i in range(0, Hl_Nmbr):
        if not( IsFloat( Hl_a[i])) or ( Hl_a[i] < 0):
            Mr_Exit('Please check the radius of HELICAL WIRE, 135 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( Hl_L[i])) or ( Hl_L[i] < 0):
            Mr_Exit('Please check the length of HELICAL WIRE, 136 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( Hl_p[i])) or ( Hl_p[i] < 0):
            Mr_Exit('Please check the pitch of HELICAL WIRE, 137 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( Hl_ppye[i])):
            Mr_Exit('Please check the Euler angles of HELICAL WIRE, 141 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( Hl_ttta[i])):
            Mr_Exit('Please check the Euler angles of HELICAL WIRE, 142 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( Hl_pphi[i])):
            Mr_Exit('Please check the Euler angles of HELICAL WIRE, 143 row, ', str(i + 4), ' column in CAE_1.xlsx')

        if not( IsFloat( Hl_I[i])):
            Mr_Exit('Please check the current density of HELICAL WIRE, 145 row, ', str(i + 4), ' column in CAE_1.xlsx')

        # degree to radian converter
        Hl_ttta[i] = deg2rad( Hl_ttta[i])
        Hl_pphi[i] = deg2rad( Hl_pphi[i])
        Hl_ppye[i] = deg2rad( Hl_ppye[i])

        # separation
        Hl_operator[i] = Euler(Hl_ttta[i], Hl_pphi[i], Hl_ppye[i])
        Hl_center[i] = Hl_cnt[:,i]

    print('   Done.')

#%%

if ((f_Nmbr == 0) and (s_Nmbr == 0) and (a_Nmbr == 0) and (ra_Nmbr == 0) and (cy_Nmbr == 0) and (Hl_Nmbr == 0) ):
    Mr_Exit('Please enter the data for atleast an electromagnet, in CAE_1.xlsx')

#%% Points at which magnetic field to be calculated :-

if int( data_2.cell_value(4,3)) == 1:
    cart_pnt = int(1)
    polar_pnt = int(0)
elif int( data_2.cell_value(4,3)) == 0:
    cart_pnt = int(0)
    polar_pnt = int(1)

print('   Imorting space point data...')

if (cart_pnt == 1):
    if ( not( IsFloat( data_2.cell_value( 7,3))) or not( IsFloat( data_2.cell_value( 8,3))) or not( IsFloat( data_2.cell_value( 7,4))) or not( IsFloat( data_2.cell_value( 8,4))) or not( IsFloat( data_2.cell_value( 7,5))) or not( IsFloat( data_2.cell_value( 8,5))) or not( IsFloat( data_2.cell_value( 10,3))) or not( IsFloat( data_2.cell_value( 10,4))) or not( IsFloat( data_2.cell_value( 10,5))) ):
        Mr_Exit('Please check the space data in CAE_2.xlsx')

    if ( ( data_2.cell_value(8,3) > data_2.cell_value(7,3)) or ( data_2.cell_value(8,4) > data_2.cell_value(7,4)) or ( data_2.cell_value(8,5) > data_2.cell_value(7,5)) ):
        Mr_Exit('Please enter the space data with min and max correctly in CAE_2.xlsx')

    if ( ( data_2.cell_value(10,3) <= 0) or ( data_2.cell_value(10,4) <= 0) or ( data_2.cell_value(10,5) <= 0)):
        Mr_Exit('Please enter the space data correctly in CAE_2.xlsx')

    xx = linspace( data_2.cell_value(8,3), data_2.cell_value(7,3), int(data_2.cell_value(10,3)) )
    yy = linspace( data_2.cell_value(8,4), data_2.cell_value(7,4), int(data_2.cell_value(10,4)) )
    zz = linspace( data_2.cell_value(8,5), data_2.cell_value(7,5), int(data_2.cell_value(10,5)) )
    #
    xx, yy, zz = meshgrid(xx, yy, zz)

if (polar_pnt ==1):
    if ( not( IsFloat( data_2.cell_value( 7,3))) or not( IsFloat( data_2.cell_value( 8,3))) or not( IsFloat( data_2.cell_value( 7,4))) or not( IsFloat( data_2.cell_value( 8,4))) or not( IsFloat( data_2.cell_value( 7,5))) or not( IsFloat( data_2.cell_value( 8,5))) or not( IsFloat( data_2.cell_value( 10,3))) or not( IsFloat( data_2.cell_value( 10,4))) or not( IsFloat( data_2.cell_value( 10,5))) ):
        Mr_Exit('Please check the space data in CAE_2.xlsx')

    if ( ( data_2.cell_value(8,3) > data_2.cell_value(7,3)) or ( data_2.cell_value(8,4) > data_2.cell_value(7,4)) or ( data_2.cell_value(8,5) > data_2.cell_value(7,5)) ):
        Mr_Exit('Please enter the space data with min and max correctly in CAE_2.xlsx')

    if ( ( data_2.cell_value(10,3) <= 0) or ( data_2.cell_value(10,4) <= 0) or ( data_2.cell_value(10,5) <= 0)):
        Mr_Exit('Please enter the space data correctly in CAE_2.xlsx')

    # degree to radian converter
    degdata1 = deg2rad(data_2.cell_value(8,4))
    degdata2 = deg2rad(data_2.cell_value(7,4))
    degpoint = int(data_2.cell_value(10,4))
    #
    rr     = linspace( data_2.cell_value(8,3), data_2.cell_value(7,3), int(data_2.cell_value(10,3)) )
    phyphy = linspace( degdata1, degdata2, degpoint)
    zz     = linspace( data_2.cell_value(8,5), data_2.cell_value(7,5), int(data_2.cell_value(10,5)) )
    #
    rr, phyphy, zz = meshgrid(rr,phyphy,zz)
    xx,   yy,   zz = Pol2Cart(rr, phyphy, zz)

print('   Done.')

#%% Contour Plot option
need_contour_plot = data_2.cell_value(30,6)

if (type(need_contour_plot) == str):
    Mr_Exit('Please check the contour plot option, 31 row, 7 column in CAE_2.xlsx')
elif (need_contour_plot != 1.0) and (need_contour_plot != 0.0):
    Mr_Exit('Please check the contour plot option, 31 row, 7 column in CAE_2.xlsx')

if need_contour_plot == 1:
    print('   Checking conditions for contour plots...')

    if cart_pnt == 1:
        logic_x = ( data_2.cell_value(8,3) == data_2.cell_value(7,3) )
        logic_y = ( data_2.cell_value(8,4) == data_2.cell_value(7,4) )
        logic_z = ( data_2.cell_value(8,5) == data_2.cell_value(7,5) )
        if ( (logic_x and logic_y) or (logic_y and logic_z) or (logic_x and logic_z) or not(logic_x or logic_y or logic_z) ):
            Mr_Exit('To draw a contour plot one parameter must be constant, 8 and 9 row in CAE_2.xlsx')

    elif polar_pnt == 1:
        logic_r = ( data_2.cell_value(8,3) == data_2.cell_value(7,3) )
        logic_t = ( data_2.cell_value(8,4) == data_2.cell_value(7,4) )
        logic_z = ( data_2.cell_value(8,5) == data_2.cell_value(7,5) )
        if not( not(logic_r) and not(logic_r) and logic_z):
            Mr_Exit('To draw a contour plot in cylindrical system z must be constant, 8 and 9 row in CAE_2.xlsx')

if need_contour_plot == 1:
    if cart_pnt == 1:

        if logic_x == True:        # Here the x is constant
            cont_xx = linspace( data_2.cell_value(8,4), data_2.cell_value(7,4), int(data_2.cell_value(10,4)) ) # y
            cont_yy = linspace( data_2.cell_value(8,5), data_2.cell_value(7,5), int(data_2.cell_value(10,5)) ) # z
            cont_x, cont_y = meshgrid(cont_xx, cont_yy)
            #
            cont_pp, cont_qq = cont_x.shape

        elif logic_y == True:      # Here the y is constant
            cont_xx = linspace( data_2.cell_value(8,5), data_2.cell_value(7,5), int(data_2.cell_value(10,5)) ) # z
            cont_yy = linspace( data_2.cell_value(8,3), data_2.cell_value(7,3), int(data_2.cell_value(10,3)) ) # x
            cont_x, cont_y = meshgrid(cont_xx, cont_yy)
            #
            cont_pp, cont_qq = cont_x.shape

        elif logic_z == True:      # Here the z is constant
            cont_xx = linspace( data_2.cell_value(8,3), data_2.cell_value(7,3), int(data_2.cell_value(10,3)) ) # x
            cont_yy = linspace( data_2.cell_value(8,4), data_2.cell_value(7,4), int(data_2.cell_value(10,4)) ) # y
            cont_x, cont_y = meshgrid(cont_xx, cont_yy)
            #
            cont_pp, cont_qq = cont_x.shape

    elif polar_pnt == 1:
        if logic_z == True:        # Here the z is constant
            cont_RR = linspace( data_2.cell_value(8,3), data_2.cell_value(7,3), int(data_2.cell_value(10,3)) )      # Same as 'rr'
            cont_TT = linspace( degdata1, degdata2, degpoint)                                                       # Same as 'phyphy'
            cont_RR, cont_TT = meshgrid(cont_RR, cont_TT)
            #
            cont_x = cont_RR * cos(cont_TT)
            cont_y = cont_RR * sin(cont_TT)
            #
            cont_pp, cont_qq = cont_x.shape

    print('   Done.')

q, w, e = xx.shape
s = q* w* e
point_1 = Cells(s)
point_2 = Cells(s)
point_3 = Cells(s)

counter = 0
for i in range(0,q):
    for j in range(0,w):
        for k in range(0,e):
            point_1[counter] = xx[i, j, k]
            point_2[counter] = yy[i, j, k]
            point_3[counter] = zz[i, j, k]
            counter = counter + 1

point = array([point_1,point_2,point_3])

#%% Export data

file = open(r'Bridge/ImportantVar.txt', 'w')
file.write( str(µ) + '\n' + str(s) + '\n' + str(f_Nmbr) + '\n' + str(s_Nmbr) + '\n' + str(a_Nmbr) \
            + '\n' + str(ra_Nmbr) + '\n' + str(cy_Nmbr) + '\n' + str(Hl_Nmbr) )
file.close()

if f_Nmbr > 0:
    
    file1 = open(r'Bridge/f_center.txt', 'w')
    file2 = open(r'Bridge/f_operator.txt', 'w')
    file3 = open(r'Bridge/f_l.txt', 'w')
    file4 = open(r'Bridge/f_I.txt', 'w')
    
    for i in range(0, f_Nmbr):
        
        savetxt(file1, f_center[i])
        savetxt(file2, f_operator[i])
        
        file3.write( str( f_l[i] ) + '\n' )
        file4.write( str( f_I[i] ) + '\n' )
    
    file1.close()
    file2.close()
    file3.close()
    file4.close()

if s_Nmbr > 0:
    
    file1 = open(r'Bridge/s_center.txt', 'w')
    file2 = open(r'Bridge/s_operator.txt', 'w')
    file3 = open(r'Bridge/s_J.txt', 'w')
    file4 = open(r'Bridge/s_ver_org_0.txt', 'w')
    file5 = open(r'Bridge/s_ver_org_6.txt', 'w')

    for i in range(0, s_Nmbr):
        
        savetxt(file1, s_center[i])
        savetxt(file2, s_operator[i])

        file3.write( str( s_J[i] ) + '\n' )

        savetxt(file4, s_ver_org[i][0])
        savetxt(file5, s_ver_org[i][6])

    file1.close()
    file2.close()
    file3.close()
    file4.close()
    file5.close()

if a_Nmbr > 0:

    file1 = open(r'Bridge/a_cnt_cart.txt', 'w')
    file2 = open(r'Bridge/a_operator.txt', 'w')
    file3 = open(r'Bridge/a_I.txt', 'w')
    file4 = open(r'Bridge/a_phy_1.txt', 'w')
    file5 = open(r'Bridge/a_phy_2.txt', 'w')
    file6 = open(r'Bridge/a_a0.txt', 'w')

    for i in range(0, a_Nmbr):
        
        savetxt(file1, a_cnt_cart[i])
        savetxt(file2, a_operator[i])

        file3.write( str( a_I[i] ) + '\n' )
        file4.write( str( a_phy_1[i] ) + '\n' )
        file5.write( str( a_phy_2[i] ) + '\n' )
        file6.write( str( a_a0[i] ) + '\n' )

    file1.close()
    file2.close()
    file3.close()
    file4.close()
    file5.close()
    file6.close()

if ra_Nmbr > 0:
    
    file1 = open(r'Bridge/ra_cnt_cart.txt', 'w')
    file2 = open(r'Bridge/ra_r_inner.txt', 'w')
    file3 = open(r'Bridge/ra_r_outer.txt', 'w')
    file4 = open(r'Bridge/ra_phy_1.txt', 'w')
    file5 = open(r'Bridge/ra_phy_2.txt', 'w')
    file6 = open(r'Bridge/ra_tckns.txt', 'w')
    file7 = open(r'Bridge/ra_operator.txt', 'w')
    file8 = open(r'Bridge/ra_J.txt', 'w')
    
    for i in range(0, ra_Nmbr):
        
        savetxt(file1, ra_cnt_cart[i])
        
        file2.write( str( ra_r_inner[i] ) + '\n' )
        file3.write( str( ra_r_outer[i] ) + '\n' )
        file4.write( str( ra_phy_1[i]   ) + '\n' )
        file5.write( str( ra_phy_2[i]   ) + '\n' )
        file6.write( str( ra_tckns[i]   ) + '\n' )
        
        savetxt(file7, ra_operator[i])
        
        file8.write( str( ra_J[i]       ) + '\n' )
    
    file1.close()
    file2.close()
    file3.close()
    file4.close()
    file5.close()
    file6.close()
    file7.close()
    file8.close()

if cy_Nmbr > 0:

    file1 = open(r'Bridge/cy_center.txt', 'w')
    file2 = open(r'Bridge/cy_l.txt', 'w')
    file3 = open(r'Bridge/cy_ro.txt', 'w')
    file4 = open(r'Bridge/cy_operator.txt', 'w')
    file5 = open(r'Bridge/cy_J.txt', 'w')

    for i in range(0, cy_Nmbr):
        
        savetxt(file1, cy_center[i])
        
        file2.write( str( cy_l[i] )  + '\n' )
        file3.write( str( cy_ro[i] ) + '\n' )
        
        savetxt(file4, cy_operator[i])
        
        file5.write( str( cy_J[i]  ) + '\n' )

    file1.close()
    file2.close()
    file3.close()
    file4.close()
    file5.close()

file = open(r'Bridge/points.txt', 'w')
savetxt(file, point)
file.close()

#%%

print('   All the input parameters are validated...')
print('   -----------------------------------------------------------------------------------   ')
