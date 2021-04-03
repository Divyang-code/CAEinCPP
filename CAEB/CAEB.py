#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 22:28:23 2021

@author: divyang
"""
import pandas as pd
import numpy as np

#%% Generate data

class Coil:
    
    def __init__(self, Total_Coils, TurnsInCoil, MajorRadius, TF_OffSet, InnerRadius, OuterRadius, thickness, I_in_1_turn):
        
        self.Total_Coils = Total_Coils
        self.TurnsInCoil = TurnsInCoil

        self.MajorRadius = MajorRadius
        self.TF_OffSet   = TF_OffSet

        self.InnerRadius = InnerRadius
        self.OuterRadius = OuterRadius
        self.thickness   = thickness

        self.I_in_1_turn = I_in_1_turn
        
        print('Number of components: ' + str(Total_Coils))
        print('Copy this number and Builder.xlsx data into CAE_1.xlsx')
    
    def Data_ring(self):
        
        # centre
        xc = [self.MajorRadius + self.TF_OffSet]
        yc = [0]
        zc = [0]

        # dimensions
        r_in  = [self.InnerRadius]
        r_out = [self.OuterRadius]
        ang_i = [0]
        ang_f = [360]
        Thkns = [self.thickness]
        
        # eular rotaion angles
        Rx = [90]
        Ry = [0]
        Rz = [0]
        
        # current density
        J = [self.TurnsInCoil* self.I_in_1_turn/ (self.OuterRadius - self.InnerRadius)]
        
        return xc, yc, zc, r_in, r_out, ang_i, ang_f, Thkns, Rx, Ry, Rz, J

    def Coil_System(self):
        
        xc, yc, zc, r_in, r_out, ang_i, ang_f, Thkns, Rx, Ry, Rz, J = self.Data_ring()

        # Unchanged variables
        zc    = zc*    self.Total_Coils
        r_in  = r_in*  self.Total_Coils
        r_out = r_out* self.Total_Coils
        ang_i = ang_i* self.Total_Coils
        ang_f = ang_f* self.Total_Coils
        Thkns = Thkns* self.Total_Coils
        Rx    = Rx*    self.Total_Coils
        Ry    = Ry*    self.Total_Coils
        J     = J*     self.Total_Coils

        ang = np.arange(0, 360, 360/ self.Total_Coils)
        Rz = [i for i in ang]
    
        # x and y
        r = xc[0]
        xc, yc = [], []
        for phi in Rz:
            xc.append(r* np.cos(phi* np.pi/180))
            yc.append(r* np.sin(phi* np.pi/180))

        return xc, yc, zc, r_in, r_out, ang_i, ang_f, Thkns, Rx, Ry, Rz, J
    
    def Write_data(self):
        
        xc, yc, zc, r_in, r_out, ang_i, ang_f, Thkns, Rx, Ry, Rz, J = self.Coil_System()
        
        Out = pd.DataFrame({'x': (['X center'] + xc),
                            'y': (['Y center'] + yc),
                            'z': (['Z center'] + zc),
                            'empty1': np.nan,
                            'r_in': (['Inner radius']  + r_in),
                            'r_out': (['Outer radius'] + r_out),
                            'ang_i': (['Initial angle']  + ang_i),
                            'ang_f': (['Final angle']  + ang_f),
                            'Thkns': (['Thickness']  + Thkns),
                            'empty2': np.nan,
                            'pye': (['X rotation'] + Rx),
                            'tta': (['Y rotation'] + Ry),
                            'phi': (['Z rotation'] + Rz),
                            'empty3': np.nan,
                            'J': (['Current density'] + J)},

                            columns = ['x','y','z',
                                       'empty1',
                                       'r_in','r_out','ang_i','ang_f','Thkns',
                                       'empty2',
                                       'pye','tta','phi',
                                       'empty3',
                                       'J'
                                      ]
                          )
        
        Out = Out.rename(columns={'empty1':'', 'empty2':'','empty3':''})
        Out = Out.transpose()

        writer = pd.ExcelWriter('Builder.xlsx', engine='xlsxwriter')
        Out.to_excel(writer, index = False, sheet_name = 'CirRact')
        writer.save()

#%% Inputs

Total_Coils = 6
TurnsInCoil = 28
MajorRadius = 30e-2
TF_OffSet   = 2e-2
InnerRadius = 12e-2
OuterRadius = 15.5e-2
thickness   = 5e-2
I_in_1_turn = 10

Make_Coils = Coil(Total_Coils, TurnsInCoil, MajorRadius, TF_OffSet, InnerRadius, OuterRadius, thickness, I_in_1_turn)
Make_Coils.Write_data()
