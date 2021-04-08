#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 22:28:23 2021

@author: divyang
"""
import pandas as pd
import numpy as np
import os

#%% Generate data

class Coil:
    
    def __init__(self):
        
        Path = os.getcwd()
        Path = Path[:-4] 				# 4 is the lenght of "CAEB" word

        data = pd.read_csv(Path + '/Inputs/Parameters.txt', sep = '=', header = None)
        data = np.array( data.loc[:, 1] )

        self.Total_Coils = int(data[0])
        self.TurnsInCoil = int(data[1])
        self.MajorRadius = data[2]
        self.TF_OffSet   = data[3]
        self.InnerRadius = data[4]
        self.OuterRadius = data[5]
        self.thickness   = data[6]
        self.I_in_1_turn = data[7]
        
        print('Number of components: ' + str(self.Total_Coils))
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
        J = [self.TurnsInCoil* self.I_in_1_turn/ ((self.OuterRadius - self.InnerRadius)* self.thickness)]
        
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

if __name__ == '__main__':

    Make_Coils = Coil()
    Make_Coils.Write_data()
