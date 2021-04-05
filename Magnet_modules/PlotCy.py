#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from numpy import zeros, ones, arange, dot, array, linspace, meshgrid, \
                  cos, sin, pi, isnan, sqrt
from Inputs import f_Nmbr, s_Nmbr, a_Nmbr, ra_Nmbr, cy_Nmbr, Hl_Nmbr, Cart2Pol, \
                   need_plot, need_contour_plot, point, need_curr_dir, \
                   fied_length_quiver, s, Cells

def Pol2Cart(rho, phi, z):
    x = rho * cos(phi)
    y = rho * sin(phi)
    return( x, y, z)

# Filaments
def draw_line(ax, center, operator, l, I, need_curr_dir):
    
    fil_1  = dot( operator, array([0,0, l/2]).reshape([3,1], order='F') ) + center.reshape([3,1], order='F')
    fil_2  = dot( operator, array([0,0,-l/2]).reshape([3,1], order='F') ) + center.reshape([3,1], order='F')

    x = array([float( fil_1[0]), float( fil_2[0]) ])
    y = array([float( fil_1[1]), float( fil_2[1]) ])
    z = array([float( fil_1[2]), float( fil_2[2]) ])
    
    ax.plot( x, y, z, color = 'fuchsia', lw = 3)
    
    if need_curr_dir == 1:
        if I > 0:
            ax.quiver(center[0], center[1], center[2], float( fil_1[0] - center[0] ), float( fil_1[1] - center[1] ), float( fil_1[2] - center[2] ) , lw = 3, Color = 'yellow', normalize = True)
        elif I < 0:
            ax.quiver(center[0], center[1], center[2], float( fil_2[0] - center[0] ), float( fil_2[1] - center[1] ), float( fil_2[2] - center[2] ) , lw = 3, Color = 'yellow', normalize = True)

# Slabs
def draw_slab(ax, mp3d, center, operator, ver_org, J, need_curr_dir):
    
    # transformation of slab
    vertices_trnf = Cells(8)
    for j in range(0,8):
        vertices_trnf[j] = dot(operator, ver_org[j].reshape(3,1) ) + center.reshape([3,1])

    s_verts1 = [( float(vertices_trnf[0][0]), float(vertices_trnf[0][1]), float(vertices_trnf[0][2])),
                ( float(vertices_trnf[1][0]), float(vertices_trnf[1][1]), float(vertices_trnf[1][2])),
                ( float(vertices_trnf[2][0]), float(vertices_trnf[2][1]), float(vertices_trnf[2][2])),
                ( float(vertices_trnf[3][0]), float(vertices_trnf[3][1]), float(vertices_trnf[3][2])) ]

    s_verts2 = [( float(vertices_trnf[0][0]), float(vertices_trnf[0][1]), float(vertices_trnf[0][2])),
                ( float(vertices_trnf[4][0]), float(vertices_trnf[4][1]), float(vertices_trnf[4][2])),
                ( float(vertices_trnf[5][0]), float(vertices_trnf[5][1]), float(vertices_trnf[5][2])),
                ( float(vertices_trnf[1][0]), float(vertices_trnf[1][1]), float(vertices_trnf[1][2])) ]

    s_verts3 = [( float(vertices_trnf[0][0]), float(vertices_trnf[0][1]), float(vertices_trnf[0][2])),
                ( float(vertices_trnf[3][0]), float(vertices_trnf[3][1]), float(vertices_trnf[3][2])),
                ( float(vertices_trnf[7][0]), float(vertices_trnf[7][1]), float(vertices_trnf[7][2])),
                ( float(vertices_trnf[4][0]), float(vertices_trnf[4][1]), float(vertices_trnf[4][2])) ]

    s_verts4 = [( float(vertices_trnf[3][0]), float(vertices_trnf[3][1]), float(vertices_trnf[3][2])),
                ( float(vertices_trnf[7][0]), float(vertices_trnf[7][1]), float(vertices_trnf[7][2])),
                ( float(vertices_trnf[6][0]), float(vertices_trnf[6][1]), float(vertices_trnf[6][2])),
                ( float(vertices_trnf[2][0]), float(vertices_trnf[2][1]), float(vertices_trnf[2][2])) ]

    s_verts5 = [( float(vertices_trnf[1][0]), float(vertices_trnf[1][1]), float(vertices_trnf[1][2])),
                ( float(vertices_trnf[2][0]), float(vertices_trnf[2][1]), float(vertices_trnf[2][2])),
                ( float(vertices_trnf[6][0]), float(vertices_trnf[6][1]), float(vertices_trnf[6][2])),
                ( float(vertices_trnf[5][0]), float(vertices_trnf[5][1]), float(vertices_trnf[5][2])) ]

    s_verts6 = [( float(vertices_trnf[4][0]), float(vertices_trnf[4][1]), float(vertices_trnf[4][2])),
                ( float(vertices_trnf[7][0]), float(vertices_trnf[7][1]), float(vertices_trnf[7][2])),
                ( float(vertices_trnf[6][0]), float(vertices_trnf[6][1]), float(vertices_trnf[6][2])),
                ( float(vertices_trnf[5][0]), float(vertices_trnf[5][1]), float(vertices_trnf[5][2])) ]
    ###
    f1 = mp3d.art3d.Poly3DCollection([ s_verts1 ], alpha=0.5, linewidth=0.7, edgecolors='k' )
    f2 = mp3d.art3d.Poly3DCollection([ s_verts2 ], alpha=0.5, linewidth=0.7, edgecolors='k' )
    f3 = mp3d.art3d.Poly3DCollection([ s_verts3 ], alpha=0.5, linewidth=0.7, edgecolors='k' )
    f4 = mp3d.art3d.Poly3DCollection([ s_verts4 ], alpha=0.5, linewidth=0.7, edgecolors='k' )
    f5 = mp3d.art3d.Poly3DCollection([ s_verts5 ], alpha=0.5, linewidth=0.7, edgecolors='k' )
    f6 = mp3d.art3d.Poly3DCollection([ s_verts6 ], alpha=0.5, linewidth=0.7, edgecolors='k' )
    #
    f1.set_facecolor((0.1, 0.5, 1, 0.5)) # Last argument is alpha = 0.5
    f2.set_facecolor((0.1, 0.5, 1, 0.5))
    f3.set_facecolor((0.1, 0.5, 1, 0.5))
    f4.set_facecolor((0.1, 0.5, 1, 0.5))
    f5.set_facecolor((0.1, 0.5, 1, 0.5))
    f6.set_facecolor((0.1, 0.5, 1, 0.5))
    #
    ax.add_collection3d(f1)
    ax.add_collection3d(f2)
    ax.add_collection3d(f3)
    ax.add_collection3d(f4)
    ax.add_collection3d(f5)
    ax.add_collection3d(f6)
    #
    if need_curr_dir == 1:
        if J > 0:
            ax.quiver( vertices_trnf[0][0], vertices_trnf[0][1], vertices_trnf[0][2], vertices_trnf[4][0] - vertices_trnf[0][0], vertices_trnf[4][1] - vertices_trnf[0][1], vertices_trnf[4][2] - vertices_trnf[0][2], normalize = True, color = 'yellow' )
        elif J < 0:
            ax.quiver( vertices_trnf[4][0], vertices_trnf[4][1], vertices_trnf[4][2], vertices_trnf[0][0] - vertices_trnf[4][0], vertices_trnf[0][1] - vertices_trnf[4][1], vertices_trnf[0][2] - vertices_trnf[4][2], normalize = True, color = 'yellow' )

# Arcs
def circle(ax, r, theta_S, theta_E, center, Euler):
    ''' r       ---> radius of the circle
        theta_S ---> starting angle of the circle
        theta_E ---> ending angle of the circle
        center  ---> center of the circle (shape 3*1)
        Euler   ---> Euler rotation matrix of the circle '''
    # Circle at origin
    t = arange( theta_S, theta_E + 0.05, 0.05)
    x, y, z = Pol2Cart( r, t, 0)
    z = zeros( len(x) )
    C = array([x, y, z])

    # Euler rotation of the circle
    cir = dot( Euler, C)
    
    # center shifting
    circle = cir + center
    X = circle[0,:]
    Y = circle[1,:]
    Z = circle[2,:]
    
    ax.plot(X, Y, Z)
     
# roll to draw the cylindre
def roll(ax, R, zi, zf, ang_i, ang_f, center, Euler):
    
    t = linspace( ang_i, ang_f, 50)
    z = array([zi, zf])
    t, z = meshgrid(t, z)
    p, q = t.shape
    r = R* ones([p,q], float)
    
    # cylindrical coordinates to Cartesian coordinate
    x, y, z = Pol2Cart(r,t,z)

    # Euler rotation
    rot = dot( Euler, array([ x.ravel(), y.ravel(), z.ravel()]) )

    x_rot = rot[0,:].reshape( x.shape)
    y_rot = rot[1,:].reshape( y.shape)
    z_rot = rot[2,:].reshape( z.shape)

    # transformation
    x = x_rot + center[0]
    y = y_rot + center[1]
    z = z_rot + center[2]
    
    ax.plot_surface( x, y, z, alpha = 0.5, color = (0.1, 0.5, 1), linewidth=0.7, edgecolors='k' )

# Disk part in cylindre
def disk(ax, Ri, Rf, ang_i, ang_f, center, Euler, d):
    
    r = array( [ Ri, Rf] )
    t = linspace( ang_i, ang_f, 50)
    t, r = meshgrid(t, r)
    p, q = t.shape
    z = zeros([p,q], float) + d
    
    # cylindrical coordinates to Cartesian coordinate
    x, y, z = Pol2Cart(r,t,z)

    # Euler rotation
    rot = dot( Euler, array([ x.ravel(), y.ravel(), z.ravel()]) )

    x_rot = rot[0,:].reshape( x.shape)
    y_rot = rot[1,:].reshape( y.shape)
    z_rot = rot[2,:].reshape( z.shape)

    # transformation
    x = x_rot + center[0]
    y = y_rot + center[1]
    z = z_rot + center[2]
    
    ax.plot_surface( x, y, z, alpha = 0.5, color = (0.1, 0.5, 1), linewidth=0.7, edgecolors='k' )

#%% Draw electro-magnets

if __name__=='__main__':

    # Read magnetic field
    import pandas as pd

    df = pd.read_csv('MAG_DATA.csv', sep=',', lineterminator='\n')
    
    Bx = array(df['Bx (T)'])
    By = array(df['By (T)'])
    Bz = array(df['Bz (T)'])
    
    if need_plot == 1:
        
        print('   Drawing Electro-Magnets and field vectors...')
        
        import matplotlib.pyplot as plt
        import mpl_toolkits.mplot3d as mp3d
        
        fig = plt.figure('Electromagnet')
        ax = fig.add_subplot( 111 , projection='3d')
        
        # Filament
        if f_Nmbr > 0:

            from Inputs import f_center, f_operator, f_l, f_I
            
            for i in range(0, f_Nmbr): 
                draw_line(ax, f_center[i], f_operator[i], f_l[i], f_I[i], need_curr_dir)

        # Slab
        if s_Nmbr > 0:
            
            from Inputs import s_center, s_operator, s_ver_org, s_J

            for i in range(0, s_Nmbr):
                draw_slab(ax, mp3d, s_center[i], s_operator[i], s_ver_org[i], s_J[i], need_curr_dir)

        # Circular Arc
        if a_Nmbr > 0:

            from Inputs import a_a0, a_phy_1, a_phy_2, a_cnt_cart, a_operator
            
            for i in range(0, a_Nmbr):
                circle(ax, a_a0[i], a_phy_1[i], a_phy_2[i], a_cnt_cart[i].reshape(3,1), a_operator[i])
        
        # Ractangular arcs
        if ra_Nmbr > 0:
            
            from Inputs import ra_r_inner, ra_r_outer, ra_operator, ra_tckns, ra_phy_1, ra_phy_2, ra_cnt_cart

            for i in range(0, ra_Nmbr):
                # roll
                roll(ax, ra_r_inner[i], -ra_tckns[i]/2, ra_tckns[i]/2, ra_phy_1[i], ra_phy_2[i], ra_cnt_cart[i], ra_operator[i])
                roll(ax, ra_r_outer[i], -ra_tckns[i]/2, ra_tckns[i]/2, ra_phy_1[i], ra_phy_2[i], ra_cnt_cart[i], ra_operator[i])
                # disk
                disk(ax, ra_r_inner[i], ra_r_outer[i], ra_phy_1[i], ra_phy_2[i], ra_cnt_cart[i], ra_operator[i], ra_tckns[i]/2)
                disk(ax, ra_r_inner[i], ra_r_outer[i], ra_phy_1[i], ra_phy_2[i], ra_cnt_cart[i], ra_operator[i],-ra_tckns[i]/2)

        # Cylindrical wires
        if cy_Nmbr > 0:

            from Inputs import cy_ro, cy_center, cy_operator, cy_l
            
            for i in range(0, cy_Nmbr):
                disk(ax, 0, cy_ro[i], 0, 2* pi, cy_center[i], cy_operator[i], cy_l[i]/2 )
                disk(ax, 0, cy_ro[i], 0, 2* pi, cy_center[i], cy_operator[i],-cy_l[i]/2 )
                roll(ax, cy_ro[i],-cy_l[i]/2, cy_l[i]/2, 0, 2* pi, cy_center[i], cy_operator[i])

        B = array([Bx, By, Bz])

        # Sketch the field vectors
        length_B_vec = Cells(s)
        for j in range(0, s):
            if ( isnan( B[0,j]) or isnan( B[1,j]) or isnan( B[1,j]) ):
                length_B_vec[j] = 1
            else:
                length_B_vec[j] = sqrt( B[0,j]**2 + B[1,j]**2 + B[2,j]**2 )
        
        length_B_vec_maxx = max(length_B_vec)
        
        for j in range(0,s):
            ax.quiver(point[0,j], point[1,j], point[2,j], B[0,j], B[1,j], B[2,j], length = fied_length_quiver*length_B_vec[j] /length_B_vec_maxx, normalize = True, color = 'forestgreen')
        
        ax.set_xlabel('X (m)', fontsize = 12, fontweight = 'bold')
        ax.set_ylabel('Y (m)', fontsize = 12, fontweight = 'bold')
        ax.set_zlabel('Z (m)', fontsize = 12, fontweight = 'bold')
        
        data_ax = plt.axis('auto')
        ax.set_xlim(min(data_ax), max(data_ax))
        ax.set_ylim(min(data_ax), max(data_ax))
        ax.set_zlim(min(data_ax), max(data_ax))
        
        plt.show()
        
        print('   Done.')
        print('   -----------------------------------------------------------------------------------   ')

    # contour plot
    if need_contour_plot == 1:
    
        print('   Generating contours...')
        
        from Inputs import cont_x, cont_y, cont_pp, cont_qq
        
        import matplotlib.pyplot as plt
        from matplotlib import cm

        B_MAGNITUDE = sqrt(Bx**2 + By**2 + Bz**2)
        
        plt.figure('Contour Plot |B|')
        cs = plt.contour( cont_x, cont_y, B_MAGNITUDE.reshape([cont_pp,cont_qq]), 100, cmap = cm.jet)
        plt.clabel(cs, inline = 1, fontsize = 10)
        plt.autoscale(enable = True, axis = 'x', tight = True)
        plt.xlabel('x (m)', fontsize = 20, fontweight = 'bold')
        plt.ylabel('y (m)', fontsize = 20, fontweight = 'bold')
        plt.colorbar()
        plt.axis('square')
        plt.grid(True)
        plt.show()
        
        plt.figure('Contour Plot (filled) |B|')
        cs = plt.contourf( cont_x, cont_y, B_MAGNITUDE.reshape([cont_pp,cont_qq]), 100, cmap = cm.jet)
        plt.autoscale(enable = True, axis = 'x', tight = True)
        plt.xlabel('x (m)', fontsize = 20, fontweight = 'bold')
        plt.ylabel('y (m)', fontsize = 20, fontweight = 'bold')
        plt.colorbar()
        plt.axis('square')
        plt.grid(False)
        plt.show()
        
        print('   Done.')
        print('   -----------------------------------------------------------------------------------   ')
