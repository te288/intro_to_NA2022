# ResSimUtils -> Reservoir Simulation Utilitiesの意味
# ここには相対浸透率，透過率行列を返す関数を記述]
import numpy as np
from scipy.sparse import spdiags

def cal_krw(Sw, connateSw, resiSo, krw_max, nw):
    # this function computes Relative Permeability of water
    
    #- Write Your code here -#
    
    return krw

def cal_kro(Sw, connateSw, resiSo, kro_max, no):
    # this function computes Relative Permeability of oil
    
    #- Write Your code here -#
    
    return kro

def cal_Ta(perm_x, perm_rel_a, h, dx, dy, vis_a, FVF_a, Pressure, nx):
    # ----- Input ----- #
    # perm_x     : Absolute perm. in x-direction
    # perm_rel_a : Relative permeability of phase-a
    # h          : hight of the reservoir
    # dx         : The Size of CV in x-direction
    # dy         : The Size of CV in y-direction
    # vis_a      : viscosity of phase-a
    # FVF_a      : Formation Volume Factor of phase-a
    # nx         : The Number of CVs in x-direction
    # ----- Output ----- #
    # Ta         : Transmissibility Matrix of phase-a
    
    T_west = np.zeros(nx) # Transmissibility between i-1 and i th block
    T_east = np.zeros(nx) # Transmissibility between i+1 and i th block
    T_o    = np.zeros(nx) # diagonal elements of Ta
    
    # T_west
    for i in range(0, nx):
        if i == 0:
            Tw = 0
            T_west[i] = 0
        else:
            Tw = #- Write Your code here -#
            if Pressure[i] < Pressure[i-1]:
                T_west[i] = #- Write Your code here -#
            elif Pressure[i] > Pressure[i-1]:
                T_west[i] = #- Write Your code here -#
            else:
                T_west[i] = #- Write Your code here -#
    # T_east
    for i in range(0, nx):
        if i == nx-1:
            Te = 0
            T_east[i] = 0
        else:
            Te = #- Write Your code here -#
            if Pressure[i] < Pressure[i+1]:
                T_east[i] = #- Write Your code here -#
            elif Pressure[i] > Pressure[i+1]:
                T_east[i] = #- Write Your code here -#
            else:
                T_east[i] = #- Write Your code here -#
    
    # To
    T_o = T_west + T_east
    
    # T_a
    T_a = spdiags([-T_east, T_o, -T_west], [-1, 0, 1], nx, nx).T
    
    return T_a