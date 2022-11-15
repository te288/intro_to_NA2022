# ResSimUtils -> Reservoir Simulation Utilitiesの意味
# ここには相対浸透率，透過率行列を返す関数を記述]
import numpy as np
from scipy.sparse import spdiags

def cal_krw_WellTest(Sw, connateSw, resiSo, krw_max, nw):
    # this function computes Relative Permeability of water
    Seff = (Sw - connateSw)/(1- connateSw - resiSo)
    krw  = krw_max*(Seff**nw)
    krw[np.where( Sw <= connateSw)] = 0
    krw[np.where( Sw >= (1-resiSo))] = krw_max
    return krw

def cal_kro_WellTest(Sw, connateSw, resiSo, kro_max, no):
    # this function computes Relative Permeability of oil
    Seff = (1-Sw - connateSw)/(1- connateSw - resiSo)
    kro  = kro_max*((Seff)**no)
    kro[np.where( Sw <= connateSw)]  = kro_max
    kro[np.where( Sw >= (1-resiSo))] = 0
    return kro

def cal_Ta_WellTest(perm_x, perm_rel_a, h, dx, dy, vis_a, FVF_a, Pressure, nx, BC_w=0, BC_e= 0, Pb_w = 0, Pb_e = 0, Perm_w = 0, Perm_e = 0, RelPerm = 0):
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
    bcs_a    = np.zeros(nx) # Boundary condition vector 
    
    # T_west
    for i in range(0, nx):
        if i == 0:
            if BC_w != 1:
                Tw = 0
                T_west[i] = 0
            elif BC_w == 1:
                Tw = h*dy*(2*(perm_x[i] * Perm_w )/(perm_x[i] + Perm_w) )/(vis_a*dx*FVF_a)
                if Pressure[i] < Pb_w:
                    T_wbc = RelPerm * Tw;
                elif Pressure[i] > Pb_w:
                    T_wbc = perm_rel_a[i] * Tw;
                else:
                    T_wbc = 0; #(perm_rel_a[i] + RelPerm)/2 * Tw;
                T_o[i]    = T_wbc
                bcs_a[i]  = T_wbc*Pb_w
                T_west[i] = 0
                
        else:
            Tw = h*dy*(2*(perm_x[i] *perm_x[i-1] )/(perm_x[i] +perm_x[i-1]) )/(vis_a*dx*FVF_a)
            if Pressure[i] < Pressure[i-1]:
                T_west[i] = perm_rel_a[i-1] * Tw;
            elif Pressure[i] > Pressure[i-1]:
                T_west[i] = perm_rel_a[i] * Tw;
            else:
                T_west[i] = 0; #(perm_rel_a[i] + perm_rel_a[i-1])/2 * Tw;
                
    # T_east
    for j in range(0, nx):
        if j == nx-1:
            if BC_e != 1:
                Te = 0
                T_east[j] = 0
            elif BC_e == 1:
                Te = h*dy*(2*(perm_x[j] * Perm_e )/(perm_x[j] + Perm_e) )/(vis_a*dx*FVF_a)
                if Pressure[j] < Pb_e:
                    T_ebc = RelPerm * Te;
                elif Pressure[j] > Pb_e:
                    T_ebc = perm_rel_a[j] * Te;
                else:
                    T_ebc = 0; #(perm_rel_a[j] + RelPerm)/2 * Te;
                T_o[j]    = T_ebc
                bcs_a[j]  = T_ebc*Pb_e
                T_east[j] = 0
        else:
            Te = h*dy*(2*(perm_x[j] *perm_x[j+1] )/(perm_x[j] +perm_x[j+1]) )/(vis_a*dx*FVF_a)
            if Pressure[j] < Pressure[j+1]:
                T_east[j] = perm_rel_a[j+1] * Te;
            elif Pressure[j] > Pressure[j+1]:
                T_east[j] = perm_rel_a[j] * Te;
            else:
                T_east[j] = 0; #(perm_rel_a[j] + perm_rel_a[j+1])/2 * Te;
    
    # To
    T_o = T_o + T_west + T_east
    
    # T_a
    T_a = spdiags([-T_east, T_o, -T_west], [-1, 0, 1], nx, nx).T
    
    return T_a, bcs_a