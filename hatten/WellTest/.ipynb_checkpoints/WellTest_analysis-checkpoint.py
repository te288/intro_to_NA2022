import numpy as np

def cal_Pwf(Q_sc, P_res, Perm, Perm_ra, Tick_res, B_a, mu_a, r_well, skin, dx):
    # this function computes bottom hole pressure
    #- Input -#
    # Q_sc : flow rate of phase a @ surface condition
    # P_res : Pressure of reservoir model grid where well is located
    # Perm : Permeability of reservoir model grid where well is located
    # Perm_ra : Relative Permeability of phase a at reservoir model grid where well is located
    # Tick_res : Tickness of reservoir
    # B_a : Formation Volume Factor of phase a
    # mu_a : viscosity of phase a
    # r_well : radius of well 
    # skin : skin factor of well
    # dx : size of reservoir model block
     
    r_eq = dx*0.2078 # simple peace man correction(1978) is used
    
    PI = 2*np.pi*Perm*Perm_ra*Tick_res/( mu_a*B_a*(np.log(r_eq/r_well)+skin) );
    
    Pwf = P_res - (Q_sc/B_a)/PI;
    
    return Pwf