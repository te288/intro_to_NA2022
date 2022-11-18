####### ----- Reservoir Properties and Simulation control ---- ######
# Simulation control Parameters
dt       = 1/24/6 # initial time step[day]. time step might be revised to meet CFL condition
t_max    = 10 #[day] simulation will continues until t becomes t_max
n_max    = 1e8 # maximum interation [-]
n_out    = 6 # 1000以上がいいと思う。0 means no output while simulation
Pmin     = 995 # Pressure lower limit for plot
Pmax     = 1000 # Pressure lower limit for plot
saveResP_contour = 0    # save png 
showResP_contour = 0    # show contour

saveResP_plot    = 0    # save png 
showResP_plot    = 0    # show plot

saveSw_contour   = 0    # save png
showSw_contour   = 0    # show contour

saveSw_plot      = 0    # save png 
showSw_plot      = 0    # show plot

fig1             = plt.figure(); # figure for pressure contour
fig2             = plt.figure(); # figure for pressure plot
fig3             = plt.figure(); # figure for saturation contour
fig4             = plt.figure(); # figure for saturation plot

# Reservoir Properties
Lx = 10000 # Reservoir Length in x direction [ft]
Ly = 10 # Reservoir Length in y direction [ft] # to avoid y-
h  = 20   # Reservoir tickness in z direction [ft]
nx = 31   # number of grid in x direction [-]
ny = 1   # number of grid in y direction [-]
dx = Lx / nx # grid size in x direction
dy = Ly / ny # grid size in y direction
perm_x  = 100*np.ones(nx) # Absolute permeability in x direction [mD]
# perm_y  = 100*np.ones(ny, nx) # Absolute permeability in y direction [mD]
phi_res = 0.3*np.ones(nx) # Porosity of reservoir [-]
P_init  = 3900*np.ones(nx) # Initial pressure in reservoir [PSI]

# define coordinate system
x    = np.zeros(nx)
x[0] = 0.5*dx
#y = np.np.zeros(ny)
#y[0] = 0.5*dy

# x coordinate
for i in range(1,len(x)):
    x[i] = x[i-1] + dx

###### ----- Fluid Properties ---- ######
# Injected water properties
vis_w        = 1.0 # viscosity of water [cp]
Bw           = 1.0 # formation volume factor of water
cw           = 1e-5 # compressibility of water [psi^-1]
Sw_i         = 0.2 # irreducible water saturation
Sw_init      = 0.2*np.ones(nx) # initial water saturation distribution
#Sw_init      = np.linspace(0, 1, nx)

perm_rw_max  = 1 # maximuam relative perm. of water [-]
nw           = 3.0 # exponent of water relative permeability curve [>1]

# Oil water Properties
vis_o        = 1.0 # viscosity of oil [cp]
Bo           = 1.1 # formation volume factor of oil
co           = 1e-5 # compressibility of oil [psi^-1]
So_r         = 0.2 # residual oil saturation
So_init      = 1 - Sw_init # initial oil saturation progile
perm_ro_max  = 1 # maximuam relative perm. of oil [-]
no           = 3.0 # exponent of oil relative permeability curve [>1]

###### ----- Boundary Condition (Aquifer) ---- ######
aquiferWest      = 1; # Aquifer exists at x = 0
aquiferEast      = 1; # Aquifer exists at x = Lx
perm_Aq_w        = 100;
perm_Aq_e        = 100;
perm_rw_Aq       = np.max(cal_krw_WellTest(np.ones(2), Sw_i, So_r, perm_rw_max, nw));
perm_ro_Aq       = np.max(cal_kro_WellTest(np.ones(2), Sw_i, So_r, perm_rw_max, nw));
Pb_w             = 3900;
Pb_e             = 3900;

###### ----- Well, Production & Injection Information ---- ######
# Injection Well Properties
# injection wells are controlled by flow rate[ft^3 /day]. Injecters are
# composed of water flow rate. Producers are only controlled by total flow
# rate.

# Injection wells (source term)
#Qin = np.zeros(ny, nx) # Injection wells
Qin = np.zeros(nx) # Injection wells
#loc_Qin  = [16]
#Qin[loc_Qin[0]]   = 2*426.5

loc_Qin  = [0]
Qin[loc_Qin[0]]   = 0

# Production Well Properties
# production wells are controlled by flow rate[ft^s /day] or constant
# Bottom Hole Pressure, Pwf [PSI]. If Pwf is larger than Pi (mean pressure
# in CV), the situation is regarded as no production.

# Production wells (sink term) controlled by Flow Rate --------------------
#Qout = np.zeros(ny, nx) # Production wells
Qout = np.zeros(nx)
#Qout_loc  = [0, 29]
#Qout_name = ['P-1', 'P-2']
#Qout[Qout_loc[0]]   = -426.5
#Qout[Qout_loc[1]]   = -426.5

Qout_loc  = [int((nx-1)/2)]
Qout_name = ['P-1']
Qout[Qout_loc[0]]   = -701 # [res cu ft / day], 125[rbbl/day]


t_open  = 0 # after t_shut [day], all of production wells are openedf 0
t_shut  = 0 # after t_shut [day], all of production wells are shut down, if 0 , wells are not shut down

# Array to hold Production History
ProdOil   = np.zeros([len(Qout_loc)]) # numpy array to hold Production History of Oil
ProdWater = np.zeros([len(Qout_loc)]) # numpy array to hold Production History of Water
ProdResP  = np.zeros([len(Qout_loc)]) # numpy array to hold Pressure Change at reservoir block where Production wells are located
t_hist    = []

###### ----- Show Relative Permeability Curve ---- ######
Sw_sample = np.linspace(0, 1, 1000);
krw_sample = cal_krw_WellTest(Sw_sample, Sw_i, So_r, perm_rw_max, nw)
kro_sample = cal_kro_WellTest(Sw_sample, Sw_i, So_r, perm_ro_max, no)

fig_RelPerm = plt.figure();
ax_rp  = fig_RelPerm.add_subplot(1,1,1)
ax_rp.plot(Sw_sample, krw_sample, color='#0000cd', label= 'Rel. Perm. Water')
ax_rp.plot(Sw_sample, kro_sample, color='#006400', label= 'Rel. Perm. Oil')
ax_rp.set_xlim(0, 1)
ax_rp.set_ylim(-0.05, 1.05)
ax_rp.grid()
plt.show()