# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 15:48:43 2021

@author: nicol
"""

import solverFunctions as solverF
import GeometryFunctions as geoF
from math import pi
from Interp_geometry import interpolation_gen
from scipy.io import loadmat
import propertyFunctions as propf 

import time

tic = time.perf_counter()

# Import a property matrix of the fluid R1233zdE created in Matlab.
# It allows to avoid the computation of a property through CoolProp of Refprop
property_matrix = loadmat('R1233zd.mat')
#%% Geometric parameters definition

# Sanden Machine

# phi_i0 = 84*pi/180
# phi_is = 318*pi/180
# phi_ie_fix = 1052*pi/180
# phi_ie_orb =  phi_ie_fix
# #phi_ie = 5*pi
# phi_o0 = 0 
# phi_os = 145*pi/180
# phi_oe_fix = phi_ie_fix
# phi_oe_orb = phi_ie_fix
# #phi_oe = 5*pi
# r_b = 3.15 #mm
# h_s = 30.6 #mm
# D_wall = 135 #mm

# # Discharge port
# x_dis0 = -2.683
# y_dis0 = -4.58
# r_dis = 5.5

# # Leakage 
# delta_r = 20e-6 #m
# delta_f = 20e-6 #m 

# tau_loss = 0.3 #Nm

# geo = GeoParam(phi_i0, phi_is, phi_ie_fix, phi_ie_orb, phi_o0, phi_os, phi_oe_fix, phi_oe_orb, r_b, h_s, D_wall, x_dis0, y_dis0, r_dis, delta_r, delta_f)
# geo.dischargeGeo( type = 'aala', r_a1 = 6.5, r_a2 = 1, angle0 = 45*pi/180, r_a0 = 8.5, dev = 24.166*pi/180, dist = 11.5)
# geo.important_values()


# Exoes Machine 

phi_i0 = 81*pi/180
phi_is = 387*pi/180 #330
phi_ie_fix = 1583*pi/180
phi_ie_orb =  phi_ie_fix
phi_o0 = 0
phi_os = phi_is - pi
phi_oe_fix = phi_ie_fix
phi_oe_orb = phi_ie_fix
r_b = 3.2 #mm
h_s = 41/2 #mm # 20.5
D_wall = 95*2 #mm

# Discharge port
x_dis0 = -1.17
y_dis0 = 3.52
r_dis = 4.3 #3.9
L_dis =  9
x_dis1 = -3.5 # -4.05

A_dis = True # if true, a function is used to compute the area as a function of the built-in volume ratio 

# Leakage 
delta_r = 20e-6 #m
delta_f = 20e-6 #m 

tau_loss = 0.3 #Nm

# GeoParam is a class impored from the 'GeometryFunctions' module. It contains all the geometry parameters defined
# for the current machine geometry.
geo = geoF.GeoParam(phi_i0, phi_is, phi_ie_fix, phi_ie_orb, phi_o0, phi_os, phi_oe_fix, phi_oe_orb, r_b, h_s, D_wall, x_dis0, y_dis0, r_dis, delta_r, delta_f, L_dis, x_dis1)
# 'dischargeGeo' is a method of the GeoParam class. It checks whether the discarge region geometry is compatible with a
# Perfect Meshing Profile (PMP). It computes the maximum r_a2 and correct the value of r_a1 to make sure that PMP is attained
# if possibile. The 'type' input refers to the kind of discharge profile ('2a' = two-arc, 'ala' = arc-line-arc).
geo.dischargeGeo( type = '2a', r_a1 = 12.405, r_a2 = 2)
# The 'important_values' method computes discharge angle, uncontact angle, maximum number of compression chamber,
# displacement volume and volume ratio for the current geometry
geo.important_values()

# In this section the geometry is loaded into a dictionary ('geo_inputs') that is called throughout the code.
# The 'interpolation_gen' function performs the interpolation of the chambers using polygons.
print('Loading geometrical model...')
(dict_V, dict_dV, dict_area) = interpolation_gen(geo)
geo_inputs = {'geo' : geo, 'phi_i0' : phi_i0, 'phi_o0' : phi_o0, 'phi_ie' : phi_ie_fix, 'phi_oe' : phi_oe_fix,
             'r_b' : r_b, 'h_s' : h_s, 'D_wall' : D_wall, 'x_dis0' : x_dis0, 'y_dis0' : y_dis0, 'r_dis' : r_dis,
             'dict_V': dict_V, 'dict_dV' : dict_dV, 'dict_area' : dict_area, 'tau_loss' : tau_loss}
print('Loading geometrical model: done')

#%% DEFINITION OF THE MODEL

# Initialization of the model
Heat_transfer = False # not yet taken into account
Leakage = True
Reed_valve = True
Bearing = False
N_rot = 1450 # Compressor speed - RPM (3100 - 1450)
omega = N_rot/60*2*pi

# Thermodynamic properties at inlet of the compressor
T_in = 90 + 273 # Inlet temperature - K (30 + 273 - 90 + 273)
Q_in =  0.466  # Inlet quality (0.39 - 0.466)
P_in = 160000 # Inlet pressure - Pa - needed if the input is not two-phase, useless otherwise
x_l_in = 0 # oil rate 

# Pressure ratio between inlet and outlet
ratio = 2.64 # (5.5625 - 2.64)

# 'property_mixture' is a function that computes a required property, given 2 independent state variables.
# This function takes as an input a matrix containing the properties of R1233zd ('property_matrix')
P_in_sat = propf.propfunction_mixture(property_matrix, x_l_in, 'P', T = T_in, Q = 0.5) 
# for i in range(len(Q_in)):
if (P_in_sat >= P_in) and Q_in >= 1:
    raise Exception('Not suitable inlet pressure for a quality of one')
    
# Definition of outlet pressure and inlet density. If the fluid is in the two-phase region, pressure is
# uniquely defined by temperature and viceversa. 
if Q_in <= 0 or Q_in >= 1:
    rho_in = propf.propfunction_mixture(property_matrix, x_l_in, 'rho', P = P_in, T = T_in) 
    P_out = P_in * ratio
else:
    P_out =  propf.propfunction_mixture(property_matrix, x_l_in, 'P', Q = 0.5, T = T_in)*ratio
    P_in =  propf.propfunction_mixture(property_matrix, x_l_in, 'P', Q = 0.5, T = T_in)
    rho_in = propf.propfunction_mixture(property_matrix, x_l_in, 'rho', Q = Q_in, T = T_in) 

# A dictionary with all the inputs is generated
inputs = {'HT' : Heat_transfer, 'LK' : Leakage, 'Q_in' : Q_in, 'T_in' : T_in, 'P_in' : P_in, 'P_out' : P_out, 'ratio' : ratio, 
          'N' : N_rot , 'x_l_in' : x_l_in}

# A dictionary with all the outputs is initialized
results = {'geo_inputs' : geo_inputs, 'inputs' : inputs, 'eta_is' : [], 'eta_v' : [], 'm_dot' : [], 'W_dot' : [], 'CV' : [], 'radial_force' : [], 'tangential_force' : [], 'axial_force' : [], 'tilting_moment' : []}

#%% Resolution      

# Definition of tollerance and crank angle step of the numerical scheme
tol_glob =  0.001 # IT MUST BE 0.001
tol_pressure = 0.0025 # IT MUST BE 0.0025
theta_step = 0.001*2

# The solver is called in the function 'full_resolution'. The output of this function are mass flow rate, isentropic
# efficiency, volumetric efficiency, absorbed power and properties within each control volumes ('CV').
# for i in range(len(ratio)):
(m_dot_final, eta_is_final, eta_v_final, W_dot_final, CV, axial_force, tangential_force, radial_force, tilting_moment) = solverF.full_resolution(geo, property_matrix, 
    dict_V, dict_dV, dict_area, T_in, rho_in, x_l_in, P_out, tol_glob, tol_pressure, theta_step, omega, tau_loss, Leakage, Heat_transfer, Reed_valve, Bearing)

# The result of the simulation are stored inside the 'results' dictionaries
results['eta_is'].append(eta_is_final)
results['eta_v'].append(eta_v_final)
results['m_dot'].append(m_dot_final)
results['W_dot'].append(W_dot_final)
results['CV'] = CV
results['radial_force'] = radial_force
results['tangential_force'] = tangential_force
results['axial_force'] = axial_force
results['tilting_moment']= tilting_moment

#%% Save results
# Name of the file that will be saved inside the specified directory. It contains the results of the simulation performed.
filename = 'test'
filename = solverF.file_name(filename)

# The results are saved
solverF.save_results(filename, results)

#%% Time taken
toc = time.perf_counter()
print('Simulation time:', str(toc - tic), 's')

