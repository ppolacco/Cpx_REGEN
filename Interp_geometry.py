# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 10:15:03 2021

@author: nicol
"""


import GeometryFunctions as GeoF
from math import pi
import numpy as np


# The 'interpolation_gen' function returns 3 dictionaries where all the coefficients of the polynomial fit
# are stored. Thanks to these coefficients the variation of all the sections and volume are known throughout a rotation
# of the shaft.
def interpolation_gen(geo):
    # AREA COEFFICIENTS
    # Discretization of the angles going from 0 to 2pi. This discretization is used to approximate section variations
    # (90 intervals)
    theta_area = np.linspace(0,2*pi,91)
    
    # Initialization of the vectors to interpolate
    area_port =  np.zeros(len(theta_area))
    area_suction =  np.zeros(len(theta_area))
    area_discharge =  np.zeros(len(theta_area))
    
    # Port area, suction area and discharge area are computed throughout the whole rotation
    # using the proper function of the geometry module. This is a crucial part in term of
    # computational effort.
    for i in range(len(theta_area)):
        (_,_,area_port[i]) = GeoF.polyShadedPortArea(geo, theta_area[i])
        (_,_,area_suction[i]) = GeoF.suctionFlowArea(geo, theta_area[i])
        (_,_,area_discharge[i]) = GeoF.area_d_dd(geo,theta_area[i])
        
    # conversion from square meters to square millimiters
    area_port = area_port/1e6    
    area_suction = area_suction/1e6
    area_discharge = area_discharge/1e6    
        
    # The variation of port, suction and discharge areas is approximated by means of polynomial approximation.
    coefs_port = np.polyfit(theta_area, area_port, 15)  
    coefs_discharge = np.polyfit(theta_area, area_discharge, 15)    
    coefs_suction = np.polyfit(theta_area, area_suction, 15)  
       
    # VOLUME AND VOLUME DERIVATIVE COEFFICIENTS
    # Discretization of the angles going from 0 to 2pi. This discretization is used to approximate variations
    # of volumes and volumes derivatives (360 intervals)
    theta_vol = np.linspace(0,2*pi,361)
    
    # Initialization of the vectors to interpolate (volumes)
    V_sa_vect = np.zeros(len(theta_vol))
    V_s1_vect = np.zeros(len(theta_vol))
    V_s2_vect = np.zeros(len(theta_vol))
    V_c1_vect = np.zeros((len(theta_vol),geo.N_c_max))
    V_c2_vect = np.zeros((len(theta_vol),geo.N_c_max))
    V_d1_vect = np.zeros(len(theta_vol))
    V_d2_vect = np.zeros(len(theta_vol))
    V_dd_vect = np.zeros(len(theta_vol))
    
    # Initialization of the vectors to interpolate (volume)
    dV_sa_vect = np.zeros(len(theta_vol))
    dV_s1_vect = np.zeros(len(theta_vol))
    dV_s2_vect = np.zeros(len(theta_vol))
    dV_c1_vect = np.zeros((len(theta_vol),geo.N_c_max))
    dV_c2_vect = np.zeros((len(theta_vol),geo.N_c_max))
    dV_d1_vect = np.zeros(len(theta_vol))
    dV_d2_vect = np.zeros(len(theta_vol))
    dV_dd_vect = np.zeros(len(theta_vol))
    

    # Volumes,volume derivatives and other geometrical quantities are computed throughout the whole rotation
    # using the proper function of the geometry module. This is a crucial part in term of
    # computational effort.
    for i in range(len(theta_vol)): 
        dict_sa = GeoF.volume_force_sa(geo,theta_vol[i],False)
        dict_s1 = GeoF.volume_force_s1(geo,theta_vol[i],False)
        dict_s2 = GeoF.volume_force_s2(geo,theta_vol[i],False)
        dict_c1 = GeoF.volume_force_c1(geo,theta_vol[i],False)
        dict_c2 = GeoF.volume_force_c2(geo,theta_vol[i],False)
        dict_d1 = GeoF.volume_force_d1(geo,theta_vol[i],False)
        dict_d2 = GeoF.volume_force_d2(geo,theta_vol[i],False)
        dict_dd = GeoF.volume_force_dd(geo,theta_vol[i],False)
        # VOLUMES
        V_sa_vect[i] = dict_sa["V"]/1e9
        V_s1_vect[i] = dict_s1["V"]/1e9
        V_s2_vect[i] = dict_s2["V"]/1e9
        for alpha in range(1,geo.N_c_max + 1,1):
            V_c1_vect[i,alpha-1] = dict_c1[alpha]["V"]/1e9
            V_c2_vect[i,alpha-1] = dict_c2[alpha]["V"]/1e9
        V_d1_vect[i] = dict_d1["V"]/1e9
        V_d2_vect[i] = dict_d2["V"]/1e9
        V_dd_vect[i] = dict_dd["V"]/1e9
        # VOLUME DERIVATIVES
        dV_sa_vect[i] = dict_sa["dVdTheta"]/1e9
        dV_s1_vect[i] = dict_s1["dVdTheta"]/1e9
        dV_s2_vect[i] = dict_s2["dVdTheta"]/1e9
        for alpha in range(1,geo.N_c_max + 1,1):
            dV_c1_vect[i,alpha-1] = dict_c1[alpha]["dVdTheta"]/1e9
            dV_c2_vect[i,alpha-1] = dict_c2[alpha]["dVdTheta"]/1e9
        dV_d1_vect[i] = dict_d1["dVdTheta"]/1e9
        dV_d2_vect[i] = dict_d2["dVdTheta"]/1e9
        dV_dd_vect[i] = dict_dd["dVdTheta"]/1e9

    # The variation of volumes and volume derivative is approximated by means of polynomials.
    # 'n' represents the grade of the polynomial used to approximate the variation
    n = 5
    coefs_Vsa = np.array([np.polyfit(theta_vol, V_sa_vect, n)])         
    coefs_Vs1 = np.array([np.polyfit(theta_vol, V_s1_vect, n)  ])  
    coefs_Vs2 = np.array([np.polyfit(theta_vol, V_s2_vect, n) ])  
    coefs_dVsa = np.array([np.polyfit(theta_vol, dV_sa_vect, n) ])  
    coefs_dVs1 = np.array([np.polyfit(theta_vol, dV_s1_vect, n)  ])  
    coefs_dVs2 = np.array([np.polyfit(theta_vol, dV_s2_vect, n)])  
    
    # The variation of the compression chamber is linear so it can be approximated using a polynomial
    # of grade 1 (a line, n = 2).
    coefs_Vc1 = np.zeros(( geo.N_c_max,2))
    coefs_Vc2 = np.zeros((geo.N_c_max,2))
    coefs_dVc1 = np.zeros(( geo.N_c_max,2))
    coefs_dVc2 = np.zeros((geo.N_c_max,2))
    
    # The variation of the compression chambers' volume and volume derivative is approximated with polynomials.
    # In general there are N expansion chambers in a scroll machine. The number of compression chamber is 
    # stored as an attribute of the 'geo' object. 
    for alpha in range(1,geo.N_c_max +1 ,1):
        coefs_Vc1[alpha-1,:] = np.polyfit(theta_vol[np.isfinite(V_c1_vect[:,alpha-1])], V_c1_vect[np.isfinite(V_c1_vect[:,alpha-1]),alpha-1], 1) 
        coefs_Vc2[alpha-1,:] = np.polyfit(theta_vol[np.isfinite(V_c2_vect[:,alpha-1])], V_c2_vect[np.isfinite(V_c2_vect[:,alpha-1]),alpha-1], 1) 
        coefs_dVc1[alpha-1,:] = np.polyfit(theta_vol[np.isfinite(dV_c1_vect[:,alpha-1])], dV_c1_vect[np.isfinite(dV_c1_vect[:,alpha-1]),alpha-1], 1) 
        coefs_dVc2[alpha-1,:] = np.polyfit(theta_vol[np.isfinite(dV_c2_vect[:,alpha-1])], dV_c2_vect[np.isfinite(dV_c2_vect[:,alpha-1]),alpha-1], 1) 
    
    # The variation of the 'ddd' chamber is trickier since it is not always defined during the rotation.
    theta_vol_d = np.linspace(0,4*pi,721)
    theta_d = geo.theta_d
    
    # Initialization of the vectors to interpolate
    V_ddd_vect_bis = np.zeros(len(theta_vol_d))
    V_d1_vect_bis = np.zeros(len(theta_vol_d))
    
    # The existence of the 'ddd' chamber is checked and a proper value for its volume is stored
    # inside the correspondent vectors.
    for i in range(len(theta_vol_d)): 
        if theta_vol_d[i] <= theta_d:
            V_ddd_vect_bis[i] = np.NaN
            V_d1_vect_bis[i] = np.NaN
        elif  theta_d < theta_vol_d[i] and theta_vol_d[i] <= 2*pi:
            V_ddd_vect_bis[i] = V_d1_vect[i] + V_d2_vect[i] + V_dd_vect[i]
            V_d1_vect_bis[i] = V_d1_vect[i]
        elif theta_vol_d[i] > 2*pi and theta_vol_d[i] < theta_d + 2*pi - 2*pi/360:
            V_ddd_vect_bis[i] = V_dd_vect[i - 360] + V_d2_vect[i - 360] + V_d1_vect[i - 360]
            V_d1_vect_bis[i] = V_d1_vect[i - 360]
        else:
            V_ddd_vect_bis[i] = V_dd_vect[i -360] 
            V_d1_vect_bis[i] = np.NaN

    # Approximation of the variation of the 'ddd' chamber using polynomial approximation.
    # (the approximation for the 'ddd' chamber is made using a polynomial of grade 2*n).
    coefs_Vddd = np.array([np.polyfit(theta_vol_d[np.isfinite(V_ddd_vect_bis)], V_ddd_vect_bis[np.isfinite(V_ddd_vect_bis)], n*2) ])  
    coefs_Vd1 = np.array([np.polyfit(theta_vol_d[np.isfinite(V_d1_vect_bis)], V_d1_vect_bis[np.isfinite(V_d1_vect_bis)], n)    ])  
    coefs_Vdd = coefs_Vddd

    # Procedure to compute coefficients of polynomial for volume derivative of chamber 'ddd'.
    coefs_dVddd = coefs_Vddd[:,0:len(coefs_Vddd[0,:])-1].copy()
    L_ddd = len((coefs_dVddd[0,:]))
    for i in range(L_ddd):
        coefs_dVddd[:,i] = coefs_dVddd[:,i]*(L_ddd-i)
    coefs_dVdd = coefs_dVddd
    
    # Procedure to compute coefficients of polynomial for volume derivative of chamber 'd1'.
    coefs_dVd1 = coefs_Vd1[:,0:len(coefs_Vd1[0,:])-1].copy()
    L_d = len((coefs_dVd1[0,:]))
    for i in range(L_d):
        coefs_dVd1[:,i] = coefs_dVd1[:,i]*(L_d-i)

    # Chamber 'd2' is symmetrical to chamber 'Vd1'.
    coefs_Vd2 = coefs_Vd1
    coefs_dVd2 = coefs_dVd1

    # The dictionaries that are the output of this function are set.
    dict_area = {"s" : coefs_suction, "d" : coefs_discharge, "port" : coefs_port}
    dict_V = {"sa" : coefs_Vsa, "s1" : coefs_Vs1, "s2" : coefs_Vs2, "c1" : coefs_Vc1, "c2" : coefs_Vc2, "d1" : coefs_Vd1, "d2" : coefs_Vd2, "dd" : coefs_Vdd, "ddd" : coefs_Vddd}
    dict_dV = {"sa" : coefs_dVsa, "s1" : coefs_dVs1, "s2" : coefs_dVs2, "c1" : coefs_dVc1, "c2" : coefs_dVc2, "d1" : coefs_dVd1, "d2" : coefs_dVd2, "dd" : coefs_dVdd, "ddd" : coefs_dVddd}

    return dict_V, dict_dV, dict_area


