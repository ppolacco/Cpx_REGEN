# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 11:23:04 2021

@author: nicol
"""

import numpy as np
import propertyFunctions as propF

# This function resolves the energy balance equation given the rotational speed and the dictionaries in which
# mass flow in and out of a given control volume are specified.
def dmCV_dTheta(omega, m_dot_dict):
    """
    Parameters
    ----------
    omega : float.
    m_dot : dict.

    Returns
    -------
    dmCVdTheta.
    """
    # Initialization of the net mass flow rate
    m_dot_sum =  0 
    
    # Computation of the net mass flow rate
    for i, j in m_dot_dict.items():
        m_dot_sum += m_dot_dict[i] 
    
    # The rate of change of mass inside a control volume with respect to the angle theta is obtained as the ratio
    # between the net mass flow rate and the rotation speed of the machine.
    dmCVdTheta = 1/omega*m_dot_sum
    
    return dmCVdTheta
    
# This function resolves the oil mass balance to compute the rate change of oil mass fraction withing a contro volume,
# given rotational speed, mass inside the control volume, the dictionaries in which mass flow in and out of a given
# control volume are specified and oil mass fraction inside a control volume.
def dxoil_dTheta(omega, m_cv, dmCVdTheta, m_dot_dict, x_l_dict, x_l):
    
    # Initialization of the oil net mass flow rate
    m_dot_x_l_sum = 0 

    # Computation of the oil net mass flow rate
    for i, j in m_dot_dict.items():
        m_dot_x_l_sum += m_dot_dict[i]* x_l_dict[i]     
   
    # The rate of change of oil mass fraction is obtained as the difference between rate of change of oil mass due to
    # oil mass flow rate and overall mass flow rate, divided by the mass inside the control volume.
    dxoildTheta = 1/m_cv*(1/omega*m_dot_x_l_sum - x_l*dmCVdTheta)

    return dxoildTheta

# Temperature derivatives with respect to the crank angle is computed by solving the energy balance inside
# the control volume. 
def dT_dTheta(property_matrix, omega, m_cv, T, rho, dVdTheta, dmCVdTheta, dxoildTheta, m_dot_dict, h_dict, x_l, Q_dot = 0):
  
    # Initialization of the rate of change of temperature
    dTdTheta = 0
  
    # The energy balance equation can be decomposed into 5 terms accounting for 5 different mechanisms through
    # which energy is transferred from and to the control volume
  
    # Term 1 - Work exerted by the change of volume
    dPdT = propF.propfunction_mixture(property_matrix, x_l, 'dPdT', T = T, rho = rho)
    term1 = -T*dPdT*(dVdTheta - 1/rho*dmCVdTheta)
    
    # Term 2 - Energy contained in the oil mass fraction
    u_l = propF.propfunction_mixture(property_matrix, 1, 'u', T = T, rho = rho)
    u_g = propF.propfunction_mixture(property_matrix, 0, 'u', T = T, rho = rho)
    term2 = -m_cv*(u_l - u_g)*dxoildTheta
    
    # Term 3 - Energy trasported by the fluid coming into the control volume and going out of it
    h =  propF.propfunction_mixture(property_matrix, x_l, 'h', T = T, rho = rho)   
    term3 = -h*dmCVdTheta
    
    # Term 4 - Heat exchange in the control volume
    term4 = Q_dot/omega
    
    # Term 5 - Rate of change of enthalpy
    m_dot_h_sum = 0 
    for i, j in m_dot_dict.items():
        m_dot_h_sum += m_dot_dict[i]* h_dict[i]     
    term5 = 1/omega* m_dot_h_sum
    
    # Rate of change of internal energy
    dudT = propF.propfunction_mixture(property_matrix, x_l, 'dudT', T = T, rho = rho)
    
    # All the terms are taken into account in a single equation. Temperature derivative is the only unknwon term
    # in the resulting equation
    dTdTheta = (term1 + term2 + term3 + term4 + term5)/dudT/m_cv
    
    return  dTdTheta

# ***???***
def supply_HT(HT_bool, T_flow_in, m_dot = 0, T_wall = 0):
    
    if HT_bool != True:
        T_flow_out = T_flow_in
        
    return T_flow_out

# The 'discharge_is_guess' function is meant to compute a first-attempt value of temperature and density in the discharge
# chambers as they are initialized. Temperature and density at machine inlet and pressure at outlet are known.
# The isentropic efficiency used is a first-attempt value that is specified as input (isentropic efficiency for this
# kind of machine is the rage 0.6-0.7)
def discharge_is_guess(property_matrix, rho_in, T_in, P_out, x_l, eta_is_guess):
    
    # Enthalpy and entropy at inlet are computed based on temperature and density at inlet
    h_in = propF.propfunction_mixture(property_matrix, x_l, 'h', T = T_in, rho = rho_in)
    s_in = propF.propfunction_mixture(property_matrix, x_l, 's', T = T_in, rho = rho_in)
    
    # Enthalpy at outlet is computed based on the isentropic assumption
    h_is_guess = propF.propfunction_mixture(property_matrix, x_l, 'h', s = s_in, P = P_out)   
    
    # Enthalpy at outlet is computed using the first-attempt value of the isentropic efficiency provided as input
    h_out_guess = (h_is_guess - h_in)/eta_is_guess + h_in
    
    # Temperature and density at outlet are computed based on the outlet enthalpy guessed
    T_out_guess = propF.propfunction_mixture(property_matrix, x_l, 'T', h = h_out_guess, P = P_out)
    rho_out_guess = propF.propfunction_mixture(property_matrix, x_l, 'rho', h = h_out_guess, P = P_out)
    
    # Solver to implement when adding the oil
    return T_out_guess, rho_out_guess


# The 'discharge_is_guess' function is meant to compute a first-attempt value of temperature and density in the compression
# chambers as they are initialized. Temperature and density at machine inlet and pressure at outlet are known.
# The isentropic efficiency used is a first-attempt value that is specified as input (isentropic efficiency for this
# kind of machine is the rage 0.6-0.7) 
def compression_is_guess(property_matrix, dict_V_c, alpha, rho_in, T_in,  x_l, eta_is_guess):
    
    # Entropy at inlet are computed based on temperature and density at inlet
    s_in = propF.propfunction_mixture(property_matrix, x_l, 's', T = T_in, rho = rho_in)
    
    # Volume of the compression chamber as the process begins (alpha = 0 refers to the outermost compression chamber)
    V_init = np.polyval(dict_V_c[0], 0) 
    # Volume of the compression chamber at the time the current function is called (alpha > 1)
    V_final = np.polyval(dict_V_c[alpha - 1], 0)
    
    # To compute the first-attempt density of the compression chamber, the assumption made is that mass inside the 
    # compression chamber is constant throughtout the whole process
    rho_comp_guess = rho_in*(V_init/V_final)
    # First-attempt temperature is computes based on entropy and density.
    T_comp_guess = propF.propfunction_mixture(property_matrix, x_l, 'T', s = s_in, rho = rho_comp_guess)
    
    # Solver to implement when adding the oil
    return T_comp_guess, rho_comp_guess

# The 'isentropic_nozzle' function computes mass flow rate, enthalpy and Reynolds number of a fluid flowing through 
# an orifice that could be asumed to behave as an isentropic nozzle. Flows though the suction port and the discarge
# port are modeled as isentropic nozzles in the current machine model.
# It is not specified at this stage what is the upstream state and downstream state.
# This function is called within a class method. Properties of the chamber that is the object of the method are
# specified using number 1, properties of the external chamber using number 2.
# i.e.: if we are computing the mass flow rate of the leakage that takes place from compression chamber 'c1' to
# suction chamber 's1', T1 represents temperature of chamber 'c1', T2 represents temperature of chamber 's1'.
def isentropic_nozzle(geo, property_matrix, area, rho1, T1, x_l1, rho2 = None, T2 = None,  x_l2 = None ,h2 = None, P2 = None, leakage_type = None):
    
    # 'flow_function' is defined within the 'isentropic_nozzle' environment and as such they share the same variables.
    # This function is used to compute the mass flow rate of two-phase flow flowing through an isentropic nozzle.
    def flow_function(area, s_up, P_up, P_down , v_up, v_down):   
        #     Parameters for two phase flow
        # area correction factor
        Xd = 1;
        # Morris correction coefficient
        Cd = 0.77;
        # homogenous flow coefficient
        psi = 1;
        # area ratio between throat and suction
        sigma = 0;
        # Number of points for integral calculation
        N_p = 10;
        
        # Discretization of the pressure variation through the nozzle
        D_p = (P_up - P_down)/(N_p-1)
        P = np.linspace(P_up, P_down, N_p, endpoint =True)   
        
        # Initialization of the specific volume
        v = np.ones(len(P))*v_up
 
        # Accordingly with the eqution presented in Ian Bell's thesis, to compute the mass flow rate it
        # is necessary to perform the integral of the specific volume over the whle range of pressure.
        integral = 0
        # The initial value of the specific volume is computed in correspondence of the upstream properties
        v[0] = 1/propF.propfunction_mixture(property_matrix, x_l_up, 'rho', s = s_up , P = P[0])
        for i in range(len(P)-1):
            # The i-th value of the specifiv volume is computed based on the dicretized pressure and assuming
            # an isentropic process
            v[i+1] = 1/propF.propfunction_mixture(property_matrix, x_l_up, 'rho', s = s_up , P = P[i+1]);
            # The value of the integral is updated after each dp considered
            integral = integral + (v[i+1] + v[i])/2*D_p
        
        # The final equation presented in Ian Bell's thesis is solved
        m_dot = Cd*Xd*area*(2*integral/(v_down**2 - sigma*v_up**2))**0.5
    
        return m_dot
   
    # area = 6e-11
    # rho1 = 23
    # rho2 = 22
    # T1 = 302
    # T2 = 303
    # x_l1 = 0
    # x_l2 = 0
    
    # If the flow ares is closed the mass flow rate is zero
    if area == 0: 
        m_dot = 0   
        h = 0
    # If the area is non-zero the mass flow rate computation begins.
    else :
        # Pressure 1 is computed based on temperature and density.
        # Pressure 2 could be computed using either enthalpy and pressure or temperature and density.
        if P2 == None and h2 == None :
            P1 = propF.propfunction_mixture(property_matrix, x_l1, 'P', T = T1 , rho = rho1)
            P2 = propF.propfunction_mixture(property_matrix, x_l2, 'P', T = T2 , rho = rho2)
        elif rho2 == None and T2 == None: 
            P1 = propF.propfunction_mixture(property_matrix, x_l1, 'P', T = T1 , rho = rho1)
            rho2 = propF.propfunction_mixture(property_matrix, x_l2, 'rho', h = h2 , P = P2)
            T2 = propF.propfunction_mixture(property_matrix, x_l2, 'T', h = h2 , P = P2)
        # If input are not coherent the function raises an error.
        else:
            raise Exception('Invalid inputs')
            
        # This is the minimum pressure ratio that is required for existence of flow
        DP_floor = 1; #Pa
    
        # If delta_P is higher than the minimum required the computation takes place.
        if abs(P1 - P2) > DP_floor:
            # Upstream and downstream properties are defined based on the value of pressure in the two chambers
            if P1 > P2:
                x_l_up = x_l1
                P_up = P1
                T_up = T1
                rho_up = rho1 
                h_up = propF.propfunction_mixture(property_matrix, x_l1, 'h', T = T_up , rho = rho_up) 
                visc_up = propF.propfunction_mixture(property_matrix, x_l1, 'mvisc', h = h_up, P = P_up)
                v_up = 1/rho_up
                s_up = propF.propfunction_mixture(property_matrix, x_l1, 's', T = T_up , rho = rho_up)
                x_down = x_l2
                P_down = P2
                T_down = T2
                rho_down = propF.propfunction_mixture(property_matrix, x_l1, 'rho', s = s_up , P = P_down)
                v_down = 1/rho_down
            else:
                x_l_up = x_l2
                P_up = P2
                T_up = T2
                rho_up = rho2
                h_up = propF.propfunction_mixture(property_matrix, x_l1, 'h', T = T_up , rho = rho_up)
                visc_up = propF.propfunction_mixture(property_matrix, x_l1, 'mvisc', h = h_up, P = P_up)
                v_up = 1/rho_up
                s_up = propF.propfunction_mixture(property_matrix, x_l1, 's', T = T_up , rho = rho_up)
                x_down = x_l1
                P_down = P1
                T_down = T1
                rho_down = rho_down = propF.propfunction_mixture(property_matrix, x_l1, 'rho', s = s_up , P = P_down)
                v_down = 1/rho_down
           
            # Computation of the critical pressure ratio. If pressure ratio exceeds this value sonic condition
            # in the nozzle are attained
            k_star = propF.propfunction_mixture(property_matrix, x_l2, 'k_star', T = T_up , rho = rho_up)
            ratio_crit = (1 + (k_star - 1)/2)**(k_star/(1 - k_star))
            # Pressure ratio between downstream chamber and upstream chamber. This ratio is always < 1. 
            ratio = P_down/P_up
            
            # If pressure ratio is lower than the critical pressure ratio (upstream pressure >> downstream pressure)
            # the flow is choked and flow properties do not depend from downstream state anymore.
            if ratio <= ratio_crit:
               P_down = P_up*ratio_crit
               v_down = 1/propF.propfunction_mixture(property_matrix, x_l2, 'rho', s = s_up , P = P_down)   
            # Once all the properties are defined, the mass flow rate is computed by means of the 'flow_function'
            m_dot = flow_function(area, s_up, P_up, P_down , v_up, v_down)       
            # The mass flux is considered to be positive if it enters the chamber on which the method is applied
            if P1 > P2:
                m_dot = -m_dot
                h = h_up            
            if P1 < P2:
                m_dot = m_dot
                h = h_up
        
        # If delta_P is lower than the minimum required the flow through the nozzle is zero.
        else:
            m_dot = 0   
            h = 0
        
    # Computation of Reynolds number of the flow.
    # Re = 0 if mass flow rate is zero
    if leakage_type == None or area == 0 or abs(P1 - P2) <= DP_floor:
        Re = 0 
    # Radial leakage Reynolds number
    elif leakage_type == 'r':
        Re = abs(m_dot)/area*2*geo.delta_r/visc_up
    # Flank leakage Reynolds number
    elif leakage_type == 'fl':
        Re = abs(m_dot)/area*2*geo.delta_f/visc_up
    # Raise error if input are not coherent
    else: 
        raise Exception('Invalid inputs')
        
        
    return (m_dot, h, Re)  

    


    
    
    
    
    
    
    
    
    
    
    
    
    
    