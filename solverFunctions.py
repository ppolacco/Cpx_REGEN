# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 11:29:26 2021

@author: nicol
"""


import GeometryFunctions as geoF
from math import pi, exp
import numpy as np
import propertyFunctions as propF
import processFunctions as procF
import pickle
import os.path
import bearing

# property_matrix = loadmat('R1233zd.mat') 
# m_dot_in = 0
# m_dot_out = 0
# x_l_in = None
# x_l_out = None
# x_l = 0
# h_in = 0
# h_out = 0
# Q_dot = 0
# T = 360
# rho = 30
# m_cv = 0.1
# dVdTheta=-10/1e6
# omega = 300

# A 'controlVolume' class is used to represent every chamber within the scroll machine. The object is used to 
# keep track of the evolving properties of the chamber as the crank rotation goes on. The following input are
# required to generate the chamber object in first place:
# - name: identifies the kind of chamber that is generated.
# - geo: object that contains all the geometric properties of the machine.
# - property_matrix: table in which all the thermodynamic properties of the fluid are stored.
# - dict_V: dictionaries in which the coefficients of the polynamial that approximate chambers' volume evolution are specified.
# - dict_dV: dictionaries in which the coefficients of the polynamial that approximate volume derivtives evolution are specified.
# - dict_area: dictionaries in which the coefficients of the polynamial that approximate surfaces evolution are specified.
# - T_in, rho_in, x_l_in: known thermodynamic properties.
class controlVolume:
    
    # Constructor function for the 'controlVolume' class
    def __init__(self, name, geo, property_matrix, dict_V, dict_dV, dict_area, T_in, rho_in, x_l_in, P_out = None, alpha = 1):
        
        self.name = name
        self.alpha = alpha
        self.dict_V = dict_V  
        self.dict_dV = dict_dV 
        self.dict_area = dict_area
        self.theta = 0
        self.existence = []
        
        # The inlet values of temperature and density are directly assigned as initial values of temperature and
        # density for suction chambers and first pair of compression chambers. The volume is evaluated through the
        # polyval function applied to the polynomial specified in 'dict_V'
        if name == "sa" or name =="s1" or name =="s2" or (name == "c1" and alpha == 1) or (name == "c2" and alpha == 1) : 
            self.theta_shift = 0
            (self.T_init, self.rho_init) = (T_in, rho_in)  
            self.V_init = np.polyval(dict_V[name][alpha-1], self.theta + self.theta_shift)
            # The volume of the suction chamber 's1' is initialized as 5e-8 (which is a convential value to represent
            # a null volume)
            if name =='s1': 
                self.V_init = 5e-8
            self.dict_init = {'T' :  [0] , 'rho' :[0] }
        
        # The initial values of temperature and density in the discharge chambers are evaluated by means of the
        # function 'discharge_is_guess'. This function use the properties at machine inlet, the pressure ratio
        # and a first-try value of the isentropic efficiency to guess the initial value of temperature and density
        # inside the discharge chambers.
        elif name == "ddd" or name == 'd1' or name == 'd2' or name == 'dd':
            # A first-attempt value of the isentropic efficiency is specified 
            eta_is = 0.6
            (self.T_init, self.rho_init) = procF.discharge_is_guess(property_matrix, rho_in, T_in, P_out, x_l_in, eta_is)
            # ***???***
            if (geo.theta_d < 3/2*pi and (name == 'ddd' or name =='dd')) or (geo.theta_d >= 3/2*pi and ( name == 'd1' or name == 'd2')):
                self.theta_shift = 2*pi
                self.x_out = x_l_in
            else:
                self.theta_shift = 0
            self.V_init = np.polyval(dict_V[name][alpha-1], self.theta + self.theta_shift)
            self.dict_init = {'T' :  [0] , 'rho' :[0] }
            
        # The initial values of temperature and density in the inner compression chambers are evaluated by means 
        # of the function 'compression_is_guess'. This function use the properties at machine inlet, the pressure ratio
        # and a first-try value of the isentropic efficiency to guess the initial value of temperature and density
        # inside the inner compression chambers.
        elif (name == 'c1' or name == 'c2') and alpha > 1:
            self.theta_shift = 0
            eta_is = 0.6
            (self.T_init, self.rho_init) = procF.compression_is_guess(property_matrix, dict_V[name], self.alpha, rho_in, T_in,  x_l_in, eta_is)
            self.V_init =  np.polyval(dict_V[name][alpha-1], self.theta + self.theta_shift)
            self.dict_init = {'T' :  [0] , 'rho' :[0] }
        else:
            raise Exception('Invalid inputs')
        
        # Other properties are computed by means of the initial values of temperature and density computed thus far.
        self.x_l_init = x_l_in
        self.P_init = propF.propfunction_mixture(property_matrix, self.x_l_init, 'P', T = self.T_init , rho = self.rho_init)   
        self.Q_init = propF.propfunction_mixture(property_matrix, self.x_l_init, 'Q', T = self.T_init , rho = self.rho_init)    
        self.s_init = propF.propfunction_mixture(property_matrix, self.x_l_init, 's', T = self.T_init , rho = self.rho_init) 
        self.h_init = propF.propfunction_mixture(property_matrix, self.x_l_init, 'h', T = self.T_init , rho = self.rho_init)
        self.m_init =  self.rho_init * self.V_init   
        
        # 
        self.h_out = propF.propfunction_mixture(property_matrix, self.x_l_init, 'h', T = self.T_init , rho = self.rho_init)   

        # 
        self.m_dot_dict_init = {'su' : 0., 'ex' : 0., 'fl1' : 0., 'fl2' : 0., 'r1' : 0., 'r2' : 0., 'r3' : 0., 'r4' : 0.}
        self.h_dict_init = {'su' : 0., 'ex' : 0., 'fl1' : 0., 'fl2' : 0., 'r1' : 0., 'r2' : 0., 'r3' : 0., 'r4' : 0.}
        self.x_l_dict_init = {'su' : 0., 'ex' : 0., 'fl1' : 0., 'fl2' : 0., 'r1' : 0., 'r2' : 0., 'r3' : 0., 'r4' : 0.}
        
        # - The 'mh_sum' attribute is used to compute the total amount of energy flowing through the considered chamber
        # - The 'm_sum' attribute is used to compute the total amount of mass flowing through the considered chamber
        # - The 'h_disc' attribute is used to compute the discharge enthalpy of the chamber.
        self.mh_sum = 0
        self.m_sum = 1e-10
        self.h_disc = 0
        
        
    # The 'next_values' method takes as input the property matrix, the crank angle and the step. Within this method the
    # functions that compute mass and temperature derivatives are called. The output of this method are the updated values
    # of all the properties within the considered control volume.
    def next_values(self, property_matrix, omega, theta_step, x_l = 0, Q_dot = 0):
        
        # Only volume is updated if the considered control volume is the suction 'sa' chamber.
        if self.name == 'sa':
            self.theta = self.theta + theta_step
            self.V = np.polyval(self.dict_V[self.name][self.alpha-1], self.theta)  
        else:
            # Computation of the mass derivative
            dmCVdTheta = procF.dmCV_dTheta(omega, self.m_dot_dict)
            # Computation of oil mass fraction derivative
            dxoildTheta = procF.dxoil_dTheta(omega, self.m, np.polyval(self.dict_dV[self.name][self.alpha-1], self.theta + self.theta_shift), self.m_dot_dict, self.x_l_dict,  self.x_l)
            # Computation of temperature derivative.
            dTdTheta = procF.dT_dTheta(property_matrix, omega, self.m, self.T, self.rho, np.polyval(self.dict_dV[self.name][self.alpha-1], self.theta + self.theta_shift),  dmCVdTheta, dxoildTheta, self.m_dot_dict, self.h_dict, self.x_l, Q_dot)
            
            # If chamber's volume is lower than 5e-8 the derivatives are set equal to 0 in order to avoid any
            # numerical issue related with division by very small numbers.
            if self.V <=  5e-8: # avoid error suction chamber 
                dTdTheta = 0
                dmCVdTheta = 0
            
            # Attributes are updated accordingly with the derivatives just computed
            self.theta = self.theta + theta_step
            self.T = self.T + theta_step*dTdTheta
            self.m = self.m + theta_step*dmCVdTheta
            self.x_l = self.x_l + theta_step*dxoildTheta       
            
            # Evaluation of chamber volume
            self.V =  np.polyval(self.dict_V[self.name][self.alpha-1], self.theta + self.theta_shift) 
            # If volume is lower than 5e-8, its value is set to 5e-8 to avoid error in the suction chamber
            if self.V <=  5e-8:
                self.V = 5e-8
            
            # Computation of fluid density within the control volume.
            self.rho = self.m/self.V
            # if self.name == 'dd':
            
            #     print(self.name,self.theta)

            # Computation of other properties based on temperature (computed thanks to energy equation) and
            # density (computed via the definition of density itself).
            self.Q = propF.propfunction_mixture(property_matrix, self.x_l, 'Q', T = self.T , rho = self.rho)
            self.P = propF.propfunction_mixture(property_matrix, self.x_l, 'P', T = self.T , rho = self.rho)
            self.s = propF.propfunction_mixture(property_matrix, self.x_l, 's', T = self.T , rho = self.rho) 
            self.h = propF.propfunction_mixture(property_matrix, self.x_l, 'h', T = self.T , rho = self.rho) 
            
    # The 'record_values' function is used to record the value of the properties within each chamber at every crank 
    # angle step during the last crank angle rotation (after stopping criteria are met).
    def record_values(self,existence):
        
        self.existence.append(existence)
        self.P_vect_final.append(self.P)
        self.T_vect_final.append(self.T)
        self.Q_vect_final.append(self.Q)
        self.s_vect_final.append(self.s)
        self.m_vect_final.append(self.m)
        self.x_l_vect_final.append(self.x_l)
        self.rho_vect_final.append(self.rho)
        self.theta_vect_final.append(self.theta)
        self.V_vect_final.append(self.V)
        # 'self.m_dot_vect_final' is a dictionary composed of different lists, one for each mass flow that occurs within
        # the scroll machine.
        for k in self.m_dot_vect_final: 
            self.m_dot_vect_final[k].append(self.m_dot_dict[k])
    
    
    def discharge(self, old_CV):
        
        self.T = old_CV.T
        self.rho = old_CV.rho
        self.x_l = old_CV.x_l
        self.m = old_CV.m
        self.V = old_CV.V
        self.Q = old_CV.Q
        self.P = old_CV.P
        self.s = old_CV.s
        self.h = old_CV.h
        
        self.theta = old_CV.theta
        
        if old_CV.name == 'ddd':
            self.mh_sum = old_CV.mh_sum 
            self.m_sum = old_CV.m_sum
            self.h_disc = old_CV.h_disc
            self.h_out = old_CV.h_out
        

    # Initialize is a method of the 'controlVolume' class. It does not require any input, 'CV' is an optional input whose default value
    # is set to 'None'
    def initialize(self, CV = None):
    
        # If a CV is passed as input to the function, its properties are attached as 
        # init properties to the 'controlVolume' object
        if CV != None:
            self.theta = 0
            self.T_init = CV.T
            self.rho_init = CV.rho
            self.x_l_init = CV.x_l
            self.P_init = CV.P  
            self.Q_init = CV.Q  
            self.s_init = CV.s
            self.h_init = CV.h
            
            # 'V_init' is evaluated thanks to 'dict_V', where all the information about volumes evolution are stored.
            self.V_init = np.polyval(self.dict_V[self.name][self.alpha-1], self.theta + self.theta_shift)
            
            # The volume of the suction chamber 's1' is always initialized as 5e-8
            if self.name == 's1': 
                self.V_init = 5e-8
                
            # Mass inside the control volume is evaluated as volume times density
            self.m_init = self.V_init*self.rho_init 
           
            # The dictionaries related to mass flows are re-initialized.
            self.m_dot_dict_init = {'su' : 0., 'ex' : 0., 'fl1' : 0., 'fl2' : 0., 'r1' : 0., 'r2' : 0., 'r3' : 0., 'r4' : 0.}
            self.h_dict_init = {'su' : 0., 'ex' : 0., 'fl1' : 0., 'fl2' : 0., 'r1' : 0., 'r2' : 0., 'r3' : 0., 'r4' : 0.}
            self.x_l_dict_init = {'su' : 0., 'ex' : 0., 'fl1' : 0., 'fl2' : 0., 'r1' : 0., 'r2' : 0., 'r3' : 0., 'r4' : 0.}
            
        # Initial properties are attached to the properties attributes.
        self.T = self.T_init
        self.rho = self.rho_init
        self.x_l = self.x_l_init
        self.P = self.P_init
        self.Q = self.Q_init
        self.s = self.s_init
        self.h = self.h_init
        self.m = self.m_init     
        self.V = self.V_init
        
        # Initial properties related to mass flows are attached to the properties attributes.
        self.m_dot_dict = self.m_dot_dict_init
        self.h_dict = self.h_dict_init
        self.x_l_dict = self.x_l_dict_init
        
        # The dictionary containing the initial properties is initialized.
        self.dict_init['T'].append(self.T_init)  
        self.dict_init['rho'].append(self.rho_init)
        
        # If discharge entalpy is different from 0 and the difference between outlet enthalpy and discharge enthalpy
        # is lower than 0.1%, outlet enthalpy is set equal to discharge enthalpy
        if self.h_disc != 0 and (self.h_out - self.h_disc)/self.h_disc > 0.001  : 
            self.h_out = self.h_disc
     
        # Initialization of the attribute in which energy flowing through the control volume is stored
        self.mh_sum = 0
        # Initialization of the attribute in which mass flowing through the control volume is stored
        self.m_sum = 1e-10
        # Initialization of the attribute in which discharge enthalpy is stored
        self.h_disc = 0
        
        # Attributes for the final properties are initialized.
        self.P_vect_final = []
        self.T_vect_final = []
        self.Q_vect_final = []
        self.s_vect_final = []
        self.m_vect_final = []
        self.x_l_vect_final = []
        self.rho_vect_final = []
        self.theta_vect_final = []
        self.V_vect_final = []
        self.m_dot_vect_final = {'su' : [], 'ex' : [], 'fl1' : [], 'fl2' : [], 'r1' : [], 'r2' : [], 'r3' : [], 'r4' : []}
    
        
    # The 'update' method of the 'controlVolume' class computes pressure, quality, entropy and enthalpy 
    # of a given control volume based on its pressure and density.
    def update(self, property_matrix):    
        self.P = propF.propfunction_mixture(property_matrix, self.x_l, 'P', T = self.T , rho = self.rho)
        self.Q = propF.propfunction_mixture(property_matrix, self.x_l, 'Q', T = self.T , rho = self.rho)
        self.s = propF.propfunction_mixture(property_matrix, self.x_l, 's', T = self.T , rho = self.rho) 
        self.h = propF.propfunction_mixture(property_matrix, self.x_l, 'h', T = self.T , rho = self.rho) 
        
    
    # The 'mass_flow_su' method computes the mass flow rate at suction based on the assumption that it behaves
    # as an isentropic nozzle 
    def mass_flow_su(self, property_matrix, geo, rho_in, T_in, x_l_in, theta_step):

        # Suction flow is non-zero if and only if the considered chamber is a suction chamber, otherwise the mass
        # flow rate is set equal to 0.         
        if self.name == 's1' or self.name == 's2':
            # Evaluation of the suction port area at the current crank angle
            A_su = np.polyval(self.dict_area['s'], self.theta)
            # Mass flow rate is computed based on the isentropic nozzle hypothesis
            (m_dot_su, h_su, _)  = procF.isentropic_nozzle(geo, property_matrix, A_su,  self.rho,  self.T, self.x_l, rho2 = rho_in,  T2 = T_in,  x_l2 = x_l_in)
            # Chamber dctionaries are updted
            self.m_dot_dict['su'] = m_dot_su
            self.h_dict['su'] = h_su            
            # Energy flow, mass flow and discharge enthalpy of the current object are updated.
            self.mh_sum = self.mh_sum + m_dot_su*h_su
            self.m_sum = self.m_sum + m_dot_su/(2*pi)*theta_step
            self.h_disc = self.mh_sum/self.m_sum/(2*pi)*theta_step
        else:
            (self.m_dot_dict['su'], self.h_dict['su']) = (0., 0.)
          
        
    # The 'mass_flow_ex' method computes the mass flow rate at discharge based on the assumption that it behaves
    # as an isentropic nozzle     
    def mass_flow_ex(self, property_matrix, geo, P_out, theta_step, Reed_valve):    
 
        # Discharge flow is non-zero if and only if the considered chamber is a discharge chamber, otherwise the mass
        # flow rate is set equal to 0. A simple model of a discharge valve is implemented: the flow is non-zero only
        # if the pressure inside the discharge chamber is higher than the pressure in the discharge pipe.
        if self.name == 'ddd' or self.name == 'dd': 
            A_port = np.polyval(self.dict_area['port'], self.theta)
            # A simple valve model is considered. It is based on a very simple model of a pre-charged valve
            # that closes the discharge port if the pressure in the discharge chamber is lower than a fixed
            # value and opens gradually the port until the pressure reaches the value of the discharge region.
            if Reed_valve == True:
                delta_P = P_out/10
                if self.P < P_out - delta_P:
                    (m_dot_ex, h_ex, _) = (0,0,0)
                elif self.P > P_out: 
                    (m_dot_ex, h_ex, _) = procF.isentropic_nozzle(geo, property_matrix, A_port,  self.rho,  self.T, self.x_l, h2 = self.h_out, P2 = P_out, x_l2 = self.x_out)
                else:
                    A_port = A_port*(self.P - P_out + delta_P)/(delta_P)
                    (m_dot_ex, h_ex, _) = procF.isentropic_nozzle(geo, property_matrix, A_port,  self.rho,  self.T, self.x_l, h2 = self.h_out, P2 = P_out, x_l2 = self.x_out)
            # If reed valve model is not included the discharge port is considered to be fully open at all times
            if Reed_valve == False:
                (m_dot_ex, h_ex, _) = procF.isentropic_nozzle(geo, property_matrix, A_port,  self.rho,  self.T, self.x_l, h2 = self.h_out, P2 = P_out, x_l2 = self.x_out)
            
            self.m_dot_dict['ex'] = m_dot_ex
            self.h_dict['ex'] = h_ex
            

            self.mh_sum = self.mh_sum + m_dot_ex*h_ex
            self.m_sum = self.m_sum + m_dot_ex/(2*pi)*theta_step 
            self.h_disc = self.mh_sum/self.m_sum/(2*pi)*theta_step
        else:
            (self.m_dot_dict['ex'], self.h_dict['ex']) = (0., 0.)
        
    
    # The 'mass_flow_dis_dd' method computes the mass flow rate from chambers 'd1' and 'd2' to chamber 'dd' based
    # on the assumption that it behaves as an isentropic nozzle
    def mass_flow_dis_dd(self, property_matrix, geo, CV_d1):
        
        # The flow from chambers 'd1' and 'd2' to chamber 'dd' is non-zero only if the chamber we are considering is
        # the 'dd' chamber.
        if self.name == 'dd':
            # Evaluation of the port area that is present between chambers 'dd' and 'd1' (or 'd2').
            A_dis = np.polyval(self.dict_area['d'], self.theta)
            # Mass flow rate is computed based on the isentropic nozzle hypothesis
            (m_dot_dis, h_dis, _) = procF.isentropic_nozzle(geo, property_matrix, A_dis,  self.rho,  self.T, self.x_l, rho2 = CV_d1.rho,  T2 = CV_d1.T,  x_l2 = CV_d1.x_l)
            # Mass flow rate is multiplied by 2 because the same flow rate comes from 'd2' and 'd1'
            self.m_dot_dict['fl2'] = 2*m_dot_dis 
            # Chamber dctionaries are updted 
            self.h_dict['fl2'] = h_dis        
            CV_d1.m_dot_dict['fl2'] = - m_dot_dis
            CV_d1.h_dict['fl2'] = h_dis 
    
    
    # The 'leakages' method computes the mass flow rate through internal chambers of the machine.
    # The 'leakages_computation' function passes to the current method the correct objects, based on the kind leakage
    # that has to be modeled.
    def leakages(self, property_matrix, geo, exist, CVs1 = None, CVc1 = None, CVd1 = None, CVddd = None):
        
        A_fl = flank_area(geo)
        
        # The number of existing compression chambers is computed accordingly to the geometry
        N_c = geoF.N_c_calculation(geo,self.theta)
        
        # The correct name of the discharge chamber is extract from the input variable 'exist'
        if exist =='d':
            CVdis = CVd1
            name_dis = 'd1'
            coef = 1
        elif exist == 'ddd':
            CVdis = CVddd
            name_dis = 'ddd'
            coef = 2
            
        # The radial leakage from chamber 's1' to chamber 'sa' is always present and for this reason is computed
        # ahead of everythink else.
        if self.name == 's1' :
                phi_ssa = geoF.get_phi_ssa(geo, self.theta)
                # Radial leakage area is computed through the 'radial_area' function from the 'GeometryFunction' module
                A_sa = radial_area(geo, geo.phi_ie , max( (geo.phi_ie - self.theta), phi_ssa))
                # Mass flow rate is computed using the isentropic nozzle assumption
                (self.m_dot_dict['r1'], self.h_dict['r1'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_sa,  self.rho,  self.T, self.x_l, rho2 = self.rho_init,  T2 = self.T_init,  x_l2 = self.x_l_init, leakage_type = 'r')
                # Leakage Mass flow rate is affected by Reynolds number. Therefore, a correction factor is computed 
                # in the 'corrected_mass_flow' function, accordingly with the empyrical model presented in Ian Bell's thesis. 
                self.m_dot_dict['r1'] = corrected_mass_flow(geo, Re, self.m_dot_dict['r1'], 'r')

        # Different leakage paths exist depending on the number of compression chambers and therefore on the crank angle 'theta'.
        # If zero compression chambers exist, the only existing leakage paths are:
            # - radial leakage 's1' - 'd2'
            # - flank leakage 's1' - 'd2'
            # - radial leakage 's2' - 'discharge'
            # - flank leakage 's2' - 'discharge'
        if N_c == 0:
            if self.name == 's1' :
                # radial leakage 's1' - 'd2'
                A_d2 = radial_area(geo, geo.phi_ie - self.theta, geo.phi_ie - self.theta - pi )
                (self.m_dot_dict['r3'], self.h_dict['r3'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_d2,  self.rho,  self.T, self.x_l, rho2 = CVdis.rho,  T2 = CVdis.T,  x_l2 = CVdis.x_l, leakage_type = 'r')
                self.m_dot_dict['r3'] = corrected_mass_flow(geo, Re, self.m_dot_dict['r3'], 'r')
                (self.m_dot_dict['r2'], self.h_dict['r2']) = (0., 0.)
                (self.m_dot_dict['r4'], self.h_dict['r4']) = (0., 0.)
                # flank leakage 's1' - 'd2'
                (self.m_dot_dict['fl2'], self.h_dict['fl2'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_fl,  self.rho,  self.T, self.x_l, rho2 = CVdis.rho,  T2 = CVdis.T,  x_l2 = CVdis.x_l, leakage_type = 'fl')
                self.m_dot_dict['fl2'] = corrected_mass_flow(geo, Re, self.m_dot_dict['fl2'], 'fl')
                (self.m_dot_dict['fl1'], self.h_dict['fl1']) = (0., 0.)                
            
            elif self.name == name_dis :
                # radial leakage 'discharge' - 's2'
                A_s2 = radial_area(geo, geo.phi_ie - self.theta, geo.phi_ie - self.theta - pi )
                (m_dot, h, Re) = procF.isentropic_nozzle(geo, property_matrix,  A_s2,  self.rho,  self.T, self.x_l, rho2 = CVs1.rho,  T2 = CVs1.T,  x_l2 = CVs1.x_l, leakage_type = 'r')
                m_dot = corrected_mass_flow(geo, Re, m_dot, 'r')
                (self.m_dot_dict['r2'], self.h_dict['r2']) = (coef*m_dot, h)
                (self.m_dot_dict['r1'], self.h_dict['r1']) = (0., 0.)
                (self.m_dot_dict['r2'], self.h_dict['r2']) = (0., 0.)
                (self.m_dot_dict['r4'], self.h_dict['r4']) = (0., 0.)            
                # flank leakage 'discharge' - 's2' 
                (m_dot, h, Re) = procF.isentropic_nozzle(geo, property_matrix,  A_fl,  self.rho,  self.T, self.x_l, rho2 = CVs1.rho,  T2 = CVs1.T,  x_l2 = CVs1.x_l, leakage_type = 'fl')
                m_dot = corrected_mass_flow(geo, Re, m_dot, 'fl')
                (self.m_dot_dict['fl1'], self.h_dict['fl1']) = (coef*m_dot, h)
                (self.m_dot_dict['fl2'], self.h_dict['fl2']) = (0., 0.)               
            
            else:
                (self.m_dot_dict['r1'], self.h_dict['r1']) = (0., 0.)
                (self.m_dot_dict['r2'], self.h_dict['r2']) = (0., 0.)
                (self.m_dot_dict['r3'], self.h_dict['r3']) = (0., 0.)
                (self.m_dot_dict['r4'], self.h_dict['r4']) = (0., 0.)     
                (self.m_dot_dict['fl1'], self.h_dict['fl1']) = (0., 0.)
                (self.m_dot_dict['fl2'], self.h_dict['fl2']) = (0., 0.)    
         
        # If a pair compression chambers exist, the existing leakage paths are:
            # - radial leakage 's1' - 'c21' (same as c11) 
            # - flank leakage 's1' - 'c21' (same as c11) 
            # - radial leakage 'c21' (same as c11) - 's1'
            # - radial leakage 'c21' (same as c11) - 's1'
            # - radial leakage 'c21' (same as c11) - 'sa'
            # - flank leakage 'c21' - 's2'
            # - flank leakage 'c21' - 'discharge'
            # - radial leakage 'discharge' - 'c2' (alpha = N_c)
            # - flank leakages 'discharge' - 'c2' (alpha = N_c)
        elif N_c == 1:
            if self.name == 's1' :
                # radial leakage 's1' - 'c21' (same as c11) 
                phi_ssa = geoF.get_phi_ssa(geo, self.theta)
                A_c21 = radial_area(geo, max( geo.phi_ie - self.theta - pi , phi_ssa) , geo.phi_ie - self.theta - pi)
                (self.m_dot_dict['r3'], self.h_dict['r3'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_c21,  self.rho,  self.T, self.x_l, rho2 = CVc1[1].rho,  T2 = CVc1[1].T,  x_l2 = CVc1[1].x_l, leakage_type = 'r')
                self.m_dot_dict['r3'] = corrected_mass_flow(geo, Re, self.m_dot_dict['r3'], 'r')
                (self.m_dot_dict['r2'], self.h_dict['r2']) = (0., 0.)
                (self.m_dot_dict['r4'], self.h_dict['r4']) = (0., 0.)                
                # flank leakage 's1' - 'c21' (same as c11) 
                (self.m_dot_dict['fl2'], self.h_dict['fl2'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_fl,  self.rho,  self.T, self.x_l, rho2 = CVc1[1].rho,  T2 = CVc1[1].T,  x_l2 = CVc1[1].x_l, leakage_type = 'fl')
                self.m_dot_dict['fl2'] = corrected_mass_flow(geo, Re, self.m_dot_dict['fl2'], 'fl')
                (self.m_dot_dict['fl1'], self.h_dict['fl1']) = (0., 0.)                
                
            elif self.name == 'c1' :  
                # radial leakage 'c21' (same as c11) - 's1'
                phi_ssa = geoF.get_phi_ssa(geo, self.theta)
                A_s2 = radial_area(geo, max( geo.phi_ie - self.theta - pi , phi_ssa) , geo.phi_ie - self.theta - pi)
                (self.m_dot_dict['r2'], self.h_dict['r2'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_s2,  self.rho,  self.T, self.x_l, rho2 = CVs1.rho,  T2 = CVs1.T,  x_l2 = CVs1.x_l, leakage_type = 'r')
                self.m_dot_dict['r2'] = corrected_mass_flow(geo, Re, self.m_dot_dict['r2'], 'r') 
                # radial leakage 'c21' (same as c11) - 's1'
                A_d2 = radial_area(geo, geo.phi_ie - self.theta - 2*pi, geo.phi_is)
                (self.m_dot_dict['r4'], self.h_dict['r4'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_d2,  self.rho,  self.T, self.x_l, rho2 = CVdis.rho,  T2 = CVdis.T,  x_l2 = CVdis.x_l, leakage_type = 'r')
                self.m_dot_dict['r4'] = corrected_mass_flow(geo, Re, self.m_dot_dict['r4'], 'r') 
                # radial leakage 'c21' (same as c11) - 'sa'
                A_sa = radial_area(geo,  max( (geo.phi_ie - self.theta), phi_ssa), phi_ssa)
                (self.m_dot_dict['r1'], self.h_dict['r1'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_sa,  self.rho,  self.T, self.x_l, rho2 = self.rho_init,  T2 = self.T_init,  x_l2 = self.x_l_init, leakage_type = 'r')
                self.m_dot_dict['r1'] = corrected_mass_flow(geo, Re, self.m_dot_dict['r1'], 'r') 
                (self.m_dot_dict['r3'], self.h_dict['r3']) = (0., 0.)                  
                # flank leakage 'c21' - 's2'
                (self.m_dot_dict['fl1'], self.h_dict['fl1'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_fl,  self.rho,  self.T, self.x_l, rho2 = CVs1.rho,  T2 = CVs1.T,  x_l2 = CVs1.x_l, leakage_type = 'fl')
                self.m_dot_dict['fl1'] = corrected_mass_flow(geo, Re, self.m_dot_dict['fl1'], 'fl') 
                # flank leakage 'c21' - 'discharge'
                (self.m_dot_dict['fl2'], self.h_dict['fl2'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_fl,  self.rho,  self.T, self.x_l, rho2 = CVdis.rho,  T2 = CVdis.T,  x_l2 = CVdis.x_l, leakage_type = 'fl')  
                self.m_dot_dict['fl2'] = corrected_mass_flow(geo, Re, self.m_dot_dict['fl2'], 'fl') 
            
            elif self.name == name_dis :
                # radial leakage 'discharge' - 'c2' (alpha = N_c)
                A_c2 = radial_area(geo, geo.phi_ie - self.theta - 2*pi, geo.phi_is)
                (m_dot, h, Re) = procF.isentropic_nozzle(geo, property_matrix,  A_c2,  self.rho,  self.T, self.x_l, rho2 = CVc1[1].rho,  T2 = CVc1[1].T,  x_l2 = CVc1[1].x_l, leakage_type = 'r')
                m_dot = corrected_mass_flow(geo, Re, m_dot, 'r') 
                (self.m_dot_dict['r1'], self.h_dict['r1']) = (coef*m_dot, h)
                (self.m_dot_dict['r2'], self.h_dict['r2']) = (0., 0.)
                (self.m_dot_dict['r3'], self.h_dict['r3']) = (0., 0.)
                (self.m_dot_dict['r4'], self.h_dict['r4']) = (0., 0.) 
                # flank leakages 'discharge' - 'c2' (alpha = N_c)
                (m_dot, h, Re) = procF.isentropic_nozzle(geo, property_matrix,  A_fl,  self.rho,  self.T, self.x_l, rho2 = CVc1[1].rho,  T2 = CVc1[1].T,  x_l2 = CVc1[1].x_l, leakage_type = 'fl')
                m_dot = corrected_mass_flow(geo, Re, m_dot, 'fl') 
                (self.m_dot_dict['fl1'], self.h_dict['fl1']) = (coef*m_dot, h)
                (self.m_dot_dict['fl2'], self.h_dict['fl2']) = (0., 0.)                               
        
            else :
                (self.m_dot_dict['r1'], self.h_dict['r1']) = (0., 0.)
                (self.m_dot_dict['r2'], self.h_dict['r2']) = (0., 0.)
                (self.m_dot_dict['r3'], self.h_dict['r3']) = (0., 0.)
                (self.m_dot_dict['r4'], self.h_dict['r4']) = (0., 0.)
                (self.m_dot_dict['fl1'], self.h_dict['fl1']) = (0., 0.)
                (self.m_dot_dict['fl2'], self.h_dict['fl2']) = (0., 0.)               
                
                
        # If more than 1 pair compression chambers exist, the existing leakage paths are:
            # - radial leakage 's1' - 'c21' (same as c11) 
            # - flank leakage 's1' - 'c21' (same as c11) 
            # - radial leakage 'c1' - 's1'
            # - radial leakage 'c1' - 'sa'
            # - radial leakage 'c1' - 'c2' (alpha = self.alpha + 1)
            # - flank leakage with 'c1' - 's2'
            # - flank leakage 'c1' - 'c2' (alpha = self.alpha + 1)
            # - radial leakage 'c1' - 'c2' (alpha = self.alpha - 1)
            # - radial leakage 'c1' - 'd2' (alpha = N_c)
            # - flank leakages 'c1' - 'c2' (alpha = self.alpha - 1)
            # - flank leakages 'c1' - 'd2'
            # - radial leakages 'c1' - 'c2' (alpha = self.alpha + 1)
            # - radial leakages 'c1' - 'c2' (alpha = self.alpha - 1)
            # - flank leakages 'c1' - 'c2' (alpha = self.alpha - 1)
            # - flank leakages 'c1' - 'c2' (alpha = self.alpha + 1)
            # - flank leakage 'c21' - 'discharge'
            # - flank leakages 'discharge' - 'c2' (alpha = N_c)
        elif N_c > 1:
            if self.name == 's1' :
                # radial leakage 's1' - 'c21' (same as c11)
                phi_ssa = geoF.get_phi_ssa(geo, self.theta)
                A_c21 = radial_area(geo, max( geo.phi_ie - self.theta - pi , phi_ssa) , geo.phi_ie - self.theta - pi)
                (self.m_dot_dict['r3'], self.h_dict['r3'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_c21,  self.rho,  self.T, self.x_l, rho2 = CVc1[1].rho,  T2 = CVc1[1].T,  x_l2 = CVc1[1].x_l, leakage_type = 'r')
                self.m_dot_dict['r3'] = corrected_mass_flow(geo, Re, self.m_dot_dict['r3'], 'r') 
                (self.m_dot_dict['r2'], self.h_dict['r2']) = (0., 0.)
                (self.m_dot_dict['r4'], self.h_dict['r4']) = (0., 0.)    
                # flank leakage 's1' - 'c21' (same as c11)
                (self.m_dot_dict['fl2'], self.h_dict['fl2'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_fl,  self.rho,  self.T, self.x_l, rho2 = CVc1[1].rho,  T2 = CVc1[1].T,  x_l2 = CVc1[1].x_l, leakage_type = 'fl')
                self.m_dot_dict['fl2'] = corrected_mass_flow(geo, Re, self.m_dot_dict['fl2'], 'fl') 
                (self.m_dot_dict['fl1'], self.h_dict['fl1']) = (0., 0.) 

            elif self.name == 'c1' :
                if self.alpha == 1:
                    # radial leakage 'c1' - 's1'
                    phi_ssa = geoF.get_phi_ssa(geo, self.theta)
                    A_s1 = radial_area(geo, max( geo.phi_ie - self.theta - pi , phi_ssa) , geo.phi_ie - self.theta - pi)
                    (self.m_dot_dict['r2'], self.h_dict['r2'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_s1,  self.rho,  self.T, self.x_l, rho2 = CVs1.rho,  T2 = CVs1.T,  x_l2 = CVs1.x_l, leakage_type = 'r')
                    self.m_dot_dict['r2'] = corrected_mass_flow(geo, Re, self.m_dot_dict['r2'], 'r') 
                    # radial leakage 'c1' - 'sa'
                    A_sa = radial_area(geo, max( geo.phi_ie - self.theta, phi_ssa) ,phi_ssa)
                    (self.m_dot_dict['r1'], self.h_dict['r1'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_sa,  self.rho,  self.T, self.x_l, rho2 = CVs1.rho_init,  T2 = CVs1.T_init,  x_l2 = CVs1.x_l_init, leakage_type = 'r')
                    self.m_dot_dict['r1'] = corrected_mass_flow(geo, Re, self.m_dot_dict['r1'], 'r') 
                    # radial leakage 'c1' - 'c2' (alpha = self.alpha + 1)
                    A_c2 = radial_area(geo, geo.phi_ie - self.theta - 2*pi*self.alpha, geo.phi_ie - self.theta - 2*pi*self.alpha - pi)
                    (self.m_dot_dict['r4'], self.h_dict['r4'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_c2,  self.rho,  self.T, self.x_l, rho2 = CVc1[self.alpha + 1].rho,  T2 = CVc1[self.alpha + 1].T,  x_l2 = CVc1[self.alpha + 1].x_l, leakage_type = 'r')
                    self.m_dot_dict['r4'] = corrected_mass_flow(geo, Re, self.m_dot_dict['r4'], 'r') 
                    (self.m_dot_dict['r3'], self.h_dict['r3']) = (0., 0.)               
                    # flank leakage with 'c1' - 's2'
                    (self.m_dot_dict['fl1'], self.h_dict['fl1'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_fl,  self.rho,  self.T, self.x_l, rho2 = CVs1.rho,  T2 = CVs1.T,  x_l2 = CVs1.x_l, leakage_type = 'fl')
                    self.m_dot_dict['fl1'] = corrected_mass_flow(geo, Re, self.m_dot_dict['fl1'], 'fl') 
                    # flank leakage 'c1' - 'c2' (alpha = self.alpha + 1)
                    (self.m_dot_dict['fl2'], self.h_dict['fl2'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_fl,  self.rho,  self.T, self.x_l, rho2 = CVc1[self.alpha + 1].rho,  T2 = CVc1[self.alpha + 1].T,  x_l2 = CVc1[self.alpha + 1].x_l, leakage_type = 'fl')
                    self.m_dot_dict['fl2'] = corrected_mass_flow(geo, Re, self.m_dot_dict['fl2'], 'fl') 
                    
                elif self.alpha == N_c: 
                    # radial leakage 'c1' - 'c2' (alpha = self.alpha - 1)
                    A_c2 = radial_area(geo, geo.phi_ie - self.theta - 2*pi*(self.alpha - 1), geo.phi_ie - self.theta - 2*pi*(self.alpha - 1) - pi)
                    (self.m_dot_dict['r1'], self.h_dict['r1'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_c2,  self.rho,  self.T, self.x_l, rho2 = CVc1[self.alpha - 1].rho,  T2 = CVc1[self.alpha - 1].T,  x_l2 = CVc1[self.alpha - 1].x_l, leakage_type = 'r')
                    self.m_dot_dict['r1'] = corrected_mass_flow(geo, Re, self.m_dot_dict['r1'], 'r') 
                    # radial leakage 'c1' - 'd2' (alpha = N_c)
                    A_d2 = radial_area(geo, geo.phi_ie - self.theta - 2*pi*N_c, geo.phi_is)
                    (self.m_dot_dict['r4'], self.h_dict['r4'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_d2,  self.rho,  self.T, self.x_l, rho2 = CVdis.rho,  T2 = CVdis.T,  x_l2 = CVdis.x_l, leakage_type = 'r')
                    self.m_dot_dict['r4'] = corrected_mass_flow(geo, Re, self.m_dot_dict['r4'], 'r') 
                    (self.m_dot_dict['r2'], self.h_dict['r2']) = (0., 0.)   
                    (self.m_dot_dict['r3'], self.h_dict['r3']) = (0., 0.)                    
                    # flank leakages 'c1' - 'c2' (alpha = self.alpha - 1)
                    (self.m_dot_dict['fl1'], self.h_dict['fl1'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_fl,  self.rho,  self.T, self.x_l, rho2 = CVc1[self.alpha - 1].rho,  T2 = CVc1[self.alpha - 1].T,  x_l2 = CVc1[self.alpha - 1].x_l, leakage_type = 'fl')
                    self.m_dot_dict['fl1'] = corrected_mass_flow(geo, Re, self.m_dot_dict['fl1'], 'fl') 
                    # flank leakages 'c1' - 'd2'
                    (self.m_dot_dict['fl2'], self.h_dict['fl2'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_fl,  self.rho,  self.T, self.x_l, rho2 = CVdis.rho,  T2 = CVdis.T,  x_l2 = CVdis.x_l, leakage_type = 'fl')
                    self.m_dot_dict['fl2'] = corrected_mass_flow(geo, Re, self.m_dot_dict['fl2'], 'fl') 
                
                else : 
                    # radial leakages 'c1' - 'c2' (alpha = self.alpha + 1)
                    A_c2p1 = radial_area(geo, geo.phi_ie - self.theta - 2*pi*self.alpha, geo.phi_ie - self.theta - 2*pi*self.alpha - pi)
                    (self.m_dot_dict['r4'], self.h_dict['r4'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_c2p1,  self.rho,  self.T, self.x_l, rho2 = CVc1[self.alpha + 1].rho,  T2 = CVc1[self.alpha + 1].T,  x_l2 = CVc1[self.alpha + 1].x_l, leakage_type = 'r')
                    self.m_dot_dict['r4'] = corrected_mass_flow(geo, Re, self.m_dot_dict['r4'], 'r') 
                    # radial leakages 'c1' - 'c2' (alpha = self.alpha - 1)
                    A_c2m1 = radial_area(geo, geo.phi_ie - self.theta - 2*pi*(self.alpha - 1), geo.phi_ie - self.theta - 2*pi*(self.alpha - 1) - pi)
                    (self.m_dot_dict['r1'], self.h_dict['r1'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_c2m1,  self.rho,  self.T, self.x_l, rho2 = CVc1[self.alpha - 1].rho,  T2 = CVc1[self.alpha - 1].T,  x_l2 = CVc1[self.alpha - 1].x_l, leakage_type = 'r')
                    self.m_dot_dict['r1'] = corrected_mass_flow(geo, Re, self.m_dot_dict['r1'], 'r') 
                    (self.m_dot_dict['r2'], self.h_dict['r2']) = (0., 0.)   
                    (self.m_dot_dict['r3'], self.h_dict['r3']) = (0., 0.)        
                    # flank leakages 'c1' - 'c2' (alpha = self.alpha - 1)
                    (self.m_dot_dict['fl1'], self.h_dict['fl1'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_fl,  self.rho,  self.T, self.x_l, rho2 = CVc1[self.alpha - 1].rho,  T2 = CVc1[self.alpha - 1].T,  x_l2 = CVc1[self.alpha - 1].x_l, leakage_type = 'fl')
                    self.m_dot_dict['fl1'] = corrected_mass_flow(geo, Re, self.m_dot_dict['fl1'], 'fl') 
                    # flank leakages 'c1' - 'c2' (alpha = self.alpha + 1)
                    (self.m_dot_dict['fl2'], self.h_dict['fl2'], Re) = procF.isentropic_nozzle(geo, property_matrix,  A_fl,  self.rho,  self.T, self.x_l, rho2 = CVc1[self.alpha + 1].rho,  T2 = CVc1[self.alpha + 1].T,  x_l2 = CVc1[self.alpha + 1].x_l, leakage_type = 'fl')
                    self.m_dot_dict['fl2'] = corrected_mass_flow(geo, Re, self.m_dot_dict['fl2'], 'fl') 
                    
            elif self.name == name_dis :
                # radial leakages 'discharge' - 'c2' (alpha = N_c)
                A_c2 = radial_area(geo, geo.phi_ie - self.theta - 2*pi*N_c, geo.phi_is)
                (m_dot, h, Re) = procF.isentropic_nozzle(geo, property_matrix,  A_c2,  self.rho,  self.T, self.x_l, rho2 = CVc1[N_c].rho,  T2 = CVc1[N_c].T,  x_l2 = CVc1[N_c].x_l, leakage_type = 'r')
                m_dot = corrected_mass_flow(geo, Re, m_dot, 'r') 
                (self.m_dot_dict['r1'], self.h_dict['r1']) = (coef*m_dot, h)
                (self.m_dot_dict['r2'], self.h_dict['r2']) = (0., 0.)
                (self.m_dot_dict['r3'], self.h_dict['r3']) = (0., 0.)
                (self.m_dot_dict['r4'], self.h_dict['r4']) = (0., 0.)              
                # flank leakages 'discharge' - 'c2' (alpha = N_c)
                (m_dot, h, Re) = procF.isentropic_nozzle(geo, property_matrix,  A_fl,  self.rho,  self.T, self.x_l, rho2 = CVc1[N_c].rho,  T2 = CVc1[N_c].T,  x_l2 = CVc1[N_c].x_l, leakage_type = 'fl')
                m_dot = corrected_mass_flow(geo, Re, m_dot, 'fl') 
                (self.m_dot_dict['fl1'], self.h_dict['fl1']) = (coef*m_dot, h)
                (self.m_dot_dict['fl2'], self.h_dict['fl2']) = (0., 0.)      
        
            else :
                (self.m_dot_dict['r2'], self.h_dict['r2']) = (0., 0.)
                (self.m_dot_dict['r2'], self.h_dict['r2']) = (0., 0.)
                (self.m_dot_dict['r3'], self.h_dict['r3']) = (0., 0.)
                (self.m_dot_dict['r4'], self.h_dict['r4']) = (0., 0.) 
                (self.m_dot_dict['fl1'], self.h_dict['fl1']) = (0., 0.)
                (self.m_dot_dict['fl2'], self.h_dict['fl2']) = (0., 0.)
                


# The isentropic nozzle model effectively captures the compressibility effects, although it does not capture the frictional 
# effects. The ratio of the isentropic nozzle mass flow rate prediction to that of the detailed model can be correlated
# with the Reynolds number. The methodology and the equations used for this correlation are thouroughly explained
# in appendix C of Ian Bell's thesis.
# The 'corrected_mass_flow' function takes as input the machine geometry, Reynolds number and leakage type. It provides
# as output the ratio of the isentropic nozzle mass flow rate prediction to that of the more detailed model.
def corrected_mass_flow(geo, Re, m_dot, leakage_type):
    
    # Non-dimensionalization parameters
    delta_0 = 10e-6
    L_0 = 0.005
    
    # Different empyrical coefficient are used if the leakage flow is radial or flank
    if leakage_type == 'r':
        a0 =  2.59321070e+4
        a1 =  9.14825434e-1
        a2 = -1.77588568e+2
        a3 = -2.37052788e-1
        a4 = -1.72347611e+05
        a5 = -1.20687600e+01
        a6 = -1.28861161e-02
        a7 = -1.51202604e+02
        a8 = -9.99674458e-01
        a9 =  1.61435039e-02
        a10 = 8.25533457e-01
        Re_star = 5.24358195e+03
        L = geo.t_s*1e-3
        delta = geo.delta_r     
    elif leakage_type == 'fl':   
        a0 = -2.63970396e+00
        a1 = -5.67164431e-01
        a2 = 8.36554999e-01
        a3 = 8.10567168e-01
        a4 = 6.17402826e+03
        a5 = -7.60907962e+00
        a6 = -5.10200923e-01
        a7 = -1.20517483e+03
        a8 = -1.02938914e+00
        a9 = 6.89497786e-01
        a10 = 1.09607735e+00
        Re_star = 8.26167178e+02
        L = geo.r_o*1e-3
        delta = geo.delta_f
    else :
        raise Exception('Invalid inputs')
    
    # If Reynolds number is equal to 0 no correction is applied to the mass flow
    if Re == 0:
        M = 1
    else:
        zeta = 1/(1 + exp(-0.01*(Re - Re_star)))   
        M = a0*(L/L_0)**a1 / (a2*(delta/delta_0) + a3) * (zeta*(a4*Re**a5 + a6) + (1-zeta)*(a7*Re**a8 + a9)) + a10
    
    return m_dot/M    
        
def radial_area(geo, phi_max, phi_min): 
    
    phi_0 = (geo.phi_i0 + geo.phi_o0)/2
    A_radial = geo.delta_r*geo.r_b/1e3* (0.5*(phi_max**2 - phi_min**2) - phi_0*(phi_max - phi_min))
    return A_radial

def flank_area(geo):
    A_flank = geo.h_s/1e3*geo.delta_f
    
    return A_flank

# The 'merging_process' function takes as input the discharge chambers objects 'CV_ddd', 'CV_dd', 'CV_d1' and
# 'CV_d2'. Its function is to compute the properties of the resulting chamber after the discharge chambers 
# 'd1', 'd2' and 'dd' have merged into chamber 'ddd'.
def merging_process(property_matrix, CV_ddd, CV_dd, CV_d1, CV_d2):
    
    # Energy of the 'ddd' chamber is set equal to the energy of the 'dd' chamber.
    CV_ddd.mh_sum = CV_dd.mh_sum 
    # Mass of the 'ddd' chamber is set equal to the mass of the 'dd' chamber. 
    CV_ddd.m_sum = CV_dd.m_sum
    # Discharge enthalpy of chamber 'ddd' is sert equal to discharge enthalpy of chamber 'dd'
    CV_ddd.h_disc = CV_dd.h_disc
    # Crank angle of chamber 'ddd' is set equal to crank angle of chambr 'dd' 
    CV_ddd.theta = CV_dd.theta
    
    # Volume and mass of chamber 'ddd' is equal to the sum of volumes and masses of chambers 'd1', 'd2' and 'dd'
    CV_ddd.V = CV_dd.V + CV_d1.V + CV_d2.V
    CV_ddd.m = CV_dd.m + CV_d1.m + CV_d2.m
 
    # Pressure of chamber 'ddd' is computed as the weigthed mean of pressure in chambers 'd1', 'd2' and 'dd'
    CV_ddd.P = (CV_dd.V * CV_dd.P + CV_d1.V * CV_d1.P + CV_d2.V * CV_d2.P)/CV_ddd.V
    # Quality of chamber 'ddd' is computed as the weigthed mean of quality in chambers 'd1', 'd2' and 'dd'
    CV_ddd.x_l = (CV_dd.m * CV_dd.x_l + CV_d1.m  * CV_d1.x_l + CV_d2.m  * CV_d2.x_l)/CV_ddd.m
    # Enthalpy of chamber 'ddd' is computed as the weigthed mean of enthalpy in chambers 'd1', 'd2' and 'dd'
    CV_ddd.h = (CV_dd.m * CV_dd.h + CV_d1.m  * CV_d1.h + CV_d2.m  * CV_d2.h)/CV_ddd.m
    #CV_ddd.rho = CV_ddd.m/CV_ddd.V
    
    # Temperature of chamber 'ddd' is computed using enthalpy and pressure
    CV_ddd.T = propF.propfunction_mixture(property_matrix, CV_ddd.x_l, 'T',  h = CV_ddd.h , P = CV_ddd.P)
    #CV_ddd.T = propfunction_mixture(property_matrix, CV_ddd.x_l, 'T',  rho = CV_ddd.rho , P = CV_ddd.P)
    # Temperature of chamber 'ddd' is computed using enthalpy and pressure
    CV_ddd.rho = propF.propfunction_mixture(property_matrix, CV_ddd.x_l, 'rho', h = CV_ddd.h , P = CV_ddd.P)
    
    # Properties of chamber 'ddd' are computed (AGAIN???) using the 'update' method
    CV_ddd.update(property_matrix)
    

# The 'stopping_criteria' function verifies whether convergence is attained at the end of a complete rotation of the
# crank angle. 
def stopping_criteria(geo, CV_c1, CV_ddd, tol)  :   
    
    
    L = len(CV_ddd.dict_init['T']) - 1 
    
    ddd_cond_T = abs(CV_ddd.dict_init['T'][L] \
        - CV_ddd.dict_init['T'][L-1])/CV_ddd.dict_init['T'][L] <= tol
    
    ddd_cond_rho = abs(CV_ddd.dict_init['rho'][L] \
        - CV_ddd.dict_init['rho'][L-1])/CV_ddd.dict_init['rho'][L] <= tol
    
    alpha = 1
    c1_cond_T = {}
    c1_cond_rho = {}
    
    stop = True and ddd_cond_T and ddd_cond_rho
    
    # The size the c1_cond_(T,rho) lists are variable depending on the geometry.
    while alpha <= geo.N_c_max:
        
        c1_cond_T[alpha] = abs(CV_c1[alpha].dict_init['T'][L] \
            - CV_c1[alpha].dict_init['T'][L-1])/CV_c1[alpha].dict_init['T'][L] <= tol
        c1_cond_rho[alpha] = abs(CV_c1[alpha].dict_init['rho'][L] \
            - CV_c1[alpha].dict_init['rho'][L-1])/CV_c1[alpha].dict_init['rho'][L] <= tol    

        stop = stop and  c1_cond_T[alpha] and c1_cond_rho[alpha]
        
        alpha += 1
        
    stop = stop == False
    
    return stop
    
 
    
# The function 'performance_computation' is responsible for computing mass flow rate through the machine, power
# absorbed, isentropic efficiency and volumetric efficiency. It only requires the properties of the 'ddd' and 's1'
# chambers, outlet pressure and mechanical power losses (expressed by means of the torque 'tau_loss').
def performance_computation(geo, property_matrix, CV_ddd, CV_s1, P_out, omega, tau_loss): 
    
    # The mass flow rate coming into the compressor is expressed as the mass flow rate passing through chamber 's1'
    # plus the mass flow rate passing through chamber 's2'. Since these 2 chambers are modeled as identical only the first one
    # is considered and then it is multiplied by 2.
    m_dot_in = CV_s1.m_sum*2
    
    # The mass flow rate coming out of the machine coincide with the mass flow rate through chamber 'ddd'
    m_dot_out = -CV_ddd.m_sum
    
    # Overall mass flow rate is computed as the arithmetic mean mass flow rates in and out of the machine.
    m_dot = (m_dot_in + m_dot_out)/2
    
    # Volumetric efficiency is computed as the ratio between the actual mass flow rate and maximum volumetric
    # efficiency that could be attained by the machine. The maximum mass flow rate is computed as the volume
    # displaced by the machine in one rotation (stored in the 'geo' object as 'V_disp') times the initial 
    # density of the fluid (density at inlet)
    eta_v = m_dot/(omega/2/pi*geo.V_disp*CV_s1.rho_init)
    h_out_is = propF.propfunction_mixture(property_matrix, CV_ddd.x_l, 'h', s = CV_s1.s_init , P = P_out)
    
    # Power loss is computed as torque loss times angular velocity
    W_dot_loss = omega*tau_loss
    # Useful power is computed as mass flow rate times enthalpy increase between inlet and outlet
    W_dot_i = m_dot*(CV_ddd.h_disc - CV_s1.h_init)
    W_dot = (W_dot_i + W_dot_loss)
    
    # Isentropic efficiency is computed as the ratio between the difference of inlet enthalpy and outlet
    # enthalpy under the assumption of isentropic process, and the actual enthalpy difference between inlet and outlet.
    # The actual enthalpy difference between inlet and outlet is computes as the power absorbed by the machine divided by
    # the mass flow rate processed.
    eta_is = ( h_out_is - CV_s1.h_init )*m_dot/W_dot
    
    print('----------------------------------')
    print('m_dot_in:' ,str(m_dot_in), 'kg/s')
    print('m_dot_out:' ,str(m_dot_out), 'kg/s')
    print('eta_is:' ,str(eta_is), '-')
    print('eta_v:' ,str(eta_v), '-')
    print('Mechanical losses:' ,str(W_dot_loss), 'W')
    print('Power consumed:' ,str(W_dot), 'W')
    
    return m_dot, eta_is, eta_v, W_dot_loss, W_dot   

    
# The 'convert_recorded_values' function use the attributes from all the control volumes to generate dictionaries
# that are used within the 'diversePlots' module to generate the plots for the post process tool.
def convert_recorded_values(geo, CV_s1, CV_c1, CV_d1, CV_dd, CV_ddd):
    
    # The 'm_dot_vect_final' dictionary contains the recorded value of the mass flow rate in all the chambers.
    # These values are now converted from lists to numpy arrays.
    for k in CV_s1.m_dot_vect_final: 
        # Suction
        CV_s1.m_dot_vect_final[k] = np.array(CV_s1.m_dot_vect_final[k])  
        # Discharge 1
        CV_d1.m_dot_vect_final[k] = np.array(CV_d1.m_dot_vect_final[k]) 
        # Discharge 'dd'
        CV_dd.m_dot_vect_final[k] = np.array(CV_dd.m_dot_vect_final[k]) 
        # Discharge 'ddd'
        CV_ddd.m_dot_vect_final[k] = np.array(CV_ddd.m_dot_vect_final[k]) 
        # Compression
        alpha = 1
        while alpha <= geo.N_c_max:
            CV_c1[alpha].m_dot_vect_final[k] = np.array(CV_c1[alpha].m_dot_vect_final[k])
            alpha += 1
    
    # Every other attributes of the 's1' chamber is converted into a numpy array
    vect_s1 = {'P' : np.array(CV_s1.P_vect_final), 'T' : np.array(CV_s1.T_vect_final), 
               'Q': np.array(CV_s1.Q_vect_final) , 's' : np.array(CV_s1.s_vect_final), 
               'm' : np.array(CV_s1.m_vect_final), 'x_l' : np.array(CV_s1.x_l_vect_final),
               'rho' : np.array(CV_s1.rho_vect_final),'theta' : np.array(CV_s1.theta_vect_final),
               'm_dot' : CV_s1.m_dot_vect_final, 'V': np.array(CV_s1.V_vect_final)}
    
    # Every other attributes of the 'c1' chambers is converted into a numpy array
    alpha = 1
    vect_c1 = {}
    while alpha <= geo.N_c_max:
        vect_c1[alpha] = {'P' : np.array(CV_c1[alpha].P_vect_final), 'T' : np.array(CV_c1[alpha].T_vect_final), 
                   'Q': np.array(CV_c1[alpha].Q_vect_final) , 's' : np.array(CV_c1[alpha].s_vect_final), 
                   'm' : np.array(CV_c1[alpha].m_vect_final), 'x_l' : np.array(CV_c1[alpha].x_l_vect_final), 
                   'rho' : np.array(CV_c1[alpha].rho_vect_final),'theta' : np.array(CV_c1[alpha].theta_vect_final),
                   'm_dot' : CV_c1[alpha].m_dot_vect_final, 'V': np.array(CV_c1[alpha].V_vect_final)}
        alpha += 1
    
    # Every other attributes of the 'd1' chamber is converted into a numpy array
    vect_d1 = {'P' : np.array(CV_d1.P_vect_final), 'T' : np.array(CV_d1.T_vect_final), 
               'Q': np.array(CV_d1.Q_vect_final) , 's' : np.array(CV_d1.s_vect_final), 
               'm' : np.array(CV_d1.m_vect_final), 'x_l' : np.array(CV_d1.x_l_vect_final), 
               'rho' : np.array(CV_d1.rho_vect_final),'theta' : np.array(CV_d1.theta_vect_final),
               'm_dot' : CV_d1.m_dot_vect_final, 'V': np.array(CV_d1.V_vect_final)}
    
    # Every other attributes of the 'dd' chamber is converted into a numpy array
    vect_dd = {'P' : np.array(CV_dd.P_vect_final), 'T' : np.array(CV_dd.T_vect_final), 
               'Q': np.array(CV_dd.Q_vect_final) , 's' : np.array(CV_dd.s_vect_final), 
               'm' : np.array(CV_dd.m_vect_final), 'x_l' : np.array(CV_dd.x_l_vect_final), 
               'rho' : np.array(CV_dd.rho_vect_final),'theta' : np.array(CV_dd.theta_vect_final),
               'm_dot' : CV_dd.m_dot_vect_final, 'V': np.array(CV_dd.V_vect_final)}
    
    # Every other attributes of the 'ddd' chamber is converted into a numpy array
    vect_ddd = {'P' : np.array(CV_ddd.P_vect_final), 'T' : np.array(CV_ddd.T_vect_final), 
               'Q': np.array(CV_ddd.Q_vect_final) , 's' : np.array(CV_ddd.s_vect_final), 
               'm' : np.array(CV_ddd.m_vect_final), 'x_l' : np.array(CV_ddd.x_l_vect_final), 
               'rho' : np.array(CV_ddd.rho_vect_final),'theta' : np.array(CV_ddd.theta_vect_final),
               'm_dot' : CV_ddd.m_dot_vect_final, 'V': np.array(CV_ddd.V_vect_final)}

    
    # The 'i' index iterates over the crank angle 'theta'
    for i in range(len(vect_s1['theta'])):
        # The items() method returns a view object. The view object contains the key-value pairs of the dictionary, 
        # as tuples in a list.
        for k, val in vect_s1.items():
            # This if statement checks return true at the angles in which the discharge chambers are merged, therefore
            # when the 'ddd' chamber is defined and the 'd1'-'dd' chambers are not.
            if ((geo.theta_d >= geo.theta_m and (vect_s1['theta'][i] > geo.theta_m and vect_s1['theta'][i] <= geo.theta_d))  
                    or (geo.theta_m >= geo.theta_d and (vect_s1['theta'][i] > geo.theta_m or vect_s1['theta'][i] <= geo.theta_d))) and k != 'theta' :
                
                # Since the 'dd'-'d1' chambers are not defined at this value of theta (chambers are merged into 'ddd') their
                # properties are substituted with 'nan' and will not be displayed in the final plot.
                if k == 'm_dot':
                    for n in vect_s1[k]:
                        vect_d1[k][n][i] = np.nan
                        vect_dd[k][n][i] = np.nan
                else:      
                    vect_d1[k][i] = np.nan
                    vect_dd[k][i] = np.nan
                
            # Since the if statements returns false, it means that at this value of theta the 'ddd' chamber is not defined
            # and the 'dd'-'d1' chambers are not merged. The 'ddd' chamber properties are substituted with 'nan'.
            elif  k != 'theta':
                if k == 'm_dot':
                     for n  in vect_d1[k]:
                         vect_ddd[k][n][i] = np.nan
                else:
                    vect_ddd[k][i] = np.nan
            
            # The innermost compression chamber becomes the discharge chamber right after the discharge angle. The 
            # properties of this chamber are substituted with 'nan' and will not be displayed in the final plot
            if vect_s1['theta'][i] > geo.theta_d and k != 'theta' :
                 if k == 'm_dot':
                     for n in vect_s1[k]:
                         vect_c1[geo.N_c_max][k][n][i] = np.nan
                 else:
                     vect_c1[geo.N_c_max][k][i] = np.nan
                     
    return vect_s1, vect_c1, vect_d1, vect_dd, vect_ddd   


# The 'leakage_computation' function wraps up the leakage flow that has to be modeled and passes the correct 
# inputs to the correct function.
def leakage_computation(property_matrix, geo, exist, CV_s1, CV_c1, CV_d1, CV_ddd, leakage):
    
    if leakage == True:
        CV_s1.leakages(property_matrix, geo, exist, CVc1 = CV_c1, CVd1 = CV_d1, CVddd = CV_ddd)
        CV_d1.leakages(property_matrix, geo, exist, CVs1 = CV_s1, CVc1 = CV_c1, CVddd = CV_ddd)
        CV_ddd.leakages(property_matrix, geo, exist, CVs1 = CV_s1, CVc1 = CV_c1, CVddd = CV_ddd)
        alpha = 1
        while alpha <= geoF.N_c_calculation(geo,CV_s1.theta):
            CV_c1[alpha].leakages(property_matrix, geo, exist, CVs1 = CV_s1, CVc1 = CV_c1, CVd1 = CV_d1, CVddd = CV_ddd)
            alpha += 1    

# Computation of dynamic effects exerted on the scroll set.
def force_computation(CV,CV_sa,geo,theta_step,dict_V):    
    # AXIAL FORCE
    # Axial force is initialized
    axial_force = []
    axial_force_aux = 0 
    
    # The length of the lists that contain the themrodynamic properties in each chamber coincide with the number of
    # crank angle step that have been simulated.
    for i in range(len(CV["s1"].P_vect_final)):
        
        # The axial force exerted on each chamber contributes on the total axial force exerted on the scroll set.
        # The first computed (outside the for loop) is the force exerted by the region between the scroll set and
        # the outer wall of the casing
        V_sa = np.polyval(dict_V["sa"][0], theta_step*i)
        axial_force_aux += (CV_sa.P_init - 101325)*(V_sa/(1e-3*geo.h_s))
        # The force exerted by each other chamber is computed in this for loop.
        for element in CV:
            # The axial force is computed as the chamber section area time the differential pressure between chamber and 
            # machine shell
            if element == 'c1':
                for alpha in range(len(CV[element])):
                    axial_force_aux += CV[element][alpha+1].existence[i]*2*(CV[element][alpha+1].P_vect_final[i] - 101325)*(CV[element][alpha+1].V_vect_final[i]/(1e-3*geo.h_s))
            elif element == 's1' or element == 'd1':
                axial_force_aux += CV[element].existence[i]*2*(CV[element].P_vect_final[i] - 101325)*(CV[element].V_vect_final[i]/(1e-3*geo.h_s))
            else:
                axial_force_aux += CV[element].existence[i]*(CV[element].P_vect_final[i] - 101325)*(CV[element].V_vect_final[i]/(1e-3*geo.h_s))
        # The exial force at the i-th theta is stored within the 'axial_force' list. 'axial_force_aux' is re-initialized
        axial_force.append(axial_force_aux)
        axial_force_aux = 0
    
    # RADIAL FORCE & TILTING MOMENT
    # Radial force and tilting moments are initialized
    tangential_force = []
    radial_force = []
    orbiting_force = []
    tilting_moment = []
    force_x_aux = 0
    force_y_aux = 0
    tilting_moment_aux_x = 0
    tilting_moment_aux_y = 0
    
    # - Tangential Force = it is computed as the ratio between the moment exerted by one chamber and its lever arm, computed
    #                  as the norm of the vector (cx,cy)
    # - Tilting moment = it is computed as the axial force exerted on a chamber multiplied by the lever arm of that force
    
    # The length of the lists that contain the themrodynamic properties in each chamber coincide with the number of
    # crank angle step that have been simulated.
    for i in range(len(CV["s1"].P_vect_final)):
        # The geometric characteristic of each chamber are computed at each crank angle step.
        dict_s1 = geoF.volume_force_s1(geo,theta_step*i,True)
        dict_s2 = geoF.volume_force_s2(geo,theta_step*i,True)
        dict_sa = geoF.volume_force_sa(geo,theta_step*i,True)
        dict_c1 = geoF.volume_force_c1(geo,theta_step*i,True)
        dict_c2 = geoF.volume_force_c2(geo,theta_step*i,True)
        dict_d1 = geoF.volume_force_d1(geo,theta_step*i,True)
        dict_d2 = geoF.volume_force_d2(geo,theta_step*i,True)
        dict_dd = geoF.volume_force_dd(geo,theta_step*i,True)
        
        # Forces acting on the region between scroll set and external wall are computed
        force_x_aux += dict_sa["fx_p"]*CV_sa.P_init*1e-6
        force_y_aux += dict_sa["fy_p"]*CV_sa.P_init*1e-6
        # Tilting moment exerted by the fluid in the region between the scroll set and the outer wall of the casing 
        # is computed
        tilting_moment_aux_x += (CV_sa.P_init - 101325)*(dict_sa["V"]/geo.h_s)*dict_sa["cy"]*1e-9
        tilting_moment_aux_y -= (CV_sa.P_init - 101325)*(dict_sa["V"]/geo.h_s)*dict_sa["cx"]*1e-9
        
        for element in CV:    
            # Force and moment exerted by the fluid on compression chambers 
            if element == 'c1':
                for alpha in range(len(CV[element])):
                    if np.isnan(dict_c1[alpha+1]["cy"]) == False:
                        force_x_aux += CV[element][alpha+1].existence[i]*dict_c1[alpha+1]["fx_p"]*CV[element][alpha+1].P_vect_final[i]*1e-6
                        force_x_aux += CV[element][alpha+1].existence[i]*dict_c2[alpha+1]["fx_p"]*CV[element][alpha+1].P_vect_final[i]*1e-6
                        force_y_aux += CV[element][alpha+1].existence[i]*dict_c1[alpha+1]["fy_p"]*CV[element][alpha+1].P_vect_final[i]*1e-6
                        force_y_aux += CV[element][alpha+1].existence[i]*dict_c2[alpha+1]["fy_p"]*CV[element][alpha+1].P_vect_final[i]*1e-6
                        tilting_moment_aux_x += CV[element][alpha+1].existence[i]*(CV[element][alpha+1].P_vect_final[i] - 101325)*(CV[element][alpha+1].V_vect_final[i]*dict_c1[alpha+1]["cy"]/(geo.h_s))
                        tilting_moment_aux_y -= CV[element][alpha+1].existence[i]*(CV[element][alpha+1].P_vect_final[i] - 101325)*(CV[element][alpha+1].V_vect_final[i]*dict_c1[alpha+1]["cx"]/(geo.h_s))
                        tilting_moment_aux_x += CV[element][alpha+1].existence[i]*(CV[element][alpha+1].P_vect_final[i] - 101325)*(CV[element][alpha+1].V_vect_final[i]*dict_c2[alpha+1]["cy"]/(geo.h_s))
                        tilting_moment_aux_y -= CV[element][alpha+1].existence[i]*(CV[element][alpha+1].P_vect_final[i] - 101325)*(CV[element][alpha+1].V_vect_final[i]*dict_c2[alpha+1]["cx"]/(geo.h_s))
            
            # Force and moment exerted by the fluid on suction chambers 
            elif element == 's1':
                force_x_aux += CV[element].existence[i]*dict_s1["fx_p"]*CV[element].P_vect_final[i]*1e-6
                force_x_aux += CV[element].existence[i]*dict_s2["fx_p"]*CV[element].P_vect_final[i]*1e-6
                force_y_aux += CV[element].existence[i]*dict_s1["fy_p"]*CV[element].P_vect_final[i]*1e-6
                force_y_aux += CV[element].existence[i]*dict_s2["fy_p"]*CV[element].P_vect_final[i]*1e-6
                tilting_moment_aux_x += CV[element].existence[i]*(CV[element].P_vect_final[i] - 101325)*(CV[element].V_vect_final[i]/(geo.h_s))*dict_s1["cy"]
                tilting_moment_aux_y -= CV[element].existence[i]*(CV[element].P_vect_final[i] - 101325)*(CV[element].V_vect_final[i]/(geo.h_s))*dict_s1["cx"]         
                tilting_moment_aux_x += CV[element].existence[i]*(CV[element].P_vect_final[i] - 101325)*(CV[element].V_vect_final[i]/(geo.h_s))*dict_s2["cy"]
                tilting_moment_aux_y -= CV[element].existence[i]*(CV[element].P_vect_final[i] - 101325)*(CV[element].V_vect_final[i]/(geo.h_s))*dict_s2["cx"]         
            
            # Force and moment exerted by the fluid on discharge chambers 
            elif element == 'd1':
                force_x_aux += CV[element].existence[i]*dict_d1["fx_p"]*CV[element].P_vect_final[i]*1e-6
                force_x_aux += CV[element].existence[i]*dict_d2["fx_p"]*CV[element].P_vect_final[i]*1e-6
                force_y_aux += CV[element].existence[i]*dict_d1["fy_p"]*CV[element].P_vect_final[i]*1e-6
                force_y_aux += CV[element].existence[i]*dict_d2["fy_p"]*CV[element].P_vect_final[i]*1e-6
                tilting_moment_aux_x += CV[element].existence[i]*(CV[element].P_vect_final[i] - 101325)*(CV[element].V_vect_final[i]/(geo.h_s))*dict_d1["cy"]
                tilting_moment_aux_y -= CV[element].existence[i]*(CV[element].P_vect_final[i] - 101325)*(CV[element].V_vect_final[i]/(geo.h_s))*dict_d1["cx"]                     
                tilting_moment_aux_x += CV[element].existence[i]*(CV[element].P_vect_final[i] - 101325)*(CV[element].V_vect_final[i]/(geo.h_s))*dict_d2["cy"]
                tilting_moment_aux_y -= CV[element].existence[i]*(CV[element].P_vect_final[i] - 101325)*(CV[element].V_vect_final[i]/(geo.h_s))*dict_d2["cx"]                     

            # Force and moment exerted by the fluid on discharge chambers 
            elif element == 'dd':
                force_x_aux += CV[element].existence[i]*dict_dd["fx_p"]*CV[element].P_vect_final[i]*1e-6
                force_y_aux += CV[element].existence[i]*dict_dd["fy_p"]*CV[element].P_vect_final[i]*1e-6
                tilting_moment_aux_x += CV[element].existence[i]*(CV[element].P_vect_final[i] - 101325)*(CV[element].V_vect_final[i]/(geo.h_s))*dict_dd["cy"]
                tilting_moment_aux_y -= CV[element].existence[i]*(CV[element].P_vect_final[i] - 101325)*(CV[element].V_vect_final[i]/(geo.h_s))*dict_dd["cx"]                                    
            
            # Force and moment exerted by the fluid on discharge chambers
            elif element == 'ddd':
                force_x_aux += CV[element].existence[i]*(dict_dd["fx_p"]*CV[element].P_vect_final[i]*1e-6 + dict_d1["fx_p"]*CV[element].P_vect_final[i]*1e-6 + dict_d2["fx_p"]*CV[element].P_vect_final[i]*1e-6)
                force_y_aux += CV[element].existence[i]*(dict_dd["fy_p"]*CV[element].P_vect_final[i]*1e-6 + dict_d1["fy_p"]*CV[element].P_vect_final[i]*1e-6 + dict_d2["fy_p"]*CV[element].P_vect_final[i]*1e-6)
                cx_ddd = (dict_d1['cx']*dict_d1['V'] + dict_d2['cx']*dict_d2['V'] + dict_dd['cx']*dict_dd['V'])/(dict_d1['V'] + dict_d2['V'] + dict_dd['V'])
                cy_ddd = (dict_d1['cy']*dict_d1['V'] + dict_d2['cy']*dict_d2['V'] + dict_dd['cy']*dict_dd['V'])/(dict_d1['V'] + dict_d2['V'] + dict_dd['V'])
                tilting_moment_aux_x += CV[element].existence[i]*(CV[element].P_vect_final[i] - 101325)*(CV[element].V_vect_final[i]/(geo.h_s))*cy_ddd
                tilting_moment_aux_y -= CV[element].existence[i]*(CV[element].P_vect_final[i] - 101325)*(CV[element].V_vect_final[i]/(geo.h_s))*cx_ddd                     
        
        # Radial force is computed by means of trigonometric considerations, given the force acting on the plane
        # and the crank angle 'theta'
        F = np.sqrt(force_x_aux**2 + force_y_aux**2)
        orbiting_force.append(F)
        alpha = np.arctan2(force_y_aux,force_x_aux)
        tangential_force.append(F*np.cos(np.pi/2 + (geo.phi_ie - np.pi/2 - theta_step*i) - alpha))
        radial_force.append(F*np.sin(np.pi/2 + (geo.phi_ie - np.pi/2 - theta_step*i) - alpha))
        
        # The Forces and moments acting at the i-th theta is stored within the 'tangential_force' list.
        tilting_moment.append(np.sqrt(tilting_moment_aux_x**2 + tilting_moment_aux_y**2))
        
        # Auxiliary variables are re-initialized
        tilting_moment_aux_x = 0
        tilting_moment_aux_y = 0
        
    # Force are converted from string to array
    axial_force = np.array(axial_force)
    tangential_force = np.array(tangential_force)
    radial_force = np.array(radial_force)
    orbiting_force = np.array(orbiting_force)
        
    return axial_force, tangential_force, radial_force, orbiting_force, tilting_moment


# The 'full_resolution' function is called by the main script. It uses as input the geometry of the machine ('geo', 'dictV',
# 'dict_dV', 'dict_area') and the known properties ('T_in', 'rho_in', 'x_l_in', 'P_out') to compute mass flow rate,
# isentropic efficiency, volumetric efficiency and absorbed power.
def full_resolution(geo, property_matrix, dict_V, dict_dV, dict_area, T_in, rho_in, x_l_in, P_out, tol_glob, 
                    tol_pressure, theta_step_base, omega, tau_loss, leakage, heat_transfer, Reed_valve, bearing_model):
    
    # Suction chambers control volumes are generated
    print('----------------------------------')
    print('Starting new simulation')
    CV_s1 = controlVolume("s1", geo, property_matrix, dict_V, dict_dV, dict_area, T_in, rho_in, x_l_in)
    CV_sa = controlVolume("sa", geo, property_matrix, dict_V, dict_dV, dict_area, T_in, rho_in, x_l_in)
    
    # Initialization of the properties within the suction chambers control volumes
    CV_s1.initialize()
    CV_sa.initialize()
    
    # Generation of the compression chambers control volumes. In general, there are more than a pair of compression
    # chambers in a scroll machine, depending on the geometry parameters that are specfied. The maximum number
    # of compression chamber is specified in the 'geo' object as the attribute 'N_c_max'.
    alpha = 1
    CV_c1 = {}
    # Compression chambers are generated until 'alpha' becomes higher that 'N_c_max'. These chambers are placed in
    # a list and number thanks to the index 'alpha'
    while alpha <= geo.N_c_max:
        # generation of the alpha-th compression chamber
        CV_c1[alpha] = controlVolume("c1", geo, property_matrix, dict_V, dict_dV, dict_area, T_in, rho_in, x_l_in, alpha = alpha) 
        # Initialization of the properties within the alpha-th compression chambers control volumes
        CV_c1[alpha].initialize()
        alpha += 1
    
    # Discharge chambers control volumes are generated
    CV_ddd = controlVolume("ddd", geo, property_matrix, dict_V, dict_dV, dict_area, T_in, rho_in, x_l_in, P_out) 
    CV_d1 = controlVolume("d1", geo, property_matrix, dict_V, dict_dV, dict_area, T_in, rho_in, x_l_in, P_out) 
    CV_dd = controlVolume("dd", geo, property_matrix, dict_V, dict_dV, dict_area, T_in, rho_in, x_l_in, P_out)
    
    # Initialization of properties within the discharge chambers control columes
    CV_ddd.initialize()
    CV_d1.initialize()
    CV_dd.initialize()
    
    # The iterative procedure begins here. The stop criterion for the while loop are defined within the 
    # 'stopping_criteria' function.
    while stopping_criteria(geo, CV_c1, CV_ddd, tol_glob):
        
        # Beginning of the rotation
        theta = 0
        # Definition of the step for the crank angles. Initially the 'theta_step_base' is used. It is updated if necessary
        theta_step = theta_step_base
        
        # Inside this while loop the process is simulated from theta = 0 to theta = discharge_angle
        while theta < geo.theta_d - theta_step:
        
            # If leakage computation is turn on, leakages are evaluated. 
            leakage_computation(property_matrix, geo, 'ddd', CV_s1, CV_c1, CV_d1, CV_ddd, leakage)
            
            # Mass flow rates through 'ddd' chamber and 's1' chambers are evaluated using their respective method.
            CV_s1.mass_flow_su(property_matrix, geo, rho_in, T_in, x_l_in, theta_step)
            CV_ddd.mass_flow_ex(property_matrix, geo, P_out, theta_step, Reed_valve) 

            # The 'next_values' method is called for the 's1' chamber. The method updates property attributes based
            # on mass and temperature derivatives computed through mass and energy balances.
            CV_s1.next_values(property_matrix, omega, theta_step)
            alpha = 1
            # The 'next_values' method is called for the compression chambers (in number of N_c_max). The method updates property attributes based
            # on mass and temperature derivatives computed through mass and energy balances.
            while alpha <= geo.N_c_max:
                CV_c1[alpha].next_values(property_matrix, omega, theta_step)
                alpha += 1
            # The 'next_values' method is called for the 'ddd' chamber. The method updates property attributes based
            # on mass and temperature derivatives computed through mass and energy balances.
            CV_ddd.next_values(property_matrix, omega, theta_step)
            
            # Crank angle is updated.
            theta += theta_step
            
        # The 'discharge' method is applied to the discharge chambers 'd1' and 'dd' after the discharge angle is overtaken.
        CV_d1.discharge(CV_c1[geo.N_c_max])
        CV_dd.discharge(CV_ddd)
        
        # The merging process of the discharge chambers is simulated within this while loop. This part of the 
        # simulation goes on until pressure difference of chambers 'd1' and 'dd' is lower than 'tol_pressure',
        # the crank angle must also reach the value (discharge_angle + 0.66pi).
        while abs((CV_d1.P - CV_dd.P)/CV_d1.P) > tol_pressure and theta < geo.theta_d + 2*pi/3:
        
            # If leakage computation is turned on, leakages are evaluated.
            leakage_computation(property_matrix, geo, 'd', CV_s1, CV_c1, CV_d1, CV_ddd, leakage)
          
            # Mass flow rates through 'dd' chamber and 's1' chambers are evaluated using their respective method.
            CV_s1.mass_flow_su(property_matrix, geo, rho_in, T_in, x_l_in, theta_step)
            CV_dd.mass_flow_ex(property_matrix, geo, P_out, theta_step, Reed_valve)
            
            # Mass flow rate through the discharge port is evaluated
            CV_dd.mass_flow_dis_dd(property_matrix, geo, CV_d1)
            
            # The 'next_values' method is called for the 's1' chamber. The method updates property attributes based
            # on mass and temperature derivatives computed through mass and energy balances.
            CV_s1.next_values(property_matrix, omega, theta_step)
            alpha = 1
            # The 'next_values' method is called for the compression chambers (in number of N_c_max). The method updates property attributes based
            # on mass and temperature derivatives computed through mass and energy balances.
            while alpha <= geo.N_c_max-1:
                CV_c1[alpha].next_values(property_matrix, omega, theta_step)
                alpha += 1
            # The 'next_values' method is called for the 'dd' and 'd1' chambers. The method updates property attributes based
            # on mass and temperature derivatives computed through mass and energy balances.
            CV_dd.next_values(property_matrix, omega, theta_step)
            CV_d1.next_values(property_matrix, omega, theta_step)
            
            # Crank angle is updated
            theta += theta_step
            
        # The 'merging_process' function is called right after merging conditions are met. This function computes the 
        # properties of the fluid within the 'ddd' chambers that is created after chambers 'd1', 'd2' and 'dd' have 
        # merged.
        merging_process(property_matrix, CV_ddd, CV_dd, CV_d1, CV_d1)
        # After the merging angle, the shift angle can be set to 0.
        CV_ddd.theta_shift = 0
        
        # In this while loop the remaining part of the crank rotation is simulated
        while theta < 2*pi:
            
            # If leakage computation is turned on, leakages are evaluated.
            leakage_computation(property_matrix, geo, 'ddd', CV_s1, CV_c1, CV_d1, CV_ddd, leakage)
           
            # Mass flow rates through 'ddd' chamber and 's1' chambers are evaluated using their respective method.
            CV_s1.mass_flow_su(property_matrix, geo, rho_in, T_in, x_l_in, theta_step)
            CV_ddd.mass_flow_ex(property_matrix, geo, P_out, theta_step, Reed_valve) 
        
            
            # The 'next_values' method is called for the 's1' chamber. The method updates property attributes based
            # on mass and temperature derivatives computed through mass and energy balances.
            CV_s1.next_values(property_matrix, omega, theta_step)
            alpha = 1
            # The 'next_values' method is called for the compression chambers (in number of N_c_max). The method updates property attributes based
            # on mass and temperature derivatives computed through mass and energy balances.
            while alpha <= geo.N_c_max-1:
                CV_c1[alpha].next_values(property_matrix, omega, theta_step)
                alpha += 1
            # The 'next_values' method is called for the 'ddd' chamber. The method updates property attributes based
            # on mass and temperature derivatives computed through mass and energy balances.
            CV_ddd.next_values(property_matrix, omega, theta_step)
            
            # Update of the crank angle
            theta += theta_step
        
        # Mass flow rate, isentropic efficiency, volumetric efficiency and absorbed power are computed after a
        # complete rotation has been simulated. 
        (m_dot, eta_is, eta_v, W_dot_loss, W_dot) = performance_computation(geo, property_matrix, CV_ddd, CV_s1, P_out, omega, tau_loss)
        
        # The shift angle used for the computation of the 'ddd' volume is set to its initial value of 2*pi.
        CV_ddd.theta_shift = 2*pi
        # The initial values of properties in the 'ddd' chambers for the next rotations are set equal to the value
        # of such properties at theta = 2*pi
        CV_ddd.initialize(CV_ddd)
        
        # The initial values of properties in the innermost compression chambers for the next rotations are set equal to the value
        # of such properties at theta = 2*pi
        alpha = geo.N_c_max
        while alpha > 1:
            CV_c1[alpha].initialize(CV_c1[alpha-1])
            alpha -= 1
        # The initial values of properties in the innermost compression chambers for the next rotations are set equal to
        # the value of these properties in the suction chamber 's1' at theta = 2*pi.
        CV_c1[1].initialize(CV_s1)
        # The initial values of properties in the suction chambers for the next rotations are set equal to
        # the value of these properties in the suction region 'sa' at theta = 2*pi. 
        CV_s1.initialize(CV_sa)
    
        
    # A final rotation is done after the convergence criteria are respected, the simulation could have been stopped
    # here but this rotation is done in order to save all the values of the properties for every chamber as a function
    # of the orbiting angle. The 'record_values' method is called after properties are computed in every chamber.
    
    theta = 0
    theta_step = theta_step_base
    
    while theta < geo.theta_d - theta_step:
    
        leakage_computation(property_matrix, geo, 'ddd', CV_s1, CV_c1, CV_d1, CV_ddd, leakage)
        CV_s1.mass_flow_su(property_matrix, geo, rho_in, T_in, x_l_in, theta_step)
        CV_ddd.mass_flow_ex(property_matrix, geo, P_out, theta_step, Reed_valve)         
        CV_s1.record_values(1)
        
        alpha = 1
        while alpha <= geo.N_c_max:
            CV_c1[alpha].record_values(1)
            alpha += 1
        
        CV_ddd.record_values(1)
        CV_d1.record_values(0)
        CV_dd.record_values(0)
            
        CV_s1.next_values(property_matrix, omega, theta_step)

        alpha = 1
        while alpha <= geo.N_c_max:
            CV_c1[alpha].next_values(property_matrix, omega, theta_step)
            alpha += 1
        CV_ddd.next_values(property_matrix, omega, theta_step)
              
        theta += theta_step
        
    CV_d1.discharge(CV_c1[geo.N_c_max])
    CV_dd.discharge(CV_ddd)
    
    while abs((CV_d1.P - CV_dd.P)/CV_d1.P) > tol_pressure and theta < geo.theta_d + 2*pi/3:
    
        leakage_computation(property_matrix, geo, 'd', CV_s1, CV_c1, CV_d1, CV_ddd, leakage)
    
        CV_s1.mass_flow_su(property_matrix, geo, rho_in, T_in, x_l_in, theta_step)
        CV_dd.mass_flow_ex(property_matrix, geo, P_out, theta_step, Reed_valve)
        CV_dd.mass_flow_dis_dd(property_matrix, geo, CV_d1)
              
    
        CV_s1.record_values(1)
        alpha = 1
        while alpha <= geo.N_c_max:
            if alpha <= geo.N_c_max - 1:
                CV_c1[alpha].record_values(1)
            else:
                CV_c1[alpha].record_values(0)
            alpha += 1
        CV_ddd.record_values(0)
        CV_d1.record_values(1)
        CV_dd.record_values(1)
        
        CV_s1.next_values(property_matrix, omega, theta_step)    
        alpha = 1
        while alpha <= geo.N_c_max-1:
            CV_c1[alpha].next_values(property_matrix, omega, theta_step)  
            alpha += 1
        CV_dd.next_values(property_matrix, omega, theta_step)
        CV_d1.next_values(property_matrix, omega, theta_step)
        
        theta += theta_step
    
    geo.theta_m = theta
    merging_process(property_matrix, CV_ddd, CV_dd, CV_d1, CV_d1)
    CV_ddd.theta_shift = 0
            
    while theta < 2*pi:
        
        leakage_computation(property_matrix, geo, 'ddd', CV_s1, CV_c1, CV_d1, CV_ddd, leakage)
    
        CV_s1.mass_flow_su(property_matrix, geo, rho_in, T_in, x_l_in, theta_step)
        CV_ddd.mass_flow_ex(property_matrix, geo, P_out, theta_step, Reed_valve)  

 
        CV_s1.record_values(1)
        alpha = 1
        while alpha <= geo.N_c_max:
            if alpha <= geo.N_c_max - 1:
                CV_c1[alpha].record_values(1)
            else:
                CV_c1[alpha].record_values(0)
            alpha += 1
        CV_ddd.record_values(1)
        CV_d1.record_values(0)
        CV_dd.record_values(0)
        
        CV_s1.next_values(property_matrix, omega, theta_step)     
        alpha = 1  
        while alpha <= geo.N_c_max-1:
            CV_c1[alpha].next_values(property_matrix, omega, theta_step)
            alpha += 1
        CV_ddd.next_values(property_matrix, omega, theta_step)
        
        theta += theta_step

    (m_dot_final, eta_is_final, eta_v_final, W_dot_loss_final, W_dot_final) =  \
    performance_computation(geo, property_matrix, CV_ddd, CV_s1, P_out, omega, tau_loss)
    
    # 'CV' is a dictionary where all the recorded value are stored.
    CV = {'s1' : CV_s1, 'c1' : CV_c1, 'd1' : CV_d1, 'dd' : CV_dd, 'ddd' : CV_ddd }
    
    # AFter the simulation has been performed it is possible to compute the dynamic effects acting on the scroll set.
    [axial_force, tangential_force, radial_force, orbiting_force, tilting_moment] = force_computation(CV,CV_sa,geo,theta_step,dict_V)
    
    if bearing_model == True:
        # Computation of the mechanical losses occurring in rolling bearings.
        # Forces acting on the rolling bearing are computed based on the machine layout provided by Exoes
        (shaft_force,pulley_force) = bearing.bearing_layout(orbiting_force, a = 35.29, b = 181.54, c = 146.25)
        # Orbiting, shaft and pulley bearing are added to the model
        orbiting_bearing = bearing.rolling_bearing(omega, 0, orbiting_force, d = 30, D = 62, B = 16, K_z = 5.2)
        shaft_bearing = bearing.rolling_bearing(omega, 0, shaft_force, d = 30, D = 62, B = 20, K_z = 5.2)
        pulley_bearing = bearing.rolling_bearing(omega, 0, pulley_force, d = 30, D = 62, B = 20, K_z = 5.2)
        # Computation of resisting moments
        tau_loss_orbiting = 0.001*orbiting_bearing.bearing_tau_model()
        tau_loss_shaft = 0.001*shaft_bearing.bearing_tau_model()
        tau_loss_pulley = 0.001*pulley_bearing.bearing_tau_model()
        
        tau_loss = np.average(tau_loss_orbiting + tau_loss_shaft + tau_loss_pulley)
    
    performance_computation(geo, property_matrix, CV_ddd, CV_s1, P_out, omega, tau_loss)
    print('----------------------------------')
    print('Simulation finished')
    
    return m_dot_final, eta_is_final, eta_v_final, W_dot_final, CV, axial_force, tangential_force, radial_force, tilting_moment, orbiting_force


#%%
# This part of the code is related to the result file managemente.

# 'save_result' function save a pickle file where results are stored.
def save_results(name, results):
    # 'list_file' specifies the precise file directory
    list_file = ['C:\\Users\\Utente1\\Documents\\Tifeo\\Python\\Compressore\\Results' , name ,'.pkl']
    
    # The complete file name is built up
    filename = "".join(list_file)       
    
    # 'file' points at the file in which the pickled object will be written
    with open(filename, "wb") as file:
        # The dump() method of the pickle module in Python, converts a Python object hierarchy into a byte stream. 
        # This process is also called as serilaization.
        pickle.dump(results, file)
  
    
def open_results(name):
    # 'list_file' specifies the precise file directory
    list_file = ['C:\\Users\\Utente1\\Documents\\Tifeo\\Python\\Compressore\\Results' , name ,'.pkl']
    # The complete file name is built up
    filename = "".join(list_file)       
    
    # 'file' points at the file to be opened
    with open(filename, "rb") as file:
        # The load() method of Python pickle module reads the pickled byte stream of one or more python objects 
        # from a file object
        output = pickle.load(file)
     
    return output


def file_name(filename):
    list_path = ['C:\\Users\\Utente1\\Documents\\Tifeo\\Python\\Compressore\\Results\\' , filename ,'.pkl']
    
    path = "".join(list_path) 
    i = 2
    while os.path.isfile(path):
        
        filename_list = [ filename , '_' ,str(i)]
        filename = "".join(filename_list) 
        
        list_path = ['C:\\Users\\Utente1\\Documents\\Tifeo\\Python\\Compressore\\Results\\' , filename ,'.pkl']
        
        path = "".join(list_path) 
        i += 1
        
    return filename



