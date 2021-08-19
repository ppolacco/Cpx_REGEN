# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 16:47:59 2021

@author: Paolo
"""

import numpy as np

class rolling_bearing:
    
    # Constructor method for bearing objects
    def __init__(self,omega,radial_force,axial_force,d,D,B,K_z):
        self.d = d
        self.D = D
        self.B = B
        self.radial_force = radial_force
        self.axial_force = axial_force
        self.d_m = 0.5*(self.d+self.D)
        self.omega = omega*60/(2*np.pi)
        self.K_z = K_z
        self.viscosity = 10 #mm^2/s
        
    # Losses on the bearing are composed of 4 different terms: 
    # - 'rolling friction'
    # - 'sliding friction'
    # - 'seal friction'
    # - 'drag friction' (not clear how to implement it)
    def bearing_tau_model(self):
        # Components of resisting moment are computed
        tau_rr = rolling_bearing.rolling_friction(self)
        tau_sl = rolling_bearing.sliding_friction(self)
        tau_seal = rolling_bearing.seal_friction(self)
        # tau_drag = rolling_bearing.drag_friction()
        tau = tau_rr + tau_sl + tau_seal
        
        return tau
    
    # SKF Model
    def rolling_friction(self):
        
        phi_ish = 1/(1+1.84*1e-9*((self.omega*self.d_m)**1.28)*self.viscosity**0.64)
        phi_rs = 1/np.exp(6*1e-8*self.viscosity*self.omega*2*self.d_m*np.sqrt(self.K_z/(2*(self.D-self.d))))
        R1 = 1.23*10e-6
        G_rr = R1*self.d_m*2.41*self.radial_force**0.31
        tau_rr = phi_ish*phi_rs*G_rr*(self.viscosity*self.omega)**0.6
        
        return tau_rr
    
    # SKF Model
    def sliding_friction(self):
        
        S1 = 0.16
        S2 = 0.0015
        G_sl = S1*self.axial_force*self.d_m**0.9 + S2*self.d_m*self.radial_force
        phi_bl = 1/np.exp(2.6*1e-8*self.d_m*(self.omega*self.viscosity)**1.4)
        mu_bl = 0.12
        mu_EHL = 0.02
        mu_sl = phi_bl*mu_bl + (1 - phi_bl)*mu_EHL 
        tau_sl = mu_sl*G_sl
        
        return tau_sl
    
    # SKF Model
    def seal_friction(self):
        
        KS1 = 0.032
        beta = 2
        KS2 = 50
        tau_seal = KS1*self.B**beta + KS2
        
        return tau_seal
    
    # def drag_friction(self):
        
    #     l_D = 5*K_L*self.B/self.d_m
    #     C_w = 2.789*1e-10*l_D**3 - 2.786*1e-4*l_D**2 + 0.0195*l_D + 0.6439
    #     tau_drag = 4*V_M*K_roll*C_w*self.B*(self.d_m**4)*self.omega**2 + 1.093*1e-7*(self.omega**2)*(self.d_m**3)*R_s*(self.omega*f_t*self.d_m**2/self.viscosity)**-1.379
    
    
def bearing_layout(orbiting_force, a, b, c):
    # Balance of moments acting on the shaft to compute radial force exerted pon shaft bearing and puley bearing
    shaft_force = orbiting_force*b/c
    pulley_force = orbiting_force*a/c
    
    return shaft_force, pulley_force