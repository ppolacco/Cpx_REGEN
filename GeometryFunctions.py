# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 11:55:52 2020

@author: Nicolas 
"""

from math import pi,sin,cos, sqrt, atan2, atan, acos, tan
import numpy as np
from shapely.geometry import Polygon

# The 'GeoParam' class contains all the informations about the machine geometry
class GeoParam:
    """
    Class containing the geometry parameters
    """
    # The class constructor function uses all the characteristics parameters of a scroll machine to generate the 'geo' object.
    def __init__(self,phi_i0, phi_is, phi_ie_fix, phi_ie_orb, phi_o0, phi_os, phi_oe_fix, phi_oe_orb, r_b, h_s, D_wall, x_dis0, y_dis0, r_dis, delta_r, delta_f, L_dis = None, x_dis1 = None, A_dis = None):
        
        # Angles are expressed in radians, radii and heights in mm
        self.phi_i0 = phi_i0
        self.phi_is = phi_is
        self.phi_ie_fix = phi_ie_fix
        self.phi_ie_orb = phi_ie_orb
        self.phi_ie = phi_ie_fix
        self.phi_o0 = phi_o0
        self.phi_os = phi_os
        self.phi_oe_fix = phi_oe_fix
        self.phi_oe_orb = phi_oe_orb
        self.phi_oe = phi_oe_fix
                 
        self.r_b = r_b
        #self.r_o_real = 6.405 
                
        self.h_s = h_s
        
        # Thickness of the scroll
        self.t_s = self.r_b*(self.phi_i0 - self.phi_o0)
        # Orbiting radius
        self.r_o = self.r_b*(pi - self.phi_i0 + self.phi_o0)
        
        # Diameter of the casing
        self.D_wall = D_wall
        
        # The discharge port is an oblonge shape characterized by the following parameters
        self.x_dis0 = x_dis0
        self.y_dis0 = y_dis0
        self.r_dis = r_dis
        self.L_dis = L_dis
        self.x_dis1 = x_dis1
        self.A_dis = A_dis
        
        # Radial and flank clearance (used for leakages computation)
        self.delta_r = delta_r
        self.delta_f = delta_f

    # Function allowing to compute the coordinate of the involute as a function of the angle. The input 'flag' specifies
    # whether the coordinates of a point on the inner curve or outer curve has to be computed 
    def getCoord(self,phi,flag):
        """ Function allowing to compute the coordinate of
        the involute as a function of the angle 
        
        i = inner curve and o = outer curve"""
        
        # Initialization of the coordinates
        x = y = 0
    
        # Coordinates are computed accordingly to what is specified in the 'flag' input
        if flag == "i":
            x = self.r_b*(np.cos(phi) + (phi - self.phi_i0)*np.sin(phi))
            y = self.r_b*(np.sin(phi) - (phi - self.phi_i0)*np.cos(phi))
        elif flag == "o":
            x = self.r_b*(np.cos(phi) + (phi - self.phi_o0)*np.sin(phi))
            y = self.r_b*(np.sin(phi) - (phi - self.phi_o0)*np.cos(phi))
        else :
            print("Error: choose between inner (i) and outer(o) curve")
        
        return (x,y)
    
    # Function allowing to compute the coordinate of the involute of the orbiting scroll as a function 
    # of the angle phi and orbiting angle theta
    def getCoordOrbiting(self,phi,flag,theta):
        """ Function allowing to compute the coordinate of
        the involute of the orbiting scroll as a function of the angle phi and orbiting angle theta
        
        i = inner curve and o = outer curve"""
        
        # Initialization of the coordinates
        x = y = 0
    
        # Coordinates are computed accordingly to what is specified in the 'flag' input
        if flag == "i":
            x = -(self.r_b*(np.cos(phi) + (phi - self.phi_i0)*np.sin(phi))) + self.r_o*cos(self.phi_ie_fix - pi/2 - theta)
            y = -(self.r_b*(np.sin(phi) - (phi - self.phi_i0)*np.cos(phi))) + self.r_o*sin(self.phi_ie_fix - pi/2 - theta)
        elif flag == "o":
            x = -(self.r_b*(np.cos(phi) + (phi - self.phi_o0)*np.sin(phi))) + self.r_o*cos(self.phi_ie_fix - pi/2 - theta)
            y = -(self.r_b*(np.sin(phi) - (phi - self.phi_o0)*np.cos(phi))) + self.r_o*sin(self.phi_ie_fix - pi/2 - theta)
        else :
            print("Error: choose between inner (i) and outer(o) curve")
        return (x,y)

    # The attributes related to the discharge geometry are computes in the 'discharGeo' method. The discharge geometry
    # is usually characterized by radii 'r_a1' and 'r_a2' and the kind of geometry ('arc-line-arc', '2 arc'). Special
    # geometries might require additional input.
    def dischargeGeo(self, r_a1 = 0.0001, r_a2 = 0.0001, type = 'None', r_a0 = 0.0001, angle0 = 0.0001,  dev = 0.0001, dist = 0.0001 ):
        """
        Parameters
        ----------        
        r1 first radius
        r2 second radius
        Type : string for the king of discharge : 2arc (2a) or arc line arc (ala)
        
        The 4 last parameters are based on the special aala geometry of Sanden (see presentation prodivided by exoes)        
        Returns
        -------
        None.
    
        """
        
        # Coordinates of the fixed scroll starting point of the spiral (inner and outer)
        (self.x_fis,self.y_fis) = self.getCoord(self.phi_is,"i")
        (self.x_fos,self.y_fos) = self.getCoord(self.phi_os,"o")
        
        
        delta_x = self.x_fis - self.x_fos
        delta_y = self.y_fis - self.y_fos
        
        # These terms are used in the computation of the maximum value of 'r_a2'.
        omega_a = cos(self.phi_os - self.phi_is) + 1
        omega_b = self.r_o*omega_a - delta_x*(sin(self.phi_os) - sin(self.phi_is)) + delta_y*(cos(self.phi_os) - cos(self.phi_is))
        omega_c = 0.5*(2*delta_x*sin(self.phi_is)*self.r_o - 2*delta_y*cos(self.phi_is)*self.r_o - delta_y**2 - delta_x**2)
        
        # Check if Perfect Meshing Profile (PMP) is achieved
        self.PMP_cond1 = 0
        r_a2max = 0
        # The first statement is related to a discharge geometry for which PMP cannot be achieved
        if self.phi_is < self.phi_os + pi - 1e-8:
            # This is the largest value of 'r_a2' that is possible without scroll-scroll collision
            self.r_a2max = (-omega_b + sqrt(omega_b**2 - 4*omega_a*omega_c))/2/omega_a
        # This instruction is carried out if perfect meshing pre-condition is met
        elif abs(self.phi_is - self.phi_os - pi) <= 1e-8:
            # Unique solution for the 'r_a2' to obtain PMP
            self.r_a2max=-omega_c/omega_b
            self.PMP_cond1 = 1
        # If none of the previous statement is satisfied the geometry is not compatible.
        else:
            print('Error: not compatible starting angles phi_os %.4f phi_is-pi %.4f' %(self.phi_os,self.phi_is-pi))

        # The code prints a warning if the chosen value of the radius 'r_a2' is higher than the maximum value.
        if r_a2 > self.r_a2max:
            print('Waring: r2 is too large, risk of collision between the wraps, max value is: %0.4f' %(r_a2max))
        
        self.r_a2 = r_a2
        
        # Geometry is computed for the '2-arc' discharge geometry
        if type == "2a":
            self.dischargeType = "2a"
            
            print('Value of r2 allowing to reach PMP: %0.4f' %self.r_a2max)
            
            # The radius of arc 1 is uniquely defined once the value of arc 2 si fixed. The two arcs must fulfill the following conditions:
                # - Curves used must be tangent to the involute at inner starting angle
                # - Curves used must be tangent to the involute at outer starting angle2
                # - Curves must be selected such that the scrolls do not contact each other during operation
                # - Curves must pass through the points on the scrolls defined by the inner and outer involute starting angles
            if r_a1 != ((1/2 * delta_y**2 + 1/2 * delta_x**2 + self.r_a2 * delta_x*sin(self.phi_os)-self.r_a2*delta_y*cos(self.phi_os)) \
                       /(self.r_a2*cos(self.phi_os-self.phi_is)+delta_x*sin(self.phi_is)-delta_y*cos(self.phi_is)+self.r_a2)): 
                # If 'r_a1' does not meet the required criteria its value is corrected
                self.r_a1 = ((1/2 * delta_y**2 + 1/2 * delta_x**2 + self.r_a2 * delta_x*sin(self.phi_os)-self.r_a2*delta_y*cos(self.phi_os)) \
                       /(self.r_a2*cos(self.phi_os-self.phi_is)+delta_x*sin(self.phi_is)-delta_y*cos(self.phi_is)+self.r_a2))
                print('Warning: r1 is not correct for 2 arcs geometry: %0.4f mm, value corrected to %0.4f mm' %(r_a1, self.r_a1))

            # PMP is checked
            if abs(self.r_a1 - self.r_a2 - self.r_o)<= 0.5e-1 and self.PMP_cond1 == 1: 
                self.PMP = 1
            else:
                self.PMP = 0  
            
            # Centers of the circles are computed
            self.x_a1 = self.x_fis - sin(self.phi_is)*self.r_a1
            self.y_a1 = self.y_fis + cos(self.phi_is)*self.r_a1
            self.x_a2 = self.x_fos - sin(self.phi_os)*self.r_a2
            self.y_a2 = self.y_fos + cos(self.phi_os)*self.r_a2
            
            # Tangency angles are computed
            self.t_a11 = atan2(self.y_a2-self.y_a1,self.x_a2-self.x_a1)
            self.t_a12 = atan2(self.y_fis-self.y_a1,self.x_fis-self.x_a1)
            self.t_a21 = atan2(self.y_a1-self.y_a2,self.x_a1-self.x_a2)
            self.t_a22 = atan2(self.y_fos-self.y_a2, self.x_fos-self.x_a2)
            
            # In order to ensure that the arcs are traversed in a counter-clockwise fashion, the values for 't_a12' 
            # and 't_a22' are increased by increments of 2 until they are greater than 't_a11' and 't_a21' respectively.
            while self.t_a12 < self.t_a11:
                self.t_a12 += 2*pi
            while self.t_a22 < self.t_a21:
                self.t_a22 += 2*pi
            
            # For '2-arc' geometries L is 0 by definition
            self.L = 0
            
        # Geometry is computed for the 'arc-line-arc' discharge geometry
        elif type == "ala":
            
            self.dischargeType = "ala"

            # The condition to avoid collision between the scroll wraps is ra1 >= ra2 + ro. 
            # If ra1 = ra2 + ro and phi_os = phi_is + pi, the perfect-meshing-profile is obtained            
            r_a1PMP = self.r_o + self.r_a2
            print('Value of r1 allowing to reach PMP: %0.4f'%r_a1PMP)
            
            # If the non-collision condition is not met a warning is printed.        
            if r_a1 < r_a1PMP:
                 print('Warning: value of r_a1 smaller than  %0.4f, risk of collision between the two wraps'%r_a1PMP)
                
            self.r_a1 = r_a1
            
            # Attainment of PMP with the current geometry is checked
            if abs(self.r_a1 - self.r_a2 - self.r_o)<= 0.5e-1 and self.PMP_cond1 == 1: 
                self.PMP = 1
            else:
                self.PMP = 0            
                
            # Centers of the circles are computed
            self.x_a1 = self.x_fis - sin(self.phi_is)*self.r_a1
            self.y_a1 = self.y_fis + cos(self.phi_is)*self.r_a1
            self.x_a2 = self.x_fos - sin(self.phi_os)*self.r_a2
            self.y_a2 = self.y_fos + cos(self.phi_os)*self.r_a2
            
            # Geometry parameters specifically related to the 'arc-line-arc' configuration
            alpha = atan((self.y_a2 - self.y_a1)/(self.x_a2 - self.x_a1))
            d = sqrt((self.x_a2 - self.x_a1)**2 + (self.y_a2 - self.y_a1)**2)
            beta = acos((self.r_a1 + self.r_a2)/d)
                
            # The length of the 'line' part is computed by means of the Pythagorean theorem
            self.L = sqrt(d**2 - (self.r_a1+self.r_a2)**2)
            
            # Coordinates of the tangency points of the discharge geometry are computed
            (self.x_a1t, self.y_a1t) = (self.x_a1 + self.r_a1*cos(beta+alpha), self.y_a1 + self.r_a1*sin(beta+alpha))
            (self.x_a2t, self.y_a2t) = (self.x_a1t + self.L*sin(beta+alpha), self.y_a1t - self.L*cos(beta+alpha))
                 
            # Tangency angles are computed
            self.t_a11 = beta + alpha
            self.t_a12 = atan2(self.y_fis - self.y_a1, self.x_fis - self.x_a1)
            self.t_a21 = atan2(self.y_a2t - self.y_a2, self.x_a2t - self.x_a2)
            self.t_a22 = atan2(self.y_fos - self.y_a2, self.x_fos - self.x_a2)
            
            # In order to ensure that the arcs are traversed in a counter-clockwise fashion, the values for 't_a12' 
            # and 't_a22' are increased by increments of 2 until they are greater than 't_a11' and 't_a21' respectively.
            while self.t_a12 < self.t_a11:
                self.t_a12 += 2*pi
            while self.t_a22 < self.t_a21:
                self.t_a22 += 2*pi
            
            # Slope and intercept of the line segment are computed
            self.m_l = -1/tan(beta+alpha)
            self.b_l = self.y_a1t - self.m_l*self.x_a1t
         
        # The Sanden off-the-shelf machine discharge geometry is not a usual one ('2-arc or 'arc-line'arc'). Indeed, 
        # the retrofitting of the machine carried out by Exoes showed that the geometry is an 'arc-arc-line-arc'
        # with the two first arcs not tangential.
        elif type == "aala":
            self.dischargeType = "aala"
            
            self.PMP = 0
            
            # Definition of the arcs used to define the 'aala' geometry
            self.r_a0 = r_a0
            self.r_a1 = r_a1
            self.r_a2 = r_a2
            
            # Firstly Arc 0 is defined. the centre of arc 0 can be easily figured out by the tangential constraint 
            # found in the link between the inner scroll involute and the arc.
            self.x_a0 = self.x_fis - sin(self.phi_is)*self.r_a0
            self.y_a0 = self.y_fis + cos(self.phi_is)*self.r_a0
            self.x_a2 = self.x_fos - sin(self.phi_os)*self.r_a2
            self.y_a2 = self.y_fos + cos(self.phi_os)*self.r_a2        
            
            # The tangency angle 't_02' is simply the four-quadrant arctangent angle, which is the anti-clockwise 
            # angle in the four-quadrant datum plane.
            self.t_a02 = atan2(self.y_fis - self.y_a0, self.x_fis - self.x_a0) 
            # The tangency angle 't_02' is directly obtained from 't_02' using the definition of 'angle0'
            self.t_a01 = self.t_a02 - angle0 
            # Point (x_t,y_t) is obtained from point (x_0,y_0).
            self.x_a0t = self.x_a0 + self.r_a0*cos(self.t_a01)
            self.y_a0t = self.y_a0 + self.r_a0*sin(self.t_a01)
            
            # The only crucial remaining parameters to figure out in order to have a complete definition of the discharge
            # geometry is the centre of arc 1, which will not be easily computable as being not tangential to the arc 0.
            # Implicit relations can be obtained by means of geometry contraints. These relations are thoroughly explained
            # in section 2.4 of Deliverable 4.2.
            # The resulting equation is a second degree algebraic equation in the unknown 'x_a1'.
            K = self.y_fis + (dist - self.r_a1)/cos(dev) - self.y_a0t
            a = -1 - tan(dev)**2
            b = 2*self.x_a0t + 2*tan(dev)**2*self.x_fis - 2* tan(dev) * K
            c = self.r_a1**2 - self.x_a0t**2 - self.x_fis**2*tan(dev)**2 - K**2 + 2*self.x_fis*tan(dev)*K
            self.x_a1 = (-b - sqrt(b**2 - 4*a*c))/2/a
            # 'y_a1' is obtained by substitution of 'x_a1' into the equation of the line parallel to the line segment
            # of the discharge geometry and passing for the center of arc 1.
            self.y_a1 = (self.x_a1 - self.x_fis)*tan(dev) + self.y_fis + (dist - self.r_a1)/cos(dev)
            
            # The procedure for the remaining part of the discharge geometry is closely related to the one seen
            # for the 'arc-line-arc' geometry.
            self.t_a12 = atan2(self.y_a0t - self.y_a1, self.x_a0t - self.x_a1)
            
            alpha = atan((self.y_a2 - self.y_a1)/(self.x_a2 - self.x_a1))
            d = sqrt((self.x_a2 - self.x_a1)**2 + (self.y_a2 - self.y_a1)**2)
            beta = acos((self.r_a1 + self.r_a2)/d)
            
            self.L = sqrt(d**2 - (self.r_a1+self.r_a2)**2)
                
            (self.x_a1t, self.y_a1t) = (self.x_a1 + self.r_a1*cos(beta+alpha), self.y_a1 + self.r_a1*sin(beta+alpha))
            (self.x_a2t, self.y_a2t) = (self.x_a1t + self.L*sin(beta+alpha), self.y_a1t - self.L*cos(beta+alpha))
            
            self.t_a11 = beta + alpha
            self.t_a21 = atan2(self.y_a2t - self.y_a2, self.x_a2t - self.x_a2)
            self.t_a22 = atan2(self.y_fos - self.y_a2, self.x_fos - self.x_a2)
            
            # In order to ensure that the arcs are traversed in a counter-clockwise fashion, the values for 't_a12' 
            # and 't_a22' are increased by increments of 2 until they are greater than 't_a11' and 't_a21' respectively. 
            while self.t_a02 < self.t_a01:
                self.t_a02 += 2*pi        
            while self.t_a12 < self.t_a11:
                self.t_a12 += 2*pi
            while self.t_a22 < self.t_a21:
                self.t_a22 += 2*pi
                
            # Slope and intercept of the line segment are computed
            self.m_l = -1/tan(beta+alpha)
            self.b_l = self.y_a1t - self.m_l*self.x_a1t
            
        # No other discharge geometries are yet considered
        else:
            print('Error: type of dischage not known:', type)
            
    # The 'discharge_port' method compute the coordinates of the port perimeter, based on the port parameters that
    # are defined within the '__init__' method of the GeoParam class.
    def discharge_port(self,start = 0):
        
        # If L_dis is lest equal to 'None' the discharge port is a circular port
        if self.L_dis == None :
            N_div = 100
            t = np.linspace(start, start + 2*pi, N_div)
            xport = self.x_dis0 + self.r_dis*np.cos(t)
            yport = self.y_dis0 + self.r_dis*np.sin(t)
        
        # If L_dis is specified the discharge port geometry is an oblonge shape proposed by Exoes.
        else:
            N_div = 50
            self.y_dis1 = self.y_dis0 - sqrt(self.L_dis**2 - (self.x_dis0 - self.x_dis1)**2)
            alpha = - atan((self.x_dis0 - self.x_dis1)/(self.y_dis0 - self.y_dis1))
            m = 1/tan(-alpha)
            
            start = alpha + pi
            if start > alpha and start < alpha + pi:
                t0 = np.linspace(start, alpha + pi, N_div)
                x0 = self.x_dis0 + self.r_dis*np.cos(t0)
                y0 = self.y_dis0 + self.r_dis*np.sin(t0)
                
                t1 = np.linspace(alpha + pi, alpha + 2*pi, N_div)
                x1 = self.x_dis1 + self.r_dis*np.cos(t1)
                y1 = self.y_dis1 + self.r_dis*np.sin(t1)

                t2 = np.linspace(alpha, start, N_div)
                x2 = self.x_dis0 + self.r_dis*np.cos(t2)
                y2 = self.y_dis0 + self.r_dis*np.sin(t2)
                
                b0 = y0[-1] - m * x0[-1]
                x_l_0 = np.linspace( x0[-1], x1[0], N_div)
                y_l_0 = m*x_l_0 + b0
                
                b1 = y1[-1] - m * x1[-1]
                x_l_1 = np.linspace(x1[-1], x2[0], N_div)
                y_l_1 = m*x_l_1 + b1
                
                xport = np.concatenate((x0, x_l_0, x1, x_l_1, x2), axis=None)
                yport = np.concatenate((y0, y_l_0, y1, y_l_1, y2), axis=None)
            else : 
                t0 = np.linspace(start, alpha+ 2*pi, N_div)
                x0 = self.x_dis1 + self.r_dis*np.cos(t0)
                y0 = self.y_dis1 + self.r_dis*np.sin(t0)
                
                t1 = np.linspace(alpha, alpha+ pi, N_div)
                x1 = self.x_dis0 + self.r_dis*np.cos(t1)
                y1 = self.y_dis0 + self.r_dis*np.sin(t1)
                
                t2 = np.linspace(alpha + pi, start, N_div)
                x2 = self.x_dis1 + self.r_dis*np.cos(t2)
                y2 = self.y_dis1 + self.r_dis*np.sin(t2)
                
                b0 = y0[-1] - m * x0[-1]
                x_l_0 = np.linspace(x0[-1], x1[0], N_div)
                y_l_0 = m*x_l_0 + b0
                
                b1 = y1[-1] - m * x1[-1]
                x_l_1 = np.linspace(x1[-1], x2[0], N_div)
                y_l_1 = m*x_l_1 + b1
                
                xport = np.concatenate((x0, x_l_0, x1, x_l_1, x2), axis=None)
                yport = np.concatenate((y0, y_l_0, y1, y_l_1, y2), axis=None)
                
        return (xport, yport)

    # The 'getCoordFullGeometry' method of the 'GeoParam' class can be used to obtain the coordinates of the full machine
    # as a function of the crank angle 'theta'
    def getCoordFullGeometry(self, theta):
        """
        Function allowing to get coordinates of scroll contours as a function of the orbiting angle theta
        """
        # The orbiting scroll coordinates are the fixed scroll coordinates mirrored through the origin and offset 
        # by the orbiting radius. The 'theta_m' attribute is defined for sake of simplicity.
        theta_m = self.phi_ie_fix - pi/2 -theta 
        N_div = 1000
        
        
        # Coordinates of the fixed scroll
        # "f" stands for fixed, the fixed involute is first composed with 6 different curves
        
        # Curve 1: inner scroll from end to start
        phi_curve1f = np.linspace(self.phi_ie_fix, self.phi_is, N_div)
        (x_curve1f, y_curve1f) = self.getCoord(phi_curve1f, 'i')
        
        # Arc 0 (if exists)
        if self.dischargeType == "aala":
            phi_curveArc0 = np.linspace(self.t_a02, self.t_a01, N_div)
            x_curveArc0 = self.x_a0 + self.r_a0*np.cos(phi_curveArc0)
            y_curveArc0 = self.y_a0 + self.r_a0*np.sin(phi_curveArc0)
        else:
            x_curveArc0 = x_curve1f[-1]
            y_curveArc0 = y_curve1f[-1]

        phi_curve2f = np.linspace(self.t_a12, self.t_a11, N_div)
        x_curve2f = self.x_a1 + self.r_a1*np.cos(phi_curve2f)
        y_curve2f = self.y_a1 + self.r_a1*np.sin(phi_curve2f)    
        
        #Curve 2: first arc 
        phi_curve2f = np.linspace(self.t_a12, self.t_a11, N_div)
        x_curve2f = self.x_a1 + self.r_a1*np.cos(phi_curve2f)
        y_curve2f = self.y_a1 + self.r_a1*np.sin(phi_curve2f)
        
        #Curve 3: line (if exists) 
        if self.L == 0:
            x_curve3f = x_curve2f[-1]
            y_curve3f = y_curve2f[-1]
        elif self.L > 0:
            x_curve3f = np.linspace(self.x_a1t, self.x_a2t, N_div)
            y_curve3f = self.m_l*x_curve3f + self.b_l
        else:
            print ('Error: negative length of the line (arc-line-arc): L = %0.4f mm' %self.L)
        
        #Curve 4: second arc
        phi_curve4f = np.linspace(self.t_a21, self.t_a22, N_div)
        x_curve4f = self.x_a2 + self.r_a2*np.cos(phi_curve4f)
        y_curve4f = self.y_a2 + self.r_a2*np.sin(phi_curve4f)
        
        #Curve 5: outer scroll from start to end 
        phi_curve5f = np.linspace(self.phi_os, self.phi_oe_fix, N_div)
        (x_curve5f, y_curve5f) = self.getCoord(phi_curve5f, 'o')
        
        #Curve 6 : closing of the involute
        x_curve6f = np.array([x_curve5f[-1],x_curve1f[0]])
        y_curve6f = np.array([y_curve5f[-1],y_curve1f[0]])
       
        # Coordinates of the fixed scroll are stored as a unique numpy array
        self.x_fscroll = np.concatenate((x_curve1f, x_curveArc0, x_curve2f, x_curve3f, x_curve4f, x_curve5f, x_curve6f), axis=None)
        self.y_fscroll = np.concatenate((y_curve1f, y_curveArc0, y_curve2f, y_curve3f, y_curve4f, y_curve5f, y_curve6f), axis=None)
        
        
        # Coordinates of the orbiting scroll
        # "o" stands for orbiting, the orbiting involute is first composed with 6 different curves
        
        # Curve 1: inner scroll from end to start
        phi_curve1o = np.linspace(self.phi_ie_orb, self.phi_is, N_div)
        (x_curve1o, y_curve1o) = self.getCoord(phi_curve1o, 'i')
        
        # Arc 0 (if exists)
        if self.dischargeType == "aala":
            phi_curveArc0o = np.linspace(self.t_a02, self.t_a01, N_div)
            x_curveArc0o = self.x_a0 + self.r_a0*np.cos(phi_curveArc0o)
            y_curveArc0o = self.y_a0 + self.r_a0*np.sin(phi_curveArc0o)
        else:
            x_curveArc0o = x_curve1o[-1]
            y_curveArc0o = y_curve1o[-1]
        
        # Curve 2: first arc 
        phi_curve2o = np.linspace(self.t_a12, self.t_a11, N_div)
        x_curve2o = self.x_a1 + self.r_a1*np.cos(phi_curve2o)
        y_curve2o = self.y_a1 + self.r_a1*np.sin(phi_curve2o)    
        
        # Curve 3: line (if exists) 
        if self.L == 0:
            x_curve3o = x_curve2o[-1]
            y_curve3o = y_curve2o[-1]
        elif self.L > 0:
            x_curve3o = np.linspace(self.x_a1t, self.x_a2t, N_div)
            y_curve3o = self.m_l*x_curve3o + self.b_l
        else:
            print ('Error: negative length of the line (arc-line-arc): L = %0.4f mm' %self.L)
            
        # Curve 4: second arc
        phi_curve4o = np.linspace(self.t_a21, self.t_a22, N_div)
        x_curve4o = self.x_a2 + self.r_a2*np.cos(phi_curve4o)
        y_curve4o = self.y_a2 + self.r_a2*np.sin(phi_curve4o)
        
        # Curve 5: outer scroll from start to end 
        phi_curve5o = np.linspace(self.phi_os, self.phi_oe_orb, N_div)
        (x_curve5o, y_curve5o) = self.getCoord(phi_curve5o, 'o')
        
        # Curve 6 : closing of the involute
        x_curve6o = np.array([x_curve5o[-1],x_curve1o[0]])
        y_curve6o = np.array([y_curve5o[-1],y_curve1o[0]])
        
        # Coordinates of the orbiting scroll are stored as a unique numpy array
        x_curveOrb = np.concatenate((x_curve1o, x_curveArc0o, x_curve2o, x_curve3o, x_curve4o, x_curve5o, x_curve6o), axis=None)
        y_curveOrb = np.concatenate((y_curve1o, y_curveArc0o, y_curve2o, y_curve3o, y_curve4o, y_curve5o, y_curve6o), axis=None)
        
        # The coordinates of the orbiting scroll are mirrored through the origin and offset by the orbiting radius.
        self.x_oscroll = -x_curveOrb + self.r_o*cos(theta_m)
        self.y_oscroll = -y_curveOrb + self.r_o*sin(theta_m)
        
        # 
        self.vect_oxy = np.array([self.x_oscroll,self.y_oscroll])
        self.vect_oxy = self.vect_oxy.transpose()
        
        # Coordinates of the external wall are computed as well.
        phi_wall = np.linspace(0,2*pi,360)
        self.x_wall =  self.D_wall/2*np.cos(phi_wall)
        self.y_wall =  self.D_wall/2*np.sin(phi_wall)
        
    
    # The 'getCoordNorm' method is used to compute the coordinate of a normal vector pointing towards the inside of 
    # the scroll, depending on the volute angle 'phi'. This vector is useful for sake of force computation.
    def getCoordNorm(self,phi,flag):
        """ Function allowing to compute the coordinate of
        the norm of an involute for the calculation of the force, the vector is therefore
        always pointing towards the inside of the scroll
        
        fi = fixed inner curve, fo = fixed outer curve, oi = orbiting inner curve 
        and oo = orbiting outer curve"""
        
        # Vectors are initialized. 
        nx = np.zeros(np.size(phi))
        ny = np.zeros(np.size(phi))
    
        # The computation of the normal vector coordinates is performed throughout the whole range of volute angle 'phi'
        for i in range(len(phi)):
                # Depending on what is specified in the 'flag' input a different zone of the scroll is considered.
                if flag=="fi":
                    nx[i] = -sin(phi[i])
                    ny[i] = cos(phi[i])
                    #nx[i] = sin(phi[i])
                    #ny[i] = -cos(phi[i])
                elif flag=="fo":
                    nx[i] = -sin(phi[i])
                    ny[i] = cos(phi[i])
                elif flag=="oi":
                    nx[i] = -sin(phi[i])
                    ny[i] = +cos(phi[i])
                elif flag=="oo":
                    nx[i] = sin(phi[i])
                    ny[i] = -cos(phi[i])
                else :
                    print("Error: choose between fi, fo ,oi or oo")
        
        return (nx,ny)
   
    # The 'important_values' method computes 5 geometry parameters that are crucial for the simulation:
        # - 'N_c_max': the maximum number of compression chambers that can exist in the machine, given the geometry of the scroll
        # - 'theta_d': the discharge angle is the crank angle at which the discharge region opens up to the compression chambers
        # - 'theta_u': ????
        # - 'V_disp': the volume displaced by the machine in a complete rotation
        # - 'V_ratio': the ratio between the volume displacement and the volume of the innermost compression chamber 
        #              at the discharge angle
    def important_values(self):
        """
        Parameters
        ----------
        self : class self
        Returns
        -------
        discharge, uncontact angle and maximum number of compression chamber
        """
        # if self.dischargeType == 'aala':
        #     self.N_c_max = np.floor((self.phi_ie_fix - self.phi_os - pi + (self.t_a02 - self.t_a01))/2/pi)
        #     self.theta_d = self.phi_ie_fix - self.phi_os -2*pi*self.N_c_max - pi + (self.t_a02 - self.t_a01)
        #     self.theta_u = (-self.t_a11 + self.phi_ie_fix - pi/2)%(2*pi)
        # else:       
        #     self.N_c_max = np.floor((self.phi_ie_fix - self.phi_os - pi )/2/pi)
        #     self.theta_d = self.phi_ie_fix - self.phi_os -2*pi*self.N_c_max - pi 
        #     self.theta_u = (-self.t_a11 + self.phi_ie_fix - pi/2)%(2*pi)
        
        self.N_c_max = int(np.floor((self.phi_ie_fix - self.phi_os - pi )/2/pi))
        self.theta_d = self.phi_ie_fix - self.phi_os -2*pi*self.N_c_max - pi 
        self.theta_u = self.theta_d + 0.1#(-self.t_a11 + self.phi_ie_fix - pi/2)%(2*pi)
        
        # self.V_disp = -2*pi*self.h_s*self.r_b*self.r_o*(3*pi - 2*self.phi_ie_fix + self.phi_i0 + self.phi_o0)
    
        self.V_disp = 2*volume_force_c1(self,0,False)[1]["V"]/1e9
        V_d1_d = volume_force_d1(self, self.theta_d + 1e-5,False)["V"]/1e9
        
        self.V_ratio = self.V_disp/2/V_d1_d
        # theta_uu = (theta_d - geo.t_a11 - geo.t_a12 )%(2*pi)
        # theta_u = (-geo.t_a11 + geo.phi_ie - pi/2)%(2*pi)
    
# The 'N_c_calculation' gives the number of compression chamber as a function of the crank angle 'theta'
def N_c_calculation(geo,theta):
    """
    Give the number of compression chamber as a function of theta
    """
    # if geo.dischargeType == 'aala':
    #     N_c = np.floor((geo.phi_ie_fix - theta - geo.phi_os - pi + (geo.t_a02 - geo.t_a01))/2/pi)
    # else:
    #     N_c = np.floor((geo.phi_ie_fix - theta - geo.phi_os - pi)/2/pi)
    N_c = np.floor((geo.phi_ie_fix - theta -1e-6 - geo.phi_os - pi)/2/pi)
    
    return N_c 

# The 'polygonArea' function computes the area of any polygon given its vertices. 
def polygonArea(x,y):
    
    # The analytical expression for the section of a polygon with (n+1) vertices is implemented.
    N = len(x) - 2
    A = 0
    for i in range(N+1):
        A = A + x[i]*y[i+1] - x[i+1]*y[i] 
    
    return abs(A/2)

# The 'polygonCentoid' function computes the centroid coordinates of any polygon given its vertices. 
def polygonCentroid(x,y):
    
    # The analytical expression for the centroid of a polygon with (n+1) vertices is implemented
    N = len(x) - 2
    sum_x = 0
    sum_y = 0
    for i in range(N+1):
        sum_x = sum_x + (x[i] + x[i+1])*( x[i]*y[i+1] - x[i+1]*y[i]) 
        sum_y = sum_y + (y[i] + y[i+1])*( x[i]*y[i+1] - x[i+1]*y[i]) 
    cx = sum_x/(6*polygonArea(x,y))
    cy = sum_y/(6*polygonArea(x,y))
    
    return cx, cy

# Function allowing to compute the coordinate of a chamber with a polygon
def polygonDefinition(geo, phi_start, phi_int1, phi_int2, phi_end, theta, flag):
    """ Function allowing to compute the coordinate of
    a chamber with a polygon
    
    if flag = "1", then it is a chamber touching the inner side of the fixed scroll
    if flag = "2", then it is a chamber touching the outer side of the fixed scroll
    """
    # The inner and outer side of the scroll are discretized by means of 50 points
    N_div = 50
    phi_outer = np.linspace(phi_start, phi_int1, N_div)
    phi_inner = np.linspace(phi_int2, phi_end, N_div)
    
    # The coordinates of the polygon are computed accordingly with what is specified in the flag input
    if flag == "1":
        # The outer side of the volume touches the inner part of the fixed scroll
        (x_outer, y_outer) = geo.getCoord(phi_outer, 'i')
        (x_inner, y_inner) = geo.getCoordOrbiting(phi_inner, 'o', theta)
    elif flag == "2":
        # The inner side of the volume touches the outer part of the fixed scroll
        (x_outer, y_outer) = geo.getCoord(phi_outer, 'o')
        (x_inner, y_inner) = geo.getCoordOrbiting(phi_inner, 'i', theta)
    
    # The polygon is defined as a concatenation of all the computed points
    x_poly = np.concatenate((x_outer, x_inner), axis=None)
    y_poly = np.concatenate((y_outer, y_inner), axis=None)

    return x_poly, y_poly, N_div

# The 'get_phi_ssa' function computes the angle delimiting the suction chamber and the admission chamber depending
# on the crank angle.
def get_phi_ssa(geo,theta,flag = 'o'):
    """
    Parameters
    ----------
    geo : class
    theta : orbiting angle
    flag : f or o

    Returns
    -------
    phi_ssa : angle delimiting the suction chamber and the admission chamber
    """
    phi_ssa =0
    # Phi_ssa from the origin 
    
    # D = geo.r_o/geo.r_b*((geo.phi_i0 - geo.phi_ie)*sin(theta) - cos(theta) + 1)/(geo.phi_ie - geo.phi_i0)
    # B = 0.5*((geo.phi_o0 - geo.phi_ie + pi) + sqrt((geo.phi_o0 - geo.phi_ie + pi)**2 - 4*D))
    # geo.phi_ssa = geo.phi_ie - pi + B
    
    # Phi_ssa from the tangent line
    if flag == 'o':
        phi_ssa = geo.phi_ie_fix - pi + geo.r_o/geo.r_b*sin(theta)/(geo.phi_ie_fix - pi - geo.phi_o0)
    elif flag =='f':
        phi_ssa = geo.phi_ie_orb - pi + geo.r_o/geo.r_b*sin(theta)/(geo.phi_ie_fix - pi - geo.phi_o0)
    else:
        print("Error flag for phi_ssa")
    
    return phi_ssa

# The function 'volume_force_s1' computes the volume, volume derivative, coordinates of the centroid and orbiting
# moment generated by the chamber 's1'. It returns a dictionary containing these values.
def volume_force_s1(geo,theta,force):
    """
    Parameters
    ----------
    geo : class geo containing the geometric parameters
    theta : Orbiting angle
    
    Returns
    -------
    A dictionnary containing the volume, volume derivative, coordinates of the centroid and orbiting moment 
    generated by the chamber s1
    """
    
    theta_0_volume = 1e-8
    
    # Computation of the angle delimiting suction chambers and admission chamber
    phi_ssa = get_phi_ssa(geo,theta,'o')
    
    # Polygon definition
    (x_polyS1, y_polyS1, N_div) = polygonDefinition(geo, geo.phi_ie_fix - theta, geo.phi_ie_fix , phi_ssa, geo.phi_ie_fix  - theta - pi, theta, "1")
    xy_polyS1 = np.array([x_polyS1,y_polyS1])
    xy_polyS1 = xy_polyS1.transpose()

    # Volume is computd passing the vertices of the polygon to  the 'polygonArea' function and multuplying that area
    # times the height of the scroll.
    V_s1 = polygonArea(x_polyS1,y_polyS1)*geo.h_s
    V_s1 += theta_0_volume
   
    # Derivative is computed by means of the definition of infinitesimal different quotient.
    # The infinitesimal increment is chosen as 2pi/360
    h = 2*pi/360
    phi_ssa_ph = get_phi_ssa(geo, theta + h,'o')
    phi_ssa_mh = get_phi_ssa(geo, theta - h,'o')
    (x_polyS1_ph, y_polyS1_ph, N_div) = polygonDefinition(geo, geo.phi_ie_fix  - theta - h, geo.phi_ie_fix , phi_ssa_ph, geo.phi_ie_fix  - theta - h - pi, theta+h,"1")
    (x_polyS1_mh, y_polyS1_mh, N_div) = polygonDefinition(geo, geo.phi_ie_fix  - theta + h, geo.phi_ie_fix , phi_ssa_mh, geo.phi_ie_fix  - theta + h - pi, theta-h,"1")
    V_s1_ph = polygonArea(x_polyS1_ph,y_polyS1_ph)*geo.h_s
    V_s1_mh = polygonArea(x_polyS1_mh,y_polyS1_mh)*geo.h_s
    dVdTheta_s1 = (V_s1_ph - V_s1_mh)/2/h
    
    poly_s1 = Polygon(xy_polyS1)
    # Centroid is computed by means of the 'centroid' method
    c_s1 = poly_s1.centroid
    cx_s1 = c_s1.x
    cy_s1 = c_s1.y

    
    if force == True:
        # Computation of the forces acting in chamber 's1' is performed
        phi_force_s1 = np.linspace(phi_ssa, geo.phi_ie_fix  - pi - theta, N_div)   
        (nx_s1,ny_s1) = geo.getCoordNorm(phi_force_s1,'oo')
        L = len(x_polyS1)
        pos_fx_s1 = x_polyS1[N_div:L]
        pos_fy_s1 = y_polyS1[N_div:L]
        L_pos = len(pos_fx_s1)
        dA_s1 = geo.h_s*np.sqrt(np.power(pos_fx_s1[1:L_pos] - pos_fx_s1[0:L_pos-1],2) + np.power(pos_fy_s1[1:L_pos] - pos_fy_s1[0:L_pos-1],2))
        dfx_p_s1 = dA_s1*(nx_s1[1:N_div]+nx_s1[0:N_div-1])/2
        dfy_p_s1 = dA_s1*(ny_s1[1:N_div]+ny_s1[0:N_div-1])/2
        fx_p_s1 = np.sum(dfx_p_s1) #mm^2 multipied by 1 Pa is multiplying by 1e-6
        fy_p_s1 = np.sum(dfy_p_s1)
        rOx_s1 = (x_polyS1[N_div:L-1] + x_polyS1[N_div+1:L])/2 - geo.r_o*cos(geo.phi_ie_fix  - pi/2 - theta)
        rOy_s1 = (y_polyS1[N_div:L-1] + y_polyS1[N_div+1:L])/2 - geo.r_o*sin(geo.phi_ie_fix  - pi/2 - theta)
        
        # Moment exerted on the chamber
        MO_p_s1 = 1e-6*np.sum(rOx_s1 * dfy_p_s1 - rOy_s1 * dfx_p_s1) #mm^3 multiply by 1Pa is multiplying by 1
        
        # Directions of the forces acting on the chamber and their application sites are stored in the 'f_vect_s1' array
        f_vect_s1 = np.vstack((pos_fx_s1, pos_fy_s1, nx_s1, ny_s1 ))
        
        # Dictionary containing al the informations about the chamber 's1'
        dict_s1 = {"V" : V_s1, "dVdTheta": dVdTheta_s1, "xy_poly": xy_polyS1, "cx" : cx_s1, "cy" : cy_s1, "MO_p" : MO_p_s1, "f_vect" : f_vect_s1, "fx_p" : fx_p_s1, "fy_p" : fy_p_s1}
    else:
        dict_s1 = {"V" : V_s1, "dVdTheta": dVdTheta_s1, "xy_poly": xy_polyS1, "cx" : cx_s1, "cy" : cy_s1, "MO_p" : [], "f_vect" : [], "fx_p" : [], "fy_p" : []}
       
    return dict_s1

# The function 'volume_force_s2' computes the volume, volume derivative, coordinates of the centroid and orbiting
# moment generated by the chamber 's2'. It returns a dictionary containing these values.
def volume_force_s2(geo,theta,force):   
    """
    Parameters
    ----------
    geo : class geo containing the geometric parameters
    theta : Orbiting angle

    Returns
    -------
    A dictionnary containing the volume, derivative of the volume, coordinate of the polynom, coordinate of the centroid and orbiting moment 
    generated by the chamber s2
    """ 
    
    theta_0_volume = 1e-8
    
    # Computation of the angle delimiting suction chambers and admission chamber
    phi_ssa = get_phi_ssa(geo,theta,'f')
    
    # Polygon definition
    (x_polyS2, y_polyS2, N_div) = polygonDefinition(geo, geo.phi_ie_fix - theta - pi,  phi_ssa, geo.phi_ie_orb, geo.phi_ie_fix - theta, theta, "2")
    xy_polyS2 = np.array([x_polyS2,y_polyS2])
    xy_polyS2 = xy_polyS2.transpose()

    # Volume is computd passing the vertices of the polygon to  the 'polygonArea' function and multuplying that area
    # times the height of the scroll.
    V_s2 = polygonArea(x_polyS2, y_polyS2)*geo.h_s
    V_s2 += theta_0_volume
    
    # Derivative is computed by means of the definition of infinitesimal different quotient.
    # The infinitesimal increment is chosen as 2pi/360
    h = 2*pi/360
    phi_ssa_ph = get_phi_ssa(geo, theta + h,'f')
    phi_ssa_mh = get_phi_ssa(geo, theta - h,'f')
    (x_polyS2_ph, y_polyS2_ph, N_div) = polygonDefinition(geo, geo.phi_ie_fix  - theta - h - pi, phi_ssa_ph , geo.phi_ie_orb, geo.phi_ie_fix  - theta - h , theta+h,"2")
    (x_polyS2_mh, y_polyS2_mh, N_div) = polygonDefinition(geo, geo.phi_ie_fix  - theta + h - pi, phi_ssa_mh , geo.phi_ie_orb, geo.phi_ie_fix  - theta + h , theta-h,"2")
    V_s2_ph = polygonArea(x_polyS2_ph, y_polyS2_ph)*geo.h_s
    V_s2_mh = polygonArea(x_polyS2_mh, y_polyS2_mh)*geo.h_s
    dVdTheta_s2 = (V_s2_ph - V_s2_mh)/2/h
    
    
    poly_s2 = Polygon(xy_polyS2)
    # Centroid is computed by means of the 'centroid' method of shapely
    c_s2 = poly_s2.centroid
    cx_s2 = c_s2.x
    cy_s2 = c_s2.y

    
    if force == True:
        # Computation of the forces acting in chamber 's2' is performed
        phi_force_s2 = np.linspace(geo.phi_ie_orb, geo.phi_ie_fix - theta,  N_div) 
        (nx_s2,ny_s2) = geo.getCoordNorm(phi_force_s2,'oi')
        L = len(x_polyS2)
        pos_fx_s2 = x_polyS2[N_div:L]
        pos_fy_s2 = y_polyS2[N_div:L]
        L_pos = len(x_polyS2[N_div:L])
        dA_s2 = geo.h_s*np.sqrt(np.power(pos_fx_s2[1:L_pos] - pos_fx_s2[0:L_pos-1],2) + np.power(pos_fy_s2[1:L_pos] - pos_fy_s2[0:L_pos-1],2))
        dfx_p_s2 = dA_s2*(nx_s2[1:N_div]+nx_s2[0:N_div-1])/2
        dfy_p_s2 = dA_s2*(ny_s2[1:N_div]+ny_s2[0:N_div-1])/2
        fx_p_s2 = np.sum(dfx_p_s2)
        fy_p_s2 = np.sum(dfy_p_s2)
        rOx_s2 = (x_polyS2[N_div:L-1] + x_polyS2[N_div+1:L])/2 - geo.r_o*cos(geo.phi_ie_fix - pi/2 - theta)
        rOy_s2 = (y_polyS2[N_div:L-1] + y_polyS2[N_div+1:L])/2 - geo.r_o*sin(geo.phi_ie_fix - pi/2 - theta)
        
        # Moment exerted on the chamber
        MO_p_s2 =  1e-6*np.sum(rOx_s2 * dfy_p_s2 - rOy_s2 * dfx_p_s2) 
        
        # Directions of the forces acting on the chamber and their application sites are stored in the 'f_vect_s2' array
        f_vect_s2 = np.vstack((pos_fx_s2, pos_fy_s2, nx_s2, ny_s2 ))
    
        # Dictionary containing al the informations about the chamber 's2'
        dict_s2 = {"V" : V_s2, "dVdTheta" : dVdTheta_s2,"xy_poly": xy_polyS2, "cx" : cx_s2, "cy" : cy_s2, "MO_p" : MO_p_s2, "f_vect" : f_vect_s2, "fx_p" : fx_p_s2, "fy_p" : fy_p_s2}
    else:
        dict_s2 = {"V" : V_s2, "dVdTheta" : dVdTheta_s2,"xy_poly": xy_polyS2, "cx" : cx_s2, "cy" : cy_s2, "MO_p" : [], "f_vect" : [], "fx_p" : [], "fy_p" : []}
        
        
    return dict_s2

# The function 'volume_force_sa' computes the volume, volume derivative, coordinates of the centroid and orbiting
# moment generated by the chamber 'sa'. It returns a dictionary containing these values.
def volume_force_sa(geo,theta,force):    
    """
    Parameters
    ----------
    geo : class geo containing the geometric parameters
    theta : Orbiting angle

    Returns
    -------
    A dictionnary containing the volume, derivative of the volume, coordinate of the polynoms and orbiting moment 
    generated by the chamber sa
    """ 
    # Number of points used to approximate the region comprised between the outer involute of the scroll and the
    # casing wall (namely the 'sa' chamber)
    N_div = 50
    
    # The 'volume_SA' function is defined within the 'volume_force_sa' function environment
    def volume_SA(geo,theta):
        
        # Polygon definition
        N_div = 50
        alpha_1 = np.arcsin(geo.r_b/geo.D_wall*2)
        x = geo.r_o * cos(3*pi/2 - theta) - geo.r_b
        alpha_2 = np.arcsin(x/geo.D_wall*2)
        
        # Definition of the angles delimiting suction chambers and admission chamber
        phi_ssa_fix = get_phi_ssa(geo,theta,'f')
        phi_ssa_orb = get_phi_ssa(geo,theta,'o')
        
        # Volute angle related to the inner border of the 'sa' region. The 'sa' region is splitted into 2 regions,
        # the first one is comprised between the outer border of the fixed scroll and the casing wall, the second
        # is comprised between the outer border of the orbiting scroll and the casing wall.
        phi_polySA1_inner = np.linspace(geo.phi_ie_fix, phi_ssa_fix,N_div)
        angle_out_1 =  (geo.phi_ie_orb - alpha_2 + pi/2)%(2*pi) 
        angle_out_2 = (geo.phi_ie_fix +  pi + alpha_1 + pi/2)%(2*pi)
        if angle_out_1 > angle_out_2:
            angle_out_2 += 2*pi
          
        # Volute angle related to the outer border of the 'sa' region
        phi_polySA1_outer = np.linspace(angle_out_1,angle_out_2, N_div) 
        
        # Coordinates of the outer border of the 'sa' region are defined assumin circular casing for the machine
        x_polySA1_outer =   geo.D_wall/2*np.cos(phi_polySA1_outer) 
        y_polySA1_outer =   geo.D_wall/2*np.sin(phi_polySA1_outer) 
        
        # Coordinates of the inner border of the 'sa' regione are computed using the outer involute of the scroll
        (x_polySA1_inner,y_polySA1_inner) = geo.getCoord(phi_polySA1_inner,"o")
        
        # The polygon is defined as a concatenation of all the computed points
        x_polySA1 = np.concatenate((x_polySA1_inner, x_polySA1_outer, x_polySA1_inner[0]), axis=None)
        y_polySA1 = np.concatenate((y_polySA1_inner, y_polySA1_outer, y_polySA1_inner[0]), axis=None)
        xy_polySA1 = np.array([x_polySA1,y_polySA1])
        xy_polySA1 = xy_polySA1.transpose()
  
        # The same procedure is repeated for the other half of the 'sa' region
        phi_polySA2_inner = np.linspace(geo.phi_ie_orb, phi_ssa_orb, N_div)
        phi_polySA2_outer = np.linspace(angle_out_2, angle_out_1+2*pi, N_div)
        
        x_polySA2_outer =   geo.D_wall/2*np.cos(phi_polySA2_outer) 
        y_polySA2_outer =   geo.D_wall/2*np.sin(phi_polySA2_outer) 
        
        (x_polySA2_inner,y_polySA2_inner) = geo.getCoordOrbiting(phi_polySA2_inner,"o",theta)
        
        x_polySA2 = np.concatenate((x_polySA2_inner, x_polySA2_outer, x_polySA2_inner[0]), axis=None)
        y_polySA2 = np.concatenate((y_polySA2_inner, y_polySA2_outer, y_polySA2_inner[0]), axis=None)
        xy_polySA2 = np.array([x_polySA2,y_polySA2])
        xy_polySA2 = xy_polySA2.transpose()

        # The volume is computed as the area of the polygons defined, multiplied by the height of the scroll.
        V_sa1 = polygonArea(x_polySA1,y_polySA1)*geo.h_s
        V_sa2 = polygonArea(x_polySA2,y_polySA2)*geo.h_s    
        V_sa = V_sa1 + V_sa2
        
        # Centroid of the 'sa' region is computed
        cx1_sa, cy1_sa = polygonCentroid(x_polySA1,y_polySA1)
        cx2_sa, cy2_sa = polygonCentroid(x_polySA2,y_polySA2)
        cx_sa = (cx1_sa*V_sa1 + cx2_sa*V_sa2)/V_sa
        cy_sa = (cy1_sa*V_sa1 + cy2_sa*V_sa2)/V_sa
    
        return V_sa, xy_polySA1, xy_polySA2, cx_sa, cy_sa
    
    # Volume of the 'sa' region is computed by means of the dedicated function 'volume_SA'
    (V_sa, xy_polySA1, xy_polySA2, cx_sa, cy_sa) = volume_SA(geo,theta)
    
    # Derivative is computed by means of the definition of infinitesimal different quotient.
    # The infinitesimal increment is chosen as 2pi/360
    h = 2*pi/360
    (V_sa_ph, _, _, _, _) = volume_SA(geo,theta+h)
    (V_sa_mh, _, _, _, _) = volume_SA(geo,theta-h)
    
    dVdTheta_sa = (V_sa_ph - V_sa_mh)/2/h
    
    if force == True:
        # Computation of the forces acting in chamber 'sa' is performed 
        # No axial forces are considered neither centroid of the chamber.
        phi_ssa_orb = get_phi_ssa(geo,theta,'o')    
        phi_force_sa = np.linspace(phi_ssa_orb, geo.phi_ie_orb, N_div)   
        (nx_sa,ny_sa) = geo.getCoordNorm(phi_force_sa,'oo')
        (x_force_sa, y_force_sa) = geo.getCoordOrbiting(phi_force_sa,"o",theta)
        L = len(x_force_sa)
        dA_sa = geo.h_s*np.sqrt(np.power(x_force_sa[1:L]-x_force_sa[0:L-1],2)+np.power(y_force_sa[1:L]-y_force_sa[0:L-1],2))
        dfx_p_sa = dA_sa*(nx_sa[1:L]+nx_sa[0:L-1])/2
        dfy_p_sa = dA_sa*(ny_sa[1:L]+ny_sa[0:L-1])/2
        fx_p_sa = np.sum(dfx_p_sa) #mm^2 multipied by 1 Pa is multiplying by 1e-6
        fy_p_sa = np.sum(dfy_p_sa)
        rOx_sa = (x_force_sa[1:L] + x_force_sa[0:L-1])/2 - geo.r_o*cos(geo.phi_ie_fix - pi/2 - theta)
        rOy_sa = (y_force_sa[1:L] + y_force_sa[0:L-1])/2 - geo.r_o*sin(geo.phi_ie_fix - pi/2 - theta)
        
        # Moment exerted on the chamber
        MO_p_sa = 1e-6*np.sum(rOx_sa * dfy_p_sa - rOy_sa * dfx_p_sa) #mm^3 multiply by 1Pa is multiplying by 1
        
        # Directions of the forces acting on the chamber and their application sites are stored in the 'f_vect_sa' array
        f_vect_sa = np.vstack((x_force_sa, y_force_sa, nx_sa, ny_sa ))
    
        # Dictionary containing al the informations about the chamber 'sa'
        dict_sa = {"V" : V_sa, "dVdTheta" : dVdTheta_sa,"xy_polySA1": xy_polySA1 ,"xy_polySA2": xy_polySA2, "cx": cx_sa, "cy": cy_sa, "MO_p" : MO_p_sa, "f_vect" : f_vect_sa, "fx_p" : fx_p_sa, "fy_p" : fy_p_sa}
    else:
        dict_sa = {"V" : V_sa, "dVdTheta" : dVdTheta_sa,"xy_polySA1": xy_polySA1 ,"xy_polySA2": xy_polySA2, "cx": cx_sa, "cy": cy_sa, "MO_p" : [], "f_vect" : [], "fx_p" : [], "fy_p" : []}
        
    
    return dict_sa

# The function 'volume_force_c1' computes the volume, volume derivative, coordinates of the centroid and orbiting
# moment generated by the chamber 'c1'. It returns a dictionary containing these values.
def volume_force_c1(geo,theta,force):
    """
    Parameters
    ----------
    geo : class geo containing the geometric parameters
    theta : Orbiting angle
    alpha : index of the compression chamber

    Returns
    -------
    A dictionnary containing the volume, derivative of the volume, coordinate of the polynom, coordinate of the centroid and orbiting moment 
    generated by the chamber c1
    """ 
   
    # The number of compression chamber is computed.
    N_c = N_c_calculation(geo,theta)
    dict_c1= {}
    alpha = 1
    
    # The computation continues as long as compression chambers exist
    while alpha <= geo.N_c_max:
        # This if statement is used to make sure that computation is not performed on chambers that do not exist.
        if alpha <= N_c:
            
            # Definition of the polygons that approximates the compression chambers
            (x_polyC1, y_polyC1, N_div) = polygonDefinition(geo, geo.phi_ie_fix -2*pi*alpha - theta, geo.phi_ie_fix - 2*pi*(alpha-1) - theta, geo.phi_ie_fix - 2*pi*(alpha-1) - pi - theta, geo.phi_ie_fix -2*pi*(alpha)-  theta - pi, theta, "1")
            xy_polyC1 = np.array([x_polyC1,y_polyC1])
            xy_polyC1 = xy_polyC1.transpose() 
            
            # Volume is computd passing the vertices of the polygon to  the 'polygonArea' function and multuplying that area
            # times the height of the scroll.
            V_c1 = polygonArea(x_polyC1, y_polyC1)*geo.h_s
            
            # Derivative is computed by means of the definition of infinitesimal different quotient.
            # The infinitesimal increment is chosen as 2pi/360
            h = 2*pi/360
            (x_polyC1_ph, y_polyC1_ph, _) = polygonDefinition(geo, geo.phi_ie_fix -2*pi*alpha - theta - h, geo.phi_ie_fix - 2*pi*(alpha-1) - theta - h, geo.phi_ie_fix - 2*pi*(alpha-1) - pi - theta - h, geo.phi_ie_fix -2*pi*(alpha)-  theta - h - pi, theta + h, "1")
            (x_polyC1_mh, y_polyC1_mh, _) = polygonDefinition(geo, geo.phi_ie_fix -2*pi*alpha - theta + h, geo.phi_ie_fix - 2*pi*(alpha-1) - theta + h, geo.phi_ie_fix - 2*pi*(alpha-1) - pi - theta + h, geo.phi_ie_fix -2*pi*(alpha)-  theta + h - pi, theta - h, "1")
            V_c1_ph = polygonArea(x_polyC1_ph, y_polyC1_ph)*geo.h_s
            V_c1_mh = polygonArea(x_polyC1_mh, y_polyC1_mh)*geo.h_s        
            dVdTheta_c1 = (V_c1_ph - V_c1_mh)/2/h
            
            poly_c1 = Polygon(xy_polyC1)
            # Centroid is computed by means of the 'polygonCentroid' function if 'theta' > 0.001.
            c_c1 = poly_c1.centroid
            cx_c1 = c_c1.x
            cy_c1 = c_c1.y
        
            if force == True:
                # Computation of the forces acting in chamber 'c1' is performed
                phi_force_c1 = np.linspace(geo.phi_ie_fix - 2*pi*(alpha-1) - pi - theta, geo.phi_ie_fix -2*pi*(alpha)-  theta - pi, N_div)   
                (nx_c1, ny_c1) = geo.getCoordNorm(phi_force_c1,'oo')
                L = len(x_polyC1)
                pos_fx_c1 = x_polyC1[N_div:L]
                pos_fy_c1 = y_polyC1[N_div:L]
                L_pos = len(pos_fx_c1)
                dA_c1 = geo.h_s*np.sqrt(np.power(pos_fx_c1[1:L_pos] - pos_fx_c1[0:L_pos-1],2) + np.power(pos_fy_c1[1:L_pos] - pos_fy_c1[0:L_pos-1],2))
                dfx_p_c1 = dA_c1*(nx_c1[1:N_div]+nx_c1[0:N_div-1])/2
                dfy_p_c1 = dA_c1*(ny_c1[1:N_div]+ny_c1[0:N_div-1])/2
                fx_p_c1 = np.sum(dfx_p_c1) #mm^2 multipied by 1 Pa is multiplying by 1e-6
                fy_p_c1 = np.sum(dfy_p_c1)
                rOx_c1 = (x_polyC1[N_div:L-1] + x_polyC1[N_div+1:L])/2 - geo.r_o*cos(geo.phi_ie_fix - pi/2 - theta)
                rOy_c1 = (y_polyC1[N_div:L-1] + y_polyC1[N_div+1:L])/2 - geo.r_o*sin(geo.phi_ie_fix - pi/2 - theta)
                f_vect_c1 = np.vstack((pos_fx_c1, pos_fy_c1, nx_c1, ny_c1 ))
                        
                # Moment exerted on the chamber 
                MO_p_c1 = 1e-6*np.sum(rOx_c1 * dfy_p_c1 - rOy_c1 * dfx_p_c1) #mm^3 multiply by 1Pa is multiplying by 1
                
                # Dictionary containing al the informations about the chamber 'c1' related to the alpha-th chamber
                dict_c1[alpha] = {"V" : V_c1, "dVdTheta": dVdTheta_c1, "xy_poly": xy_polyC1, "cx" : cx_c1, "cy" : cy_c1, "MO_p" : MO_p_c1, "f_vect" : f_vect_c1, "fx_p" : fx_p_c1, "fy_p" : fy_p_c1}
            else:
                # Dictionary containing al the informations about the chamber 'c1' related to the alpha-th chamber
                dict_c1[alpha] = {"V" : V_c1, "dVdTheta": dVdTheta_c1, "xy_poly": xy_polyC1, "cx" : cx_c1, "cy" : cy_c1, "MO_p" : [], "f_vect" : [], "fx_p" : [], "fy_p" : []}
                
        else:
           # If the chamber does not exist all the values are set to 'nan'
           V_c1 = dVdTheta_c1 = cx_c1 = cy_c1 = MO_p_c1 = pos_fx_c1 = pos_fy_c1 = nx_c1 = ny_c1 = np.nan  
           f_vect_c1 = np.vstack((pos_fx_c1, pos_fy_c1, nx_c1, ny_c1 ))
           xy_polyC1 = np.array([[np.nan, np.nan],[np.nan, np.nan]])
           xy_polyC1 = xy_polyC1.transpose() 
           
           # Dictionary containing al the informations about the chamber 'c1' related to the alpha-th chamber
           dict_c1[alpha] = {"V" : V_c1, "dVdTheta": dVdTheta_c1, "xy_poly": xy_polyC1, "cx" : cx_c1, "cy" : cy_c1, "MO_p" : MO_p_c1, "f_vect" : f_vect_c1}
      
        alpha += 1 

    return dict_c1     

# The same procedure seen for 'c1' chamber is repeated for the compression chamber 'c2'.
def volume_force_c2(geo,theta,force):
    """
    Parameters
    ----------
    geo : class geo containing the geometric parameters
    theta : Orbiting angle
    alpha : index of the compression chamber

    Returns
    -------
    A dictionnary containing the volume, derivative of the volume, coordinate of the polynom, coordinate of the centroid and orbiting moment 
    generated by the chamber c2

    """ 
   
    N_c = N_c_calculation(geo,theta)
    dict_c2= {}
    alpha = 1
    
    while alpha <= geo.N_c_max:
        # This if statement is used to make sure that computation is not performed on chambers that do not exist
        if alpha<= N_c:
             # Definition of the polygons that approximates the compression chambers
             (x_polyC2, y_polyC2, N_div) = polygonDefinition(geo, geo.phi_ie_orb -2*pi*alpha - theta - pi, geo.phi_ie_orb - 2*pi*(alpha-1) - theta - pi, geo.phi_ie_orb - 2*pi*(alpha-1)  - theta, geo.phi_ie_orb -2*pi*(alpha) - theta , theta, "2")
             xy_polyC2 = np.array([x_polyC2,y_polyC2])
             xy_polyC2 = xy_polyC2.transpose() 
             
             # Results of the 'c1' chamber are used for the 'c2' chamber as well.
             dict_c1 = volume_force_c1(geo,theta,False)
        
             # Volume
             V_c2 = dict_c1[alpha]["V"]
            
             # Derivative
             dVdTheta_c2 =  dict_c1[alpha]["dVdTheta"]
            
             # Centroid of e 'c1' chamber
             cx_c1 = dict_c1[alpha]["cx"]
             cy_c1 = dict_c1[alpha]["cy"]
            
             # Centroid of e 'c2' chamber
             (cx_c2,cy_c2) = (-cx_c1+geo.r_o*cos(geo.phi_ie_fix-pi/2-theta), -cy_c1+geo.r_o*sin(geo.phi_ie_fix-pi/2-theta))
             
             # Orbiting moment and forces
             phi_force_c2 = np.linspace( geo.phi_ie_orb - 2*pi*(alpha-1) - theta , geo.phi_ie_orb -2*pi*alpha - theta ,  N_div)   
             
             if force == True:
                 # Computation of the forces acting in chamber 'c2' is performed
                 (nx_c2, ny_c2) = geo.getCoordNorm(phi_force_c2,'oi')
                 L = len(x_polyC2)
                 pos_fx_c2 = x_polyC2[N_div:L]
                 pos_fy_c2 = y_polyC2[N_div:L]
                 L_pos = len(pos_fx_c2)
                 dA_c2 = geo.h_s*np.sqrt(np.power(pos_fx_c2[1:L_pos] - pos_fx_c2[0:L_pos-1],2) + np.power(pos_fy_c2[1:L_pos] - pos_fy_c2[0:L_pos-1],2))
                 dfx_p_c2 = dA_c2*(nx_c2[1:N_div]+nx_c2[0:N_div-1])/2
                 dfy_p_c2 = dA_c2*(ny_c2[1:N_div]+ny_c2[0:N_div-1])/2
                 fx_p_c2 = np.sum(dfx_p_c2) #mm^2 multipied by 1 Pa is multiplying by 1e-6
                 fy_p_c2 = np.sum(dfy_p_c2)
                 rOx_c2 = (x_polyC2[N_div:L-1] + x_polyC2[N_div+1:L])/2 - geo.r_o*cos(geo.phi_ie_fix - pi/2 - theta)
                 rOy_c2 = (y_polyC2[N_div:L-1] + y_polyC2[N_div+1:L])/2 - geo.r_o*sin(geo.phi_ie_fix - pi/2 - theta)
                 f_vect_c2 = np.vstack((pos_fx_c2, pos_fy_c2, nx_c2, ny_c2 ))   
                 MO_p_c2 = 1e-6*np.sum(rOx_c2 * dfy_p_c2 - rOy_c2 * dfx_p_c2) #mm^3 multiply by 1Pa is multiplying by 1
             
                 # Dictionary containing al the informations about the chamber 'c2' related to the alpha-th chamber
                 dict_c2[alpha] = {"V" : V_c2, "dVdTheta": dVdTheta_c2, "xy_poly": xy_polyC2, "cx" : cx_c2, "cy" : cy_c2, "MO_p" : MO_p_c2, "f_vect" : f_vect_c2, "fx_p" : fx_p_c2, "fy_p" : fy_p_c2}
             else:
                 dict_c2[alpha] = {"V" : V_c2, "dVdTheta": dVdTheta_c2, "xy_poly": xy_polyC2, "cx" : cx_c2, "cy" : cy_c2, "MO_p" : [], "f_vect" :[], "fx_p" : [], "fy_p" : []}
                
                
        else: 
             # Chamber does not exist
             V_c2 = dVdTheta_c2 = cx_c2 = cy_c2 = MO_p_c2 = pos_fx_c2 = pos_fy_c2 = nx_c2 = ny_c2 = np.nan 
             f_vect_c2 = np.vstack((pos_fx_c2, pos_fy_c2, nx_c2, ny_c2 ))   
             xy_polyC2 = np.array([[np.nan, np.nan],[np.nan, np.nan]])
             xy_polyC2 = xy_polyC2.transpose() 
             
             # Dictionary containing al the informations about the chamber 'c1' related to the alpha-th chamber
             dict_c2[alpha] = {"V" : V_c2, "dVdTheta": dVdTheta_c2, "xy_poly": xy_polyC2, "cx" : cx_c2, "cy" : cy_c2, "MO_p" : MO_p_c2, "f_vect" : f_vect_c2}
        
        alpha += 1 
    
    return dict_c2  

# The function 'volume_force_d1' computes the volume, volume derivative, coordinates of the centroid and orbiting
# moment generated by the chamber 'd1'. It returns a dictionary containing these values.
def volume_force_d1(geo,theta,force):
    """
    Parameters
    ----------
    geo : class geo containing the geometric parameters
    theta : Orbiting angle

    Returns
    -------
    A dictionnary containing the volume, derivative of the volume, coordinate of the polynom, coordinate of the centroid and orbiting moment 
    generated by the chamber d1
    """    
    # Computation of the number of discharge chamber     
    N_c = N_c_calculation(geo,theta)

    # Definition of the polygons that approximates the discharge chambers
    (x_polyD1, y_polyD1, N_div) = polygonDefinition(geo, geo.phi_ie_fix - 2*pi*N_c - theta, geo.phi_is, geo.phi_os, geo.phi_ie_fix -2*pi*N_c - theta - pi , theta, "1")
    xy_polyD1 = np.array([x_polyD1,y_polyD1])
    xy_polyD1 = xy_polyD1.transpose()

    # Volume is computd passing the vertices of the polygon to  the 'polygonArea' function and multuplying that area
    # times the height of the scroll.
    V_d1 = polygonArea(x_polyD1, y_polyD1)*geo.h_s

    # Derivative is computed by means of the definition of infinitesimal different quotient.
    # The infinitesimal increment is chosen as 2pi/360
    h = 2*pi/360
    (x_polyD1_ph, y_polyD1_ph, _) = polygonDefinition(geo, geo.phi_ie_fix - 2*pi*N_c - theta - h, geo.phi_is, geo.phi_os, geo.phi_ie_fix -2*pi*N_c - theta - h - pi , theta + h, "1")
    (x_polyD1_mh, y_polyD1_mh, _) = polygonDefinition(geo, geo.phi_ie_fix - 2*pi*N_c - theta + h, geo.phi_is, geo.phi_os, geo.phi_ie_fix -2*pi*N_c - theta + h - pi , theta - h, "1")
    V_d1_ph = polygonArea(x_polyD1_ph, y_polyD1_ph)*geo.h_s
    V_d1_mh = polygonArea(x_polyD1_mh, y_polyD1_mh)*geo.h_s        
    dVdTheta_d1 = (V_d1_ph - V_d1_mh)/2/h
    
    poly_d1 = Polygon(xy_polyD1)
    # Centroid is computed by means of the 'polygonCentroid' function if 'theta' > 0.001.
    c_d1 = poly_d1.centroid
    cx_d1 = c_d1.x
    cy_d1 = c_d1.y
    
    if force == True:
        # Computation of the forces acting in chamber 'd1' is performed
        phi_force_d1 = np.linspace(geo.phi_os, geo.phi_ie_fix -2*pi*N_c - theta - pi, N_div)   
        (nx_d1, ny_d1) = geo.getCoordNorm(phi_force_d1,'oo')
        L = len(x_polyD1)
        pos_fx_d1 = x_polyD1[N_div:L]
        pos_fy_d1 = y_polyD1[N_div:L]
        L_pos = len(pos_fx_d1)
        dA_d1 = geo.h_s*np.sqrt(np.power(pos_fx_d1[1:L_pos] - pos_fx_d1[0:L_pos-1],2) + np.power(pos_fy_d1[1:L_pos] - pos_fy_d1[0:L_pos-1],2))
        dfx_p_d1 = dA_d1*(nx_d1[1:N_div] + nx_d1[0:N_div-1])/2
        dfy_p_d1 = dA_d1*(ny_d1[1:N_div] + ny_d1[0:N_div-1])/2
        fx_p_d1 = np.sum(dfx_p_d1) #mm^2 multipied by 1 Pa is multiplying by 1e-6
        fy_p_d1 = np.sum(dfy_p_d1)
        rOx_d1 = (x_polyD1[N_div:L-1] + x_polyD1[N_div+1:L])/2 - geo.r_o*cos(geo.phi_ie_fix - pi/2 - theta)
        rOy_d1 = (y_polyD1[N_div:L-1] + y_polyD1[N_div+1:L])/2 - geo.r_o*sin(geo.phi_ie_fix - pi/2 - theta)
        
        # Moment exerted on the chamber
        MO_p_d1 = 1e-6*np.sum(rOx_d1 * dfy_p_d1 - rOy_d1 * dfx_p_d1) #mm^3 multiply by 1Pa is multiplying by 1
        
        f_vect_d1 = np.vstack((pos_fx_d1, pos_fy_d1, nx_d1, ny_d1 )) 
        
        # Dictionary containing al the informations about the chamber 'd1'. 
        dict_d1 = {"V" : V_d1, "dVdTheta": dVdTheta_d1, "xy_poly": xy_polyD1, "cx" : cx_d1, "cy" : cy_d1, "MO_p" : MO_p_d1, "f_vect" : f_vect_d1, "fx_p" : fx_p_d1, "fy_p" : fy_p_d1}
    else:
        dict_d1 = {"V" : V_d1, "dVdTheta": dVdTheta_d1, "xy_poly": xy_polyD1, "cx" : cx_d1, "cy" : cy_d1, "MO_p" : [], "f_vect" : [], "fx_p" : [], "fy_p" : []}        
    
    return dict_d1

# The computation done for chamber 'd1' is repeated for chamber 'd2'
def volume_force_d2(geo,theta,force):
    """
    Parameters
    ----------
    geo : class geo containing the geometric parameters
    theta : Orbiting angle

    Returns
    -------
    A dictionnary containing the volume, derivative of the volume, coordinate of the polynom, coordinate of the centroid and orbiting moment 
    generated by the chamber d2
    """    
    
    N_c = N_c_calculation(geo,theta)

    dict_d1 = volume_force_d1(geo,theta,False)
    
    # Volume
    V_d2 = dict_d1["V"]
        
    # Derivative
    dVdTheta_d2 =  dict_d1["dVdTheta"]
        
    # Centroid
    cx_d1 = dict_d1["cx"]
    cy_d1 = dict_d1["cy"]
            
    (cx_d2,cy_d2) = (-cx_d1+geo.r_o*cos(geo.phi_ie_fix-pi/2-theta), -cy_d1+geo.r_o*sin(geo.phi_ie_fix-pi/2-theta))
    
    # Polygon definition
    (x_polyD2, y_polyD2, N_div) = polygonDefinition(geo, geo.phi_ie_orb - 2*pi*N_c - theta - pi, geo.phi_os, geo.phi_is, geo.phi_ie_orb -2*pi*N_c - theta, theta, "2")
    xy_polyD2 = np.array([x_polyD2, y_polyD2])
    xy_polyD2 = xy_polyD2.transpose()

    if force == True:
        # Orbiting moment and forces
        phi_force_d2 = np.linspace(geo.phi_is, geo.phi_ie_fix -2*pi*N_c - theta, N_div)   
        (nx_d2, ny_d2) = geo.getCoordNorm(phi_force_d2,'oi')
        L = len(x_polyD2)
        pos_fx_d2 = x_polyD2[N_div:L]
        pos_fy_d2 = y_polyD2[N_div:L]
        L_pos = len(pos_fx_d2)
        dA_d2 = geo.h_s*np.sqrt(np.power(pos_fx_d2[1:L_pos] - pos_fx_d2[0:L_pos-1],2) + np.power(pos_fy_d2[1:L_pos] - pos_fy_d2[0:L_pos-1],2))
        dfx_p_d2 = dA_d2*(nx_d2[1:N_div]+nx_d2[0:N_div-1])/2
        dfy_p_d2 = dA_d2*(ny_d2[1:N_div]+ny_d2[0:N_div-1])/2
        fx_p_d2 = np.sum(dfx_p_d2) #mm^2 multipied by 1 Pa is multiplying by 1e-6
        fy_p_d2 = np.sum(dfy_p_d2)
        rOx_d2 = (x_polyD2[N_div:L-1] + x_polyD2[N_div+1:L])/2 - geo.r_o*cos(geo.phi_ie_fix - pi/2 - theta)
        rOy_d2 = (y_polyD2[N_div:L-1] + y_polyD2[N_div+1:L])/2 - geo.r_o*sin(geo.phi_ie_fix - pi/2 - theta)
        
        MO_p_d2 = 1e-6*np.sum(rOx_d2 * dfy_p_d2 - rOy_d2* dfx_p_d2) #mm^3 multiply by 1Pa is multiplying by 1
        
        f_vect_d2 = np.vstack((pos_fx_d2, pos_fy_d2, nx_d2, ny_d2 )) 
         
        dict_d2 = {"V" : V_d2, "dVdTheta": dVdTheta_d2, "xy_poly": xy_polyD2, "cx" : cx_d2, "cy" : cy_d2, "MO_p" : MO_p_d2, "f_vect" : f_vect_d2, "fx_p" : fx_p_d2, "fy_p" : fy_p_d2}    
    else:
        dict_d2 = {"V" : V_d2, "dVdTheta": dVdTheta_d2, "xy_poly": xy_polyD2, "cx" : cx_d2, "cy" : cy_d2, "MO_p" : [], "f_vect" : [], "fx_p" : [], "fy_p" : []}    
        
    
    return dict_d2 
  

def volume_force_dd(geo,theta,force):
    """
    Parameters
    ----------
    geo : class geo containing the geometric parameters
    theta : Orbiting angle

    Returns
    -------
    A dictionnary containing the volume, derivative of the volume, coordinate of the polynom, coordinate of the centroid and orbiting moment 
    generated by the chamber dd
    """    

    N_div = 50
    
    # The 'volume_dd' function is defined within the 'volume_force_dd' function environment
    def volume_dd(geo, theta):
        
        # Definition of the polygons that approximates the discharge chambers 
        phi_start = geo.phi_is
        (x_start, y_start) = geo.getCoord(phi_start, 'i')
        
        # First arc 
        phi_arc1 = np.linspace(geo.t_a12, geo.t_a11, N_div)
        x_arc1 = geo.x_a1 + geo.r_a1*np.cos(phi_arc1)
        y_arc1 = geo.y_a1 + geo.r_a1*np.sin(phi_arc1)

        # Second arc    
        phi_arc2 = np.linspace(geo.t_a21,geo.t_a22, N_div)
        x_arc2 = geo.x_a2 + geo.r_a2*np.cos(phi_arc2) 
        y_arc2 = geo.y_a2 + geo.r_a2*np.sin(phi_arc2) 
        
        # The three different discharge geometries are considered for the computation of the volume.
        # 'x_polyDDf' and 'y_polyDDf' are arrays containing the coordinates of the polygon that approximates
        # the discharge region of the fixed scroll
        if geo.dischargeType == 'aala':
             phi_curveArc0 = np.linspace(geo.t_a02, geo.t_a01, N_div)
             x_arc0 = geo.x_a0 + geo.r_a0*np.cos(phi_curveArc0)
             y_arc0 = geo.y_a0 + geo.r_a0*np.sin(phi_curveArc0)

             x_L = np.linspace(geo.x_a1t, geo.x_a2t, N_div)
             y_L = geo.m_l*x_L + geo.b_l   
             
             x_polyDDf = np.concatenate((x_arc0, x_arc1,  x_L, x_arc2), axis = None)
             y_polyDDf = np.concatenate((y_arc0, y_arc1,  y_L, y_arc2), axis = None) 
             
             k = 4
             
        elif geo.dischargeType == 'ala':
            
             x_L = np.linspace(geo.x_a1t, geo.x_a2t, N_div)
             y_L = geo.m_l*x_L + geo.b_l  
             
             x_polyDDf = np.concatenate((x_arc1,  x_L, x_arc2), axis = None)
             y_polyDDf = np.concatenate((y_arc1,  y_L, y_arc2), axis = None)                   

             k = 3     
             
        elif geo.dischargeType == '2a':     
             x_polyDDf = np.concatenate((x_arc1, x_arc2), axis = None)
             y_polyDDf = np.concatenate((y_arc1, y_arc2), axis = None)  
             
             k = 2             
        
        # 'x_polyDDo' and 'y_polyDDo' are arrays containing the coordinates of the polygon that approximates
        # the discharge region of the orbiting scroll 
        x_polyDDo = - x_polyDDf + geo.r_o*cos(geo.phi_ie_fix - pi/2 - theta)
        y_polyDDo = - y_polyDDf + geo.r_o*sin(geo.phi_ie_fix - pi/2 - theta)
        
        # Last link 
        x_end = x_polyDDf[0]
        y_end = y_polyDDf[0]
        
        # 'x_polyDD' and 'y_polyDD' are arrays containing the coordinates of the polygon that approximates
        # the whole discharge region.
        x_polyDD = np.concatenate((x_polyDDf, x_polyDDo, x_end), axis=None)
        y_polyDD = np.concatenate((y_polyDDf, y_polyDDo, y_end), axis=None)
           
        xy_polyDD = np.array([x_polyDD,y_polyDD])
        xy_polyDD = xy_polyDD.transpose()
        
        # Volume is computd passing the vertices of the polygon to  the 'polygonArea' function and multuplying that area
        # times the height of the scroll.
        V_dd = polygonArea(x_polyDD, y_polyDD)*geo.h_s
        
        return V_dd, xy_polyDD, x_polyDD, y_polyDD
    
    # Computation of the volume is called.
    (V_dd, xy_polyDD, x_polyDD, y_polyDD) = volume_dd(geo,theta)  
    
    # Derivative is computed by means of the definition of infinitesimal different quotient.
    # The infinitesimal increment is chosen as 2pi/360
    h = 2*pi/360
    (V_dd_ph, _, _, _) = volume_dd(geo,theta+h)
    (V_dd_mh, _, _, _) = volume_dd(geo,theta-h)
    dVdTheta_dd = (V_dd_ph - V_dd_mh)/2/h 
   
    poly_dd = Polygon(xy_polyDD)
    # Centroid is computed by means of the 'polygonCentroid' function if 'theta' > 0.001.
    c_dd = poly_dd.centroid
    cx_dd = c_dd.x
    cy_dd = c_dd.y
   
    if force == True:
        # Computation of the forces acting in chamber 'dd' is performed
        phi_force_arc1 = np.linspace((geo.t_a12 - pi) , (geo.t_a11 - pi), N_div)
        phi_force_arc2 = np.linspace((geo.t_a21 - pi) , (geo.t_a22 - pi), N_div)
        (nx_arc1, ny_arc1) = (np.cos(phi_force_arc1), np.sin(phi_force_arc1))
        (nx_arc2, ny_arc2) = (- np.cos(phi_force_arc2), - np.sin(phi_force_arc2))
        # The computation of forces is performed based on the discharge geometry chosen.     
        if geo.dischargeType == 'aala':
             phi_force_arc0 = np.linspace((geo.t_a02 - pi) , (geo.t_a01 - pi), N_div)
             (nx_arc0, ny_arc0) = (np.cos(phi_force_arc0), np.sin(phi_force_arc0))
             (nx_L, ny_L) = (np.full((1, N_div),  geo.m_l/sqrt((1 + geo.m_l**2))  ), np.full((1, N_div), -1/sqrt((1 + geo.m_l**2))))
             nx_dd = np.concatenate((nx_arc0,  nx_arc1, nx_L, nx_arc2), axis=None)   
             ny_dd = np.concatenate((ny_arc0,  ny_arc1, ny_L, ny_arc2), axis=None) 
             k = 4
        elif geo.dischargeType == 'ala':
             (nx_L, ny_L) = (np.full((1, N_div),  geo.m_l/sqrt((1 + geo.m_l**2))  ), np.full((1, N_div), -1/sqrt((1 + geo.m_l**2))))
             nx_dd = np.concatenate((nx_arc1, nx_L, nx_arc2), axis=None)   
             ny_dd = np.concatenate((ny_arc1, ny_L, ny_arc2), axis=None) 
             k = 3     
        elif geo.dischargeType == '2a':     
             nx_dd = np.concatenate((nx_arc1,  nx_arc2), axis=None)   
             ny_dd = np.concatenate((ny_arc1,  ny_arc2), axis=None) 
             k = 2             
        L = len(x_polyDD)
        pos_fx_dd = x_polyDD[k*N_div:L-1]
        pos_fy_dd = y_polyDD[k*N_div:L-1]
        L_pos = len(pos_fx_dd)
        dA_dd = geo.h_s*np.sqrt(np.power(pos_fx_dd[1:L_pos] - pos_fx_dd[0:L_pos-1],2) + np.power(pos_fy_dd[1:L_pos] - pos_fy_dd[0:L_pos-1],2))
        dfx_p_dd = dA_dd*(nx_dd[1:k*N_div]+nx_dd[0:k*N_div-1])/2
        dfy_p_dd = dA_dd*(ny_dd[1:k*N_div]+ny_dd[0:k*N_div-1])/2
        fx_p_dd = np.sum(dfx_p_dd) #mm^2 multipied by 1 Pa is multiplying by 1e-6
        fy_p_dd = np.sum(dfy_p_dd)
        rOx_dd = (x_polyDD[k*N_div:L-2] + x_polyDD[k*N_div+1:L-1])/2 - geo.r_o*cos(geo.phi_ie_fix - pi/2 - theta)
        rOy_dd = (y_polyDD[k*N_div:L-2] + y_polyDD[k*N_div+1:L-1])/2 - geo.r_o*sin(geo.phi_ie_fix - pi/2 - theta)
       
        # Moment exerted on the chamber
        MO_p_dd = 1e-6*np.sum(rOx_dd * dfy_p_dd - rOy_dd * dfx_p_dd) #mm^3 multiply by 1Pa is multiplying by 1
        
        f_vect_dd = np.vstack((pos_fx_dd, pos_fy_dd, nx_dd, ny_dd )) 
        
        # Dictionary containing al the informations about the chamber 'd1'. 
        dict_dd = {"V" : V_dd, "dVdTheta": dVdTheta_dd, "xy_poly": xy_polyDD, "cx" : cx_dd, "cy" : cy_dd, "MO_p" : MO_p_dd, "f_vect" : f_vect_dd, "fx_p" : fx_p_dd, "fy_p" : fy_p_dd}
    else: 
        dict_dd = {"V" : V_dd, "dVdTheta": dVdTheta_dd, "xy_poly": xy_polyDD, "cx" : cx_dd, "cy" : cy_dd, "MO_p" : [], "f_vect" : [], "fx_p" : [], "fy_p" : []}

    return dict_dd  



def volumeRatio(geo, theta):
   
    # The volume displacement is the volume displaced by the machine in one turn of the machine.
    V_disp = -2*pi*geo.h_s*geo.r_b*geo.r_o*(3*pi-2*geo.phi_ie_fix+geo.phi_i0+geo.phi_o0)
    
    # VOlume do te discharge right after the discharge angle.
    V_d1_d = volume_force_d1(geo, geo.theta_d + 0.01)["V"]
    
    # The volume ratio is defined as the volume displaced by the machine divided by the volume of the discharge chamber
    # at the discharge angle.
    V_ratio = V_disp/2/V_d1_d
    
    # theta_uu = (theta_d - geo.t_a11 - geo.t_a12 )%(2*pi)
    # theta_u = (-geo.t_a11 + geo.phi_ie - pi/2)%(2*pi)

    return V_ratio, V_disp

# The 'suctionFlowArea' function computes suction area section and its coordinates.
def suctionFlowArea(geo,theta):
    
    # The 'getCoord' method computes the coordinates of the ending part of the scroll
    (x_e, y_e) = geo.getCoord( geo.phi_ie_fix, "i")
    # The 'get_phi_ssa' function computes the angle delimiting the suction chamber and the admission chamber depending
    # on the crank angle.
    phi_ssa_orb = get_phi_ssa(geo,theta,'o') 
    (x_ssa, y_ssa) = geo.getCoordOrbiting(phi_ssa_orb, "o", theta)
    
    # Coordinates of the suction section are stored within two arrays
    x_suction = np.concatenate((x_e, x_ssa), axis=None)
    y_suction = np.concatenate((y_e, y_ssa), axis=None)
    
    # Suction area is computed based on scroll heigth and coordinate of the ending part of the scroll
    A_suction = geo.h_s*sqrt((x_e - x_ssa)**2 + (y_e - y_ssa)**2)
    
    return x_suction, y_suction, A_suction
    
# This function computes the area that divides that 'dd' chamber from the 'd1' and 'd2' chambers before the merging
# process is completed
def area_d_dd(geo,theta):
    
    N_div = 100000
    
    # Fixed curve definition
    
    # First curve - 
    phi_curve_fi = np.linspace(geo.phi_ie_fix - theta, geo.phi_is, N_div)
    (x_curve_fi, y_curve_fi) = geo.getCoord( phi_curve_fi, 'i')    
    
    # 
    phi_arc1 = np.linspace(geo.t_a12, geo.t_a11, N_div)
    x_arc1 = geo.x_a1 + geo.r_a1*np.cos(phi_arc1)
    y_arc1 = geo.y_a1 + geo.r_a1*np.sin(phi_arc1)
    
    if geo.dischargeType == 'aala':
        phi_curveArc0 = np.linspace(geo.t_a02, geo.t_a01, N_div)
        x_arc0 = geo.x_a0 + geo.r_a0*np.cos(phi_curveArc0)
        y_arc0 = geo.y_a0 + geo.r_a0*np.sin(phi_curveArc0)

        x_L = np.linspace(geo.x_a1t, geo.x_a2t, N_div)
        y_L = geo.m_l*x_L + geo.b_l   
        
        # Fixed scroll curve
        x_fscroll = np.concatenate((x_curve_fi, x_arc0, x_arc1), axis=None)
        y_fscroll = np.concatenate((y_curve_fi, y_arc0, y_arc1), axis=None)  
    
    else:
        x_fscroll = np.concatenate((x_curve_fi, x_arc1), axis=None)
        y_fscroll = np.concatenate((y_curve_fi, y_arc1), axis=None)     
    
    # Definition of the point on the orbiting scroll
    alpha_contact_2 = (geo.phi_ie_fix - theta - pi/2 - pi)%(2*pi)

    x_arc2 = -(geo.x_a2 + geo.r_a2*cos(alpha_contact_2)) + geo.r_o*cos(geo.phi_ie_fix - pi/2 - theta)
    y_arc2 = -(geo.y_a2 + geo.r_a2*sin(alpha_contact_2)) + geo.r_o*sin(geo.phi_ie_fix - pi/2 - theta)
    
    
    # (x_arc2, y_arc2) = geo.getCoordOrbiting( geo.phi_os, "o", theta)
    

    # Minimization of the distance from the point on the fixed scroll to the point on the orbiting scroll
    # Definition of the function to minimise
    dist = np.sqrt((x_fscroll - x_arc2)**2 + (y_fscroll - y_arc2)**2)
    dist_min = min(dist)
    index_min = np.argmin(dist)
    
    # Definition of the points that delimit the section that divides the 'dd' chamber from the 'd1' and 'd2' chambers.
    x_d_dd = np.concatenate((x_fscroll[index_min], x_arc2), axis=None)
    y_d_dd = np.concatenate((y_fscroll[index_min], y_arc2), axis=None)
    
    # The section area is computed based on scroll height and 
    A_d_dd = geo.h_s*dist_min
    
    return x_d_dd, y_d_dd, A_d_dd


def area_d_dd_bis(geo,theta):
    
    N_div = 100
    
    N_c_max = np.floor((geo.phi_ie - geo.phi_os - pi)/2/pi)
    theta_d = geo.phi_ie - geo.phi_os -2*pi*N_c_max - pi
    theta_u = (-geo.t_a11 + geo.phi_ie - pi/2)%(2*pi)
    
    # Fixed curve definition
    
    # First curve
    phi_curve_fi = np.linspace(geo.phi_ie  - theta - pi, geo.phi_is, N_div)
    (x_curve_fi, y_curve_fi) = geo.getCoord( phi_curve_fi, 'i')    
    
    # Arc1  
    phi_arc1 = np.linspace(geo.t_a12, geo.t_a11, N_div)
    x_arc1 = geo.x_a1 + geo.r_a1*np.cos(phi_arc1)
    y_arc1 = geo.y_a1 + geo.r_a1*np.sin(phi_arc1)
    
    # Fixed scroll curve
    
    x_fix = np.concatenate((x_curve_fi, x_arc1), axis=None)
    y_fix = np.concatenate((y_curve_fi, y_arc1), axis=None)   
    geo.x_fix = x_fix
    geo.y_fix = y_fix
    
    # Orbiting scroll curve definition
    # Arc2 
    phi_arc2 = np.linspace(geo.t_a22, geo.t_a21, N_div)
    x_arc2 = -geo.x_a2 - geo.r_a2*np.cos(phi_arc2) + geo.r_o*cos(geo.phi_ie - pi/2 - theta)
    y_arc2 = -geo.y_a2 - geo.r_a2*np.sin(phi_arc2) + geo.r_o*sin(geo.phi_ie - pi/2 - theta)    
    
    # Orbiting curge
    phi_curve_oo = np.linspace( geo.phi_ie - theta - 2*pi,geo.phi_os, N_div)
    (x_curve_oo, y_curve_oo) = geo.getCoordOrbiting( phi_curve_oo, 'o', theta)   
    
    # Merging
    x_orb = np.concatenate((x_curve_oo, x_arc2), axis=None)
    y_orb = np.concatenate((y_curve_oo, y_arc2), axis=None)  
    geo.x_orb = x_orb
    geo.y_orb = y_orb
    
    
    dist = np.zeros((len(x_fix),len(x_orb)))
    for i in range(len(x_fix)):
        for j in range(len(x_orb)):
            dist[i,j] = np.sqrt((x_fix[i] - x_orb[j])**2 + (y_fix[i] - y_orb[j])**2)
    
    dist_min = np.min(dist)
    index_min = np.argmin(dist)
    index_min = np.unravel_index(dist.argmin(), dist.shape)
    geo.x_d_dd_bis = np.concatenate((x_fix[index_min[0]], x_orb[index_min[1]]), axis=None)
    geo.y_d_dd_bis = np.concatenate((y_fix[index_min[0]], y_orb[index_min[1]]), axis=None)
    
    A_d_dd = geo.h_s*dist_min
    
    return A_d_dd

       
def area_d_dd_bell(geo,theta):
    
    N_div = 10000
    
    # Fixed curve definition
    
    # First curve
    phi_curve_fi = np.linspace(geo.phi_ie  - theta, geo.phi_is, N_div)
    (x_curve_fi, y_curve_fi) = geo.getCoord( phi_curve_fi, 'i')    
    
      
    phi_arc1 = np.linspace(geo.t_a12, geo.t_a11, N_div)
    x_arc1 = geo.x_a1 + geo.r_a1*np.cos(phi_arc1)
    y_arc1 = geo.y_a1 + geo.r_a1*np.sin(phi_arc1)
    
    # Fixed scroll curve
    
    x_fscroll = np.concatenate((x_curve_fi, x_arc1), axis=None)
    y_fscroll = np.concatenate((y_curve_fi, y_arc1), axis=None)     
    
    # Point from the orbiting scroll
    
    omega = geo.phi_os - geo.phi_is + pi
    
    
    alpha_contact_2 = ((geo.phi_ie - theta - pi/2 - pi + omega) )
    if alpha_contact_2 > geo.t_a22 + 2*pi:
        alpha_contact_2 = geo.t_a22 + 2*pi
    
    (x_arc2, y_arc2) = geo.getCoordOrbiting(geo.phi_os, "o", theta)
    
    # Minimization of the distance 
    dist = np.sqrt((x_fscroll - x_arc2)**2 + (y_fscroll - y_arc2)**2)
    dist_min = min(dist)
    index_min = np.argmin(dist)
    
    geo.x_d_dd_bell = np.concatenate((x_fscroll[index_min], x_arc2), axis=None)
    geo.y_d_dd_bell = np.concatenate((y_fscroll[index_min], y_arc2), axis=None)
    
    A_d_dd = geo.h_s*dist_min
    
    return A_d_dd    
    
    
#     #to test  next functions
# import matplotlib.pyplot as plt
# from math import pi,sin,cos, sqrt, atan2, atan, acos, tan
# import numpy as np
# t1 = np.linspace(0, 2*pi, 100)
# x1 = 0 + 2*np.cos(t1)
# y1 = 0 + 2*np.sin(t1) 


# t2 = np.linspace(0, 2*pi, 100)  
# x2 = 2 + 1*np.cos(t2)
# y2 = 0 + 1*np.sin(t2)     


# def distance(x1,y1,x2,y2):   
#     dist = np.sqrt((x1 - y2)**2 + (y1 - y2)**2)
#     return dist  

# x1 = geo.x_odis
# y1 = geo.y_odis
# x2 = geo.xport
# y2 = geo.yport

# Cross product between vectors is implemented
def crossProduct(p,q):
    product = p[0]*q[1] - p[1]*q[0]
    return product


# This function verifies whether two segments intersect or not. It return '1' if they do, otherwise it returns '0'
def crossSegment(p1, p2, p3, p4):
    
    d1 = crossProduct((p1 - p3),(p4 - p3))
    d2 = crossProduct((p2 - p3),(p4 - p3))
    d3 = crossProduct((p3 - p1),(p2 - p1))
    d4 = crossProduct((p4 - p1),(p2 - p1))
    
    if ((d1<0 and d2>0) or (d1>0 and d2<0)) and ((d3 >0 and d4 <0) or (d3<0 and d4>0)):
        intersection = 1  
    else:
        intersection = 0
    return intersection   
         

def findIntersection(x1,y1,x2,y2):  
    
    """
    Function allowing to find the intersection points between two polyngons (x1,y1) and (x2, y2)
    The distance steps of the two polygons must be small enough 
    (we must not see the lines) for the algorithms to be accurate 
    
    Return a matrix containing n (number of intersections) lines with the index of the points 
    from polygon 1 and poloygon 2
    
    """
    L1 = len(x1)
    L2 = len(x2)
    dist = np.zeros((L1,L2))
    
    radius1 = np.sqrt(np.power(x1[1:L1]-x1[0:L1-1],2) + np.power(y1[1:L1]-y1[0:L1-1],2))
    radius2 = np.sqrt(np.power(x2[1:L2]-x2[0:L2-1],2) + np.power(y2[1:L2]-y2[0:L2-1],2))

    radius = np.concatenate((radius1, radius2), axis=None)
    radius_max = np.max(radius)
    
    for i in range(len(x1)):
        for j in range(len(x2)):
            dist[i,j] = np.sqrt((x1[i] - x2[j])**2 + (y1[i] - y2[j])**2)
    index_min = []

    
    while np.min(dist) < radius_max:
        new_point = 0
        index_min_new = np.unravel_index(dist.argmin(), dist.shape)
        
        for k in range(len(index_min)):
            for i in range(-1, 2):
                if (index_min_new[0] + i == index_min[k][0] or index_min_new[1] + i == index_min[k][1]):
                   new_point = 1
                              
        try: p1 = np.array([x1[index_min_new[0] + 1],y1[index_min_new[0] + 1]]) 
        except IndexError: p1 = np.array([x1[0],y1[0]]) 
        
        p2 = np.array([x1[index_min_new[0] - 1],y1[index_min_new[0] - 1]]) 
        
        try: p3 = np.array([x2[index_min_new[1] + 1],y2[index_min_new[1] + 1]]) 
        except IndexError: p3 = np.array([x2[0],y2[0]]) 
        
        p4 = np.array([x2[index_min_new[1] - 1],y2[index_min_new[1] - 1]]) 
        
        
        if new_point==0 and crossSegment(p1, p2, p3, p4) ==1:            
            index_min.append(index_min_new)
        
        dist[index_min_new[0],index_min_new[1]] = 2*radius_max

    index_min = np.array(index_min)
    try: index_min = index_min[index_min[:,0].argsort()] 
    
    except IndexError: pass
                     
    return index_min                       



def polyShadedPortArea(geo, theta):
   
    # Orbiting scroll curve definition
    # Arc2 
    N_div = 100
    area_port = 0     
    phi_arc2 = np.linspace(geo.t_a22, geo.t_a21, N_div)
    x_arc2 = -geo.x_a2 - geo.r_a2*np.cos(phi_arc2) + geo.r_o*cos(geo.phi_ie_fix - pi/2 - theta)
    y_arc2 = -geo.y_a2 - geo.r_a2*np.sin(phi_arc2) + geo.r_o*sin(geo.phi_ie_fix - pi/2 - theta)    
    
    # Orbiting curve outer
    phi_curve_oo = np.linspace( geo.phi_ie_fix - 2*pi*(geo.N_c_max)   , geo.phi_os, N_div)
    phi_curve_oo = np.linspace( geo.phi_os  + pi, geo.phi_os, N_div)
    (x_curve_oo, y_curve_oo) = geo.getCoordOrbiting( phi_curve_oo, 'o', theta)   
    
    # Arc1
    phi_arc1 = np.linspace(geo.t_a11, geo.t_a12, N_div)
    x_arc1 = -geo.x_a1 - geo.r_a1*np.cos(phi_arc1) + geo.r_o*cos(geo.phi_ie_fix - pi/2 - theta)
    y_arc1 = -geo.y_a1 - geo.r_a1*np.sin(phi_arc1) + geo.r_o*sin(geo.phi_ie_fix - pi/2 - theta)
    
    # Merging
    if geo.dischargeType == 'aala':
         phi_curveArc0 = np.linspace(geo.t_a01, geo.t_a02, N_div)
         x_arc0 = -geo.x_a0 - geo.r_a0*np.cos(phi_curveArc0) + geo.r_o*cos(geo.phi_ie_fix - pi/2 - theta)
         y_arc0 = -geo.y_a0 - geo.r_a0*np.sin(phi_curveArc0) + geo.r_o*sin(geo.phi_ie_fix - pi/2 - theta)

         x_L_f =  np.linspace(geo.x_a2t, geo.x_a1t, N_div)
         x_L = - x_L_f + geo.r_o*cos(geo.phi_ie_fix - pi/2 - theta)
         y_L = - geo.m_l*x_L_f - geo.b_l   + geo.r_o*sin(geo.phi_ie_fix - pi/2 - theta)
         
         x_odis = np.concatenate((x_curve_oo, x_arc2,  x_L, x_arc1, x_arc0), axis=None)
         y_odis = np.concatenate((y_curve_oo, y_arc2,  y_L, y_arc1, y_arc0), axis=None) 
        
    elif geo.dischargeType == 'ala':
         x_L_f =  np.linspace(geo.x_a2t, geo.x_a1t, N_div)
         x_L = - x_L_f + geo.r_o*cos(geo.phi_ie_fix - pi/2 - theta)
         y_L = - geo.m_l*x_L_f - geo.b_l   + geo.r_o*sin(geo.phi_ie_fix - pi/2 - theta)
         
         x_odis = np.concatenate((x_curve_oo, x_arc2,  x_L, x_arc1), axis=None)
         y_odis = np.concatenate((y_curve_oo, y_arc2,  y_L, y_arc1), axis=None)                   

           
    elif geo.dischargeType == '2a':     
         x_odis = np.concatenate((x_curve_oo, x_arc2, x_arc1), axis=None)
         y_odis = np.concatenate((y_curve_oo, y_arc2, y_arc1), axis=None)  
         
               
    x_odis = np.concatenate((x_odis, x_odis[0]), axis=None)
    y_odis = np.concatenate((y_odis, y_odis[0]), axis=None)
       
    (xport, yport) = geo.discharge_port(start = 0)
    P_port = []
    P_odis = []
    
    for i in range(len(xport)):
        P_port.append( (xport[i], yport[i]) )
        
    for i in range(len(x_odis)):
        P_odis.append( (x_odis[i], y_odis[i]) )
    
    
    poly_port = Polygon(P_port)
    poly_odis = Polygon(P_odis)
    shaded_area = poly_port.intersection(poly_odis)
    
    area_port = poly_port.area - shaded_area.area
    
    # if geo.A_dis != False:
    #     area_port = 98304*exp(-2.259*geo.V_ratio)
    #     area_port = 205.7
        
    x_polyShade = 0
    y_polyShade = 0
    
    return x_polyShade, y_polyShade, area_port
