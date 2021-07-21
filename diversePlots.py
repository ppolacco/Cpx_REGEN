# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 16:44:14 2021

@author: nicolas
"""


import GeometryFunctions as geoF
from math import pi
import numpy as np

import matplotlib.pyplot as plt
import solverFunctions as solverF
from scipy.interpolate import interp1d



# The 'scrollPlot' function takes as input the object 'geo' that contains all the information related to the 
# machine geometry and the set of angles 'theta'. The output of the fucntion is the plot of the geometry specified 
# in 'geo'
def scrollPlot(geo,theta):

    geo.getCoordFullGeometry(theta)
    
    plt.figure()
    ax = plt.axes(xlim = (-70,70),ylim = (-70,70))
    
    (xport, yport) = geo.discharge_port()
    plt.plot(xport, yport,'k', zorder=1)
    plt.plot(geo.x_fscroll,geo.y_fscroll,'k',geo.x_oscroll, geo.y_oscroll, 'k')
    plt.fill(geo.x_fscroll, geo.y_fscroll, color='orange')
    plt.fill(geo.x_oscroll, geo.y_oscroll, color='grey', zorder=2)
    plt.plot(geo.x_wall,geo.y_wall,'k')

    # Plot forces 

    # scale = 10
    # f_s1 = geoF.volume_force_s1(geo,theta)["f_vect"]
    # f_s2 = geoF.volume_force_s2(geo,theta)["f_vect"]
    # f_sa = geoF.volume_force_sa(geo,theta)["f_vect"]
    
    # c1 = geoF.volume_force_c1(geo,theta)
    # c2 = geoF.volume_force_c2(geo,theta)
    # f_d1 = geoF.volume_force_d1(geo,theta)["f_vect"]
    # f_d2 = geoF.volume_force_d2(geo,theta)["f_vect"]
    # f_dd = geoF.volume_force_dd(geo,theta)["f_vect"]    
    
    # for i in range(0,len(f_s1[0,:]),3):
    #     plt.arrow(f_s1[0,i], f_s1[1,i], f_s1[2,i]*scale, f_s1[3,i]*scale, width = 0.3,zorder=3)   
        
    # for i in range(0,len(f_s2[0,:]),3):
    #     plt.arrow(f_s2[0,i], f_s2[1,i], f_s2[2,i]*scale, f_s2[3,i]*scale, width = 0.3,zorder=3)     
        
    # for i in range(0,len(f_sa[0,:]),3):
    #     plt.arrow(f_sa[0,i], f_sa[1,i], f_sa[2,i]*scale, f_sa[3,i]*scale, width = 0.3,zorder=3) 
        
    # for j in range(len(c1)):
    #     for i in range(0,len(c1[j+1]["f_vect"][0,:]),3):
    #         plt.arrow(c1[j+1]["f_vect"][0,i], c1[j+1]["f_vect"][1,i], c1[j+1]["f_vect"][2,i]*scale, c1[j+1]["f_vect"][3,i]*scale, width = 0.3,zorder=3) 
            
    #     for i in range(0,len(c2[j+1]["f_vect"][0,:]),3):
    #         plt.arrow(c2[j+1]["f_vect"][0,i], c2[j+1]["f_vect"][1,i], c2[j+1]["f_vect"][2,i]*scale, c2[j+1]["f_vect"][3,i]*scale, width = 0.3,zorder=3) 
        
        
    # for i in range(0,len(f_d1[0,:]),3):
    #     plt.arrow(f_d1[0,i], f_d1[1,i], f_d1[2,i]*scale, f_d1[3,i]*scale, width = 0.3,zorder=3)  
        
    # for i in range(0,len(f_d2[0,:]),3):
    #     plt.arrow(f_d2[0,i], f_d2[1,i], f_d2[2,i]*scale, f_d2[3,i]*scale, width = 0.3,zorder=3) 
        
    # for i in range(0,len(f_dd[0,:]),15):
    #     plt.arrow(f_dd[0,i], f_dd[1,i], f_dd[2,i]*scale, f_dd[3,i]*scale, width = 0.3,zorder=3) 
    
    
    
    ## Plot chambers
    
    
    # geo.important_values()
    # xy_polyS1 = geoF.volume_force_s1(geo,theta)["xy_poly"]
    # xy_polyS2 = geoF.volume_force_s2(geo,theta)["xy_poly"]
    # xy_polySA1 = geoF.volume_force_SA(geo,theta)["xy_polySA1"]
    # xy_polySA2 = geoF.volume_force_SA(geo,theta)["xy_polySA2"]
    
    # xy_polyC11 = geoF.volume_force_c1(geo,theta,geo.N_c_max)["xy_poly"]
    # xy_polyC12 = geoF.volume_force_c1(geo,theta,1)["xy_poly"]
    # xy_polyC21 = geoF.volume_force_c2(geo,theta,geo.N_c_max)["xy_poly"]
    # xy_polyC22 = geoF.volume_force_c2(geo,theta,1)["xy_poly"]
    
    #xy_polyD1 = geoF.volume_force_d1(geo,theta)["xy_poly"]
    # xy_polyD2 = geoF.volume_force_d2(geo,theta)["xy_poly"]
    # xy_polyDD = geoF.volume_force_dd(geo,theta)["xy_poly"]

    
    # plt.plot(xy_polyS1[:,0], xy_polyS1[:,1],'b')
    # plt.plot(xy_polyS2[:,0], xy_polyS2[:,1],'b')
    # plt.plot(xy_polySA1[:,0], xy_polySA1[:,1],'b')
    # plt.plot(xy_polySA2[:,0], xy_polySA2[:,1],'b')
    # plt.plot(xy_polyC11[:,0], xy_polyC11[:,1],'b')
    # plt.plot(xy_polyC21[:,0], xy_polyC21[:,1],'b') 
    # plt.plot(xy_polyC12[:,0], xy_polyC12[:,1],'b')
    # plt.plot(xy_polyC22[:,0], xy_polyC22[:,1],'b') 
    #plt.plot(xy_polyD1[:,0], xy_polyD1[:,1],'b')
    # plt.plot(xy_polyD2[:,0], xy_polyD2[:,1],'b')
    # plt.plot(xy_polyDD[:,0], xy_polyDD[:,1],'b')
    
    ## Plot discharge port areas 
    
    # (x_suction,y_suction,_) = suctionFlowArea(geo, theta)
    # (x_d_dd,y_d_dd,_) = area_d_dd(geo,theta)
    # (x_polyShade,y_polyShade,_) = polyShadedPortArea(geo, theta)
    
    # plt.plot(x_suction, y_suction,'b')
    # plt.plot(x_d_dd, y_d_dd,'b')
    # plt.plot(x_polyShade, y_polyShade,'b',linewidth=2.5)
    

# The function 'plot_V_dV' takes as input the object 'geo' that contains all the information related to the 
# machine geometry and it plots chambers' volume and volume derivatives versus the orbithing angle.
def plot_V_dV(geo):
     
     # The orbiting angle [0,2pi] is dicretized using 90 steps
     theta = np.linspace(0,2*pi,90)
     
     # Initialization of the vectors in which volume size will be collected
     V_s1_vect = np.zeros(len(theta))
     V_s2_vect = np.zeros(len(theta))
     V_c1_vect = np.zeros((len(theta),geo.N_c_max))
     V_c2_vect = np.zeros((len(theta),geo.N_c_max))
     V_d1_vect = np.zeros(len(theta))
     V_d2_vect = np.zeros(len(theta))
     V_dd_vect = np.zeros(len(theta))
     V_dd_vect = np.zeros(len(theta))
     V_ddd_vect = np.zeros(len(theta))
     
     # Initialization of the vectors in which volume derivatives will be collected
     dV_s1_vect = np.zeros(len(theta))
     dV_s2_vect = np.zeros(len(theta))
     dV_c1_vect = np.zeros((len(theta),geo.N_c_max))
     dV_c2_vect = np.zeros((len(theta),geo.N_c_max))
     dV_d1_vect = np.zeros(len(theta))
     dV_d2_vect = np.zeros(len(theta))
     dV_dd_vect = np.zeros(len(theta))
     dV_dd_vect = np.zeros(len(theta))
     dV_ddd_vect = np.zeros(len(theta))
     
     # The 'volume_force_*' functions returns A dictionnary containing the volume, volume derivative,
     # coordinate of the centroid and orbiting moment generated by the chamber *. There is a function for 
     # any type of chambers inside a scroll machine.
     # Chambers are evalutad at all the considered orbiting angles 'theta'
     for i in range(len(theta)): 
         # COMPUTATION OF VOLUMES
         V_s1_vect[i] = geoF.volume_force_s1(geo,theta[i])["V"]
         V_s2_vect[i] = geoF.volume_force_s2(geo,theta[i])["V"]
         # The compression chambers are in number of 'alpha'
         for alpha in range(1,geo.N_c_max + 1,1):
             V_c1_vect[i,alpha-1] = geoF.volume_force_c1(geo,theta[i])[alpha]["V"]
             V_c2_vect[i,alpha-1] = geoF.volume_force_c2(geo,theta[i])[alpha]["V"]
         V_d1_vect[i] = geoF.volume_force_d1(geo,theta[i])["V"]
         V_d2_vect[i] = geoF.volume_force_d2(geo,theta[i])["V"]
         V_dd_vect[i] = geoF.volume_force_dd(geo,theta[i])["V"]
         # COMPUTATION OF VOLUME DERIVATIVES
         dV_s1_vect[i] = geoF.volume_force_s1(geo,theta[i])["dVdTheta"]
         dV_s2_vect[i] = geoF.volume_force_s2(geo,theta[i])["dVdTheta"]
         for alpha in range(1,geo.N_c_max + 1,1):        
             dV_c1_vect[i,alpha-1] = geoF.volume_force_c1(geo,theta[i])[alpha]["dVdTheta"]
             dV_c2_vect[i,alpha-1] = geoF.volume_force_c2(geo,theta[i])[alpha]["dVdTheta"]
         dV_d1_vect[i] = geoF.volume_force_d1(geo,theta[i])["dVdTheta"]
         dV_d2_vect[i] = geoF.volume_force_d2(geo,theta[i])["dVdTheta"]
         dV_dd_vect[i] = geoF.volume_force_dd(geo,theta[i])["dVdTheta"]
         
         if (geo.theta_d >= geo.theta_u and (theta[i] >= geo.theta_u and theta[i] <= geo.theta_d)) or (geo.theta_u >= geo.theta_d and (theta[i] >= geo.theta_u or theta[i] <= geo.theta_d)): 
            V_ddd_vect[i] = V_d1_vect[i] + V_d2_vect[i] + V_dd_vect[i]
            V_d1_vect[i] = np.nan
            V_d2_vect[i] = np.nan
            V_dd_vect[i] = np.nan
            
            dV_ddd_vect[i] = dV_d1_vect[i] + dV_d2_vect[i] + dV_dd_vect[i]
            dV_d1_vect[i] = np.nan
            dV_d2_vect[i] = np.nan
            dV_dd_vect[i] = np.nan 
         else:
            V_ddd_vect[i] = np.nan
            dV_ddd_vect[i] = np.nan
    
     
     line_vs = (V_s1_vect + V_s2_vect)/1000
     line_vc1 = (V_c1_vect[:,0] + V_c2_vect[:,0])/1000
     line_vc2 = (V_c1_vect[:,1] + V_c2_vect[:,1])/1000
     line_vd = (V_d1_vect + V_d2_vect)/1000
     line_vdd = (V_dd_vect )/1000
     line_vddd = (V_ddd_vect )/1000
     
     line_dvs = (dV_s1_vect + dV_s2_vect)/1000
     line_dvc1 = (dV_c1_vect[:,0]  + dV_c2_vect[:,0] )/1000
     line_dvc2 = (dV_c1_vect[:,1] + dV_c2_vect[:,1])/1000
     line_dvd = (dV_d1_vect + dV_d2_vect)/1000
     line_dvdd = (dV_dd_vect )/1000
     line_dvddd = (dV_ddd_vect )/1000
     
     # Definition of the volume plot parameters and frame
     #Volume plot
     plt.figure(figsize=(6,4.5),constrained_layout=True)
     plt.xlim(theta.min(), theta.max())
     plt.xticks([0, pi/2, pi, 3*pi/2, 2*pi],
       [r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
     
     # Plot of the volumes
     plt.ylim(0,line_vs.max()*1.1)
     #plt.yticks([-1, 0, +1],
     #   [r'$-1$', r'$0$', r'$+1$'])
     plt.grid()
     plt.plot(theta, line_vs, 'c' , linewidth = 2 )
     plt.plot(theta, line_vc1, 'r', linewidth = 2 )
     plt.plot(theta, line_vc2, 'r', linewidth = 2 )
     plt.plot(theta, line_vd, 'g', linewidth = 2)
     plt.plot(theta, line_vdd, 'k', linewidth = 2 )
     plt.plot(theta, line_vddd, 'y', linewidth = 2 )
     plt.text(5, 90,'$V_{s1} + V_{s2}$',fontsize=14)
     plt.text(1, 80,'$V_{c1} + V_{c2}$',fontsize=14)
     plt.text(1.5, 38,'$V_{d1} + V_{d2}$',fontsize=14)
     plt.text(2, 6,'$V_{dd}$',fontsize=14)
     plt.text(5, 15,'$V_{ddd}$',fontsize=14)
     # Definition of plot lables
     plt.xlabel('Orbiting angle [rad]',fontsize=16,  fontname="Times New Roman")
     plt.ylabel('Volume [cm$^3$]',fontsize=16, fontname="Times New Roman")   
         
     # Definition of the volume derivative plot parameters and frame
     plt.figure(figsize=(6,4.5))
     plt.xlim(theta.min(), theta.max())
     plt.xticks([0, pi/2, pi, 3*pi/2, 2*pi],
       [r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
     
     # Plot of volume derivatives
     plt.ylim(-20,line_dvs.max()*1.1)
     #plt.yticks([-1, 0, +1],
     #   [r'$-1$', r'$0$', r'$+1$'])
     plt.grid()
     plt.plot(theta, line_dvs, 'c' , linewidth = 2 )
     plt.plot(theta, line_dvc1, 'r', linewidth = 2 )
     plt.plot(theta, line_dvc2, 'r', linewidth = 2 )
     plt.plot(theta, line_dvd, 'g', linewidth = 2)
     plt.plot(theta, line_dvdd, 'k', linewidth = 2 )
     plt.plot(theta, line_dvddd, 'y', linewidth = 2 )
     plt.text(0.4, 23, r'$ \dfrac{\mathrm{d}V_{s1}}{\mathrm{d} \theta} +  \dfrac{\mathrm{d}V_{s2}}{\mathrm{d} \theta}$',fontsize=14)
     plt.text(0.1, -12, r'$ \dfrac{\mathrm{d}V_{c1}}{\mathrm{d} \theta} +  \dfrac{\mathrm{d}V_{c2}}{\mathrm{d} \theta}$',fontsize=14)
     plt.text(1.8, -15, r'$ \dfrac{\mathrm{d}V_{d1}}{\mathrm{d} \theta} +  \dfrac{\mathrm{d}V_{d2}}{\mathrm{d} \theta}$',fontsize=14)
     plt.text(1.5, 4.6, r'$ \dfrac{\mathrm{d}V_{dd}}{\mathrm{d} \theta} $',fontsize=14)
     plt.text(4.2, -2, r'$ \dfrac{\mathrm{d}V_{ddd}}{\mathrm{d} \theta} $',fontsize=14)
     # Definition of plot lables
     plt.xlabel('Orbiting angle [rad]',fontsize=16,  fontname="Times New Roman")
     plt.ylabel('Volume Derivative [cm$^3$/rad]',fontsize=16, fontname="Times New Roman")       
         
     
def plot_evolution(geo, CV, P_out, P_in):
    
    
     CV_s1   = CV['s1']
     CV_c1   = CV['c1']
     CV_d1   = CV['d1']
     CV_dd   = CV['dd']
     CV_ddd   = CV['ddd']
    
     (vect_s1, vect_c1, vect_d1, vect_dd, vect_ddd) = solverF.convert_recorded_values(geo, CV_s1, CV_c1, CV_d1, CV_dd, CV_ddd)
     
     
     
     xdim = 5.5
     ydim = 4.5
     #%% Pressure 
     prop = 'P'
     
     plt.figure(dpi=100,figsize=(xdim,ydim),constrained_layout=True)
     plt.xlim(vect_s1['theta'].min(), vect_s1['theta'].max())
     plt.xticks([0, pi/2, pi, 3*pi/2, 2*pi],
       [r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])

     plt.ylim(np.nanmin(vect_s1[prop]/1e5)/1.1, np.nanmax(vect_ddd[prop]/1e5)*1.1)
     #plt.yticks([-1, 0, +1],
     #   [r'$-1$', r'$0$', r'$+1$'])
     plt.grid()
     
     
     plt.plot(vect_s1['theta'], vect_s1[prop]/1e5 , linewidth = 2, label ='s1' )
     alpha = 1
     while alpha <= geo.N_c_max:
         if alpha ==  geo.N_c_max:
             plt.plot(vect_s1['theta'], vect_c1[alpha][prop]/1e5 , 'r', linewidth = 2, label ='c1' )    
         else :
             plt.plot(vect_s1['theta'], vect_c1[alpha][prop]/1e5 , 'r', linewidth = 2 )   
         alpha += 1
        

     plt.plot(vect_s1['theta'], vect_d1[prop]/1e5 , 'g',  linewidth = 2, label ='d1' )
     plt.plot(vect_s1['theta'], vect_dd[prop]/1e5 , 'k', linewidth = 2, label ='dd' )
     plt.plot(vect_s1['theta'], vect_ddd[prop]/1e5 ,'y', linewidth = 2, label ='ddd' )

     plt.plot( vect_s1['theta']  ,  P_out/1e5*np.ones(len(vect_s1['theta'])), 'r--' , label = r'$P_{out}$')
     plt.plot( vect_s1['theta']  ,  P_in/1e5*np.ones(len(vect_s1['theta'])), 'b--' , label = r'$P_{in}$')
     
     plt.xlabel('Orbiting angle [rad]',fontsize=16,  fontname="Times New Roman")
     plt.ylabel('Pressure [bar]',fontsize=16, fontname="Times New Roman")            
     
     plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3, fontsize=14)
     
     #%% Temperature
     
     prop = 'T'
     
     plt.figure(dpi=100,figsize=(xdim,ydim),constrained_layout=True)
     plt.xlim(vect_s1['theta'].min(), vect_s1['theta'].max())
     plt.xticks([0, pi/2, pi, 3*pi/2, 2*pi],
       [r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])

     plt.ylim(np.nanmin(vect_s1[prop])/1.01, np.nanmax(vect_ddd[prop])*1.01)
     #plt.yticks([-1, 0, +1],
     #   [r'$-1$', r'$0$', r'$+1$'])
     plt.grid()
     
     
     plt.plot(vect_s1['theta'], vect_s1[prop] , linewidth = 2, label ='s1' )
     alpha = 1
     while alpha <= geo.N_c_max:
         if alpha ==  geo.N_c_max:
             plt.plot(vect_s1['theta'], vect_c1[alpha][prop] , 'r', linewidth = 2, label ='c1' )    
         else :
             plt.plot(vect_s1['theta'], vect_c1[alpha][prop] , 'r', linewidth = 2 )   
         alpha += 1
        

     plt.plot(vect_s1['theta'], vect_d1[prop] , 'g',  linewidth = 2, label ='d1' )
     plt.plot(vect_s1['theta'], vect_dd[prop] , 'k', linewidth = 2, label ='dd' )
     plt.plot(vect_s1['theta'], vect_ddd[prop] ,'y', linewidth = 2, label ='ddd' )


     
     plt.xlabel('Orbiting angle [rad]',fontsize=16,  fontname="Times New Roman")
     plt.ylabel('Temperature [K]',fontsize=16, fontname="Times New Roman")            
     
     plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3, fontsize=14)     
           
     #%% entropy
     
     prop = 's'
     
     plt.figure(figsize=(xdim,ydim),constrained_layout=True)
     plt.xlim(vect_s1['theta'].min(), vect_s1['theta'].max())
     plt.xticks([0, pi/2, pi, 3*pi/2, 2*pi],
       [r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])

     plt.ylim(np.nanmin(vect_s1[prop])/1.001, np.nanmax(vect_dd[prop])*1.001)
     #plt.yticks([-1, 0, +1],
     #   [r'$-1$', r'$0$', r'$+1$'])
     plt.grid()
     
     
     plt.plot(vect_s1['theta'], vect_s1[prop] , linewidth = 2, label ='s1' )
     alpha = 1
     while alpha <= geo.N_c_max:
         if alpha ==  geo.N_c_max:
             plt.plot(vect_s1['theta'], vect_c1[alpha][prop] , 'r', linewidth = 2, label ='c1' )    
         else :
             plt.plot(vect_s1['theta'], vect_c1[alpha][prop] , 'r', linewidth = 2 )   
         alpha += 1
        

     plt.plot(vect_s1['theta'], vect_d1[prop] , 'g',  linewidth = 2, label ='d1' )
     plt.plot(vect_s1['theta'], vect_dd[prop] , 'k', linewidth = 2, label ='dd' )
     plt.plot(vect_s1['theta'], vect_ddd[prop] ,'y', linewidth = 2, label ='ddd' )


     
     plt.xlabel('Orbiting angle [rad]',fontsize=16,  fontname="Times New Roman")
     plt.ylabel('Entropy [J/(kgK)]',fontsize=16, fontname="Times New Roman")            
     
     plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3, fontsize=14)          
     
     #%% Quality
     
     prop = 'Q'
     
     plt.figure(figsize=(xdim,ydim),constrained_layout=True)
     plt.xlim(vect_s1['theta'].min(), vect_s1['theta'].max())
     plt.xticks([0, pi/2, pi, 3*pi/2, 2*pi],
       [r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])

     plt.ylim(0, 1)
     #plt.yticks([-1, 0, +1],
     #   [r'$-1$', r'$0$', r'$+1$'])
     plt.grid()
     
     
     plt.plot(vect_s1['theta'], vect_s1[prop] , linewidth = 2, label ='s1' )
     alpha = 1
     while alpha <= geo.N_c_max:
         if alpha ==  geo.N_c_max:
             plt.plot(vect_s1['theta'], vect_c1[alpha][prop] , 'r', linewidth = 2, label ='c1' )    
         else :
             plt.plot(vect_s1['theta'], vect_c1[alpha][prop] , 'r', linewidth = 2 )   
         alpha += 1
        

     plt.plot(vect_s1['theta'], vect_d1[prop] , 'g',  linewidth = 2, label ='d1' )
     plt.plot(vect_s1['theta'], vect_dd[prop] , 'k', linewidth = 2, label ='dd' )
     plt.plot(vect_s1['theta'], vect_ddd[prop] ,'y', linewidth = 2, label ='ddd' )


     
     plt.xlabel('Orbiting angle [rad]',fontsize=16,  fontname="Times New Roman")
     plt.ylabel('Quality [-]',fontsize=16, fontname="Times New Roman")            
     
     plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3, fontsize=14)          

     #%% Mass 
     
     prop = 'm'
     
     plt.figure(figsize=(xdim,ydim),constrained_layout=True)
     plt.xlim(vect_s1['theta'].min(), vect_s1['theta'].max())
     plt.xticks([0, pi/2, pi, 3*pi/2, 2*pi],
       [r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])

     plt.ylim(np.nanmin(vect_s1[prop])/1.1, np.nanmax(vect_ddd[prop])*1.1)
     #plt.yticks([-1, 0, +1],
     #   [r'$-1$', r'$0$', r'$+1$'])
     plt.grid()
     
     
     plt.plot(vect_s1['theta'], vect_s1[prop] , linewidth = 2, label ='s1' )
     alpha = 1
     while alpha <= geo.N_c_max:
         if alpha ==  geo.N_c_max:
             plt.plot(vect_s1['theta'], vect_c1[alpha][prop] , 'r', linewidth = 2, label ='c1' )    
         else :
             plt.plot(vect_s1['theta'], vect_c1[alpha][prop] , 'r', linewidth = 2 )   
         alpha += 1
        

     plt.plot(vect_s1['theta'], vect_d1[prop] , 'g',  linewidth = 2, label ='d1' )
     plt.plot(vect_s1['theta'], vect_dd[prop] , 'k', linewidth = 2, label ='dd' )
     plt.plot(vect_s1['theta'], vect_ddd[prop] ,'y', linewidth = 2, label ='ddd' )


     
     plt.xlabel('Orbiting angle [rad]',fontsize=16,  fontname="Times New Roman")
     plt.ylabel('Mass [kg]',fontsize=16, fontname="Times New Roman")            
     
     plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3, fontsize=14)               


     #%% x_l
     
     prop = 'x_l'
     
     plt.figure(figsize=(xdim,ydim),constrained_layout=True)
     plt.xlim(vect_s1['theta'].min(), vect_s1['theta'].max())
     plt.xticks([0, pi/2, pi, 3*pi/2, 2*pi],
       [r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])

     plt.ylim(-0.1, 0.3)
     #plt.yticks([-1, 0, +1],
     #   [r'$-1$', r'$0$', r'$+1$'])
     plt.grid()
     
     
     plt.plot(vect_s1['theta'], vect_s1[prop] , linewidth = 2, label ='s1' )
     alpha = 1
     while alpha <= geo.N_c_max:
         if alpha ==  geo.N_c_max:
             plt.plot(vect_s1['theta'], vect_c1[alpha][prop] , 'r', linewidth = 2, label ='c1' )    
         else :
             plt.plot(vect_s1['theta'], vect_c1[alpha][prop] , 'r', linewidth = 2 )   
         alpha += 1
        

     plt.plot(vect_s1['theta'], vect_d1[prop] , 'g',  linewidth = 2, label ='d1' )
     plt.plot(vect_s1['theta'], vect_dd[prop] , 'k', linewidth = 2, label ='dd' )
     plt.plot(vect_s1['theta'], vect_ddd[prop] ,'y', linewidth = 2, label ='ddd' )


     
     plt.xlabel('Orbiting angle [rad]',fontsize=16,  fontname="Times New Roman")
     plt.ylabel('Oil rate [-]',fontsize=16, fontname="Times New Roman")            
     
     plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3, fontsize=14)                 

#%% rho
     
     prop = 'rho'
     
     plt.figure(figsize=(xdim,ydim),constrained_layout=True)
     plt.xlim(vect_s1['theta'].min(), vect_s1['theta'].max())
     plt.xticks([0, pi/2, pi, 3*pi/2, 2*pi],
       [r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
     
     plt.ylim(np.nanmin(vect_s1[prop])/1.1, np.nanmax(vect_dd[prop])*1.1)
     #plt.yticks([-1, 0, +1],
     #   [r'$-1$', r'$0$', r'$+1$'])
     plt.grid()
     
     
     plt.plot(vect_s1['theta'], vect_s1[prop] , linewidth = 2, label ='s1' )
     alpha = 1
     while alpha <= geo.N_c_max:
         if alpha ==  geo.N_c_max:
             plt.plot(vect_s1['theta'], vect_c1[alpha][prop] , 'r', linewidth = 2, label ='c1' )    
         else :
             plt.plot(vect_s1['theta'], vect_c1[alpha][prop] , 'r', linewidth = 2 )   
         alpha += 1
        

     plt.plot(vect_s1['theta'], vect_d1[prop] , 'g',  linewidth = 2, label ='d1' )
     plt.plot(vect_s1['theta'], vect_dd[prop] , 'k', linewidth = 2, label ='dd' )
     plt.plot(vect_s1['theta'], vect_ddd[prop] ,'y', linewidth = 2, label ='ddd' )

     
     plt.xlabel('Orbiting angle [rad]',fontsize=16,  fontname="Times New Roman")
     plt.ylabel('Density [kg/m$^3$]',fontsize=16, fontname="Times New Roman")            
     
     plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3, fontsize=14)      


#%% m_dot_su 
     prop = 'm_dot'
     
     plt.figure(dpi=100,figsize=(xdim,ydim),constrained_layout=True)
     plt.xlim(vect_s1['theta'].min(), vect_s1['theta'].max())
     plt.xticks([0, pi/2, pi, 3*pi/2, 2*pi],
       [r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])

     plt.ylim(-0.025, np.nanmax(vect_s1[prop]['su'])*1.1*2)
     #plt.yticks([-1, 0, +1],
     #   [r'$-1$', r'$0$', r'$+1$'])
     plt.grid()
     
     
     plt.plot(vect_s1['theta'], vect_s1[prop]['su']*2 , linewidth = 2,  )
     

    
     
     plt.xlabel('Orbiting angle [rad]',fontsize=16,  fontname="Times New Roman")
     plt.ylabel('Inlet mass flow rate [kg/s]',fontsize=16, fontname="Times New Roman")            



  #%% m_dot_ex
     prop = 'm_dot'
     
     plt.figure(dpi=100,figsize=(xdim,ydim),constrained_layout=True)
     plt.xlim(vect_s1['theta'].min(), vect_s1['theta'].max())
     plt.xticks([0, pi/2, pi, 3*pi/2, 2*pi],
       [r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])

     plt.ylim(np.nanmin(-vect_ddd[prop]['ex'])*1.1, np.nanmax(-vect_ddd[prop]['ex'])*1.1)
     #plt.yticks([-1, 0, +1],
     #   [r'$-1$', r'$0$', r'$+1$'])
     plt.grid()
     
     
     plt.plot(vect_s1['theta'], -vect_ddd['m_dot']['ex'] , linewidth = 2, label ='ddd' )
     plt.plot(vect_s1['theta'], -vect_dd['m_dot']['ex'] , linewidth = 2, label ='dd' )

     
     plt.xlabel('Orbiting angle [rad]', fontsize=16,  fontname="Times New Roman")
     plt.ylabel('Outlet mass flow rate [kg/s]', fontsize=16, fontname="Times New Roman")            
     

     plt.legend(fontsize=16, loc='best')
     
  # #%% m_dot_d_dd
     prop = 'm_dot'
     
     plt.figure(dpi=100,figsize=(xdim,ydim),constrained_layout=True)
     plt.xlim(vect_s1['theta'].min(), vect_s1['theta'].max())
     plt.xticks([0, pi/2, pi, 3*pi/2, 2*pi],
       [r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])

     plt.ylim(np.nanmin(-vect_d1[prop]['fl2'])*1.1, np.nanmax(-vect_d1[prop]['fl2'])*1.1)
      #plt.yticks([-1, 0, +1],
      #   [r'$-1$', r'$0$', r'$+1$'])
     plt.grid()
     
     
     plt.plot(vect_s1['theta'], -vect_d1['m_dot']['fl2'] , linewidth = 2, label ='ddd' )
    

     
     plt.xlabel('Orbiting angle [rad]', fontsize=16,  fontname="Times New Roman")
     plt.ylabel('Merging d-dd flow rate [kg/s]', fontsize=16, fontname="Times New Roman")            
     

     plt.legend(fontsize=16, loc='best')

def basic_plot(x, y, name_x, name_y, x2 = [None], y2 = [None]):    
    
    
    f = interp1d(x, y, kind='cubic')
    x_new = np.linspace(x.min(),x.max(),100, endpoint=True)
    plt.figure(figsize=(7,4.5),constrained_layout=True)
    d_x = x.max() - x.min()
    
    
    if x2[0] != None: 
        d_y = y.max() - y2.min()
        plt.ylim(y2.min() - 0.1*d_y, y.max() + 0.1*d_y)
    else : 
        d_y = y.max() - y.min() 
        plt.ylim(y.min() - 0.1*d_y, y.max() + 0.1*d_y)
    plt.xlim(x.min() - 0.1*d_x, x.max() + 0.1*d_x)
    
    
    
    plt.grid()
    
    plt.xlabel(name_x, fontsize=16,  fontname="Times New Roman")
    plt.ylabel(name_y ,fontsize=16, fontname="Times New Roman")    
    plt.plot(x, y ,  'o')
    plt.plot(x_new, f(x_new) , linewidth = 2)
    
    if x2[0] != None:
        f2 = interp1d(x2, y2, kind='cubic')
        x2_new = np.linspace(x2.min(),x2.max(),100, endpoint=True)
        plt.plot(x2, y2 ,  'o')
        plt.plot(x2_new, f2(x2_new) , linewidth = 2)
    
    plt.legend(['data', 'curve'],fontsize=16, loc='best')
    
    


 
    
def plot_volumes(geo, CV):
    
   CV_s1   = CV['s1']
   CV_c1   = CV['c1']
   CV_d1   = CV['d1']
   CV_dd   = CV['dd']
   CV_ddd  = CV['ddd']
    
   (vect_s1, vect_c1, vect_d1, vect_dd, vect_ddd) = solverF.convert_recorded_values(geo, CV_s1, CV_c1, CV_d1, CV_dd, CV_ddd)
        
    
   theta_vol = np.linspace(0,2*pi,361)
   

   V_s1_vect = np.zeros(len(theta_vol))
   V_c1_vect = np.zeros((len(theta_vol),geo.N_c_max))
   V_d1_vect = np.zeros(len(theta_vol))
   V_dd_vect = np.zeros(len(theta_vol))
   V_ddd_vect = np.zeros(len(theta_vol))

    
   for i in range(len(theta_vol)): 
       
      
       V_s1_vect[i] = geoF.volume_force_s1(geo,theta_vol[i])["V"]/1e3
      
       
       for alpha in range(1,geo.N_c_max + 1,1):
           V_c1_vect[i,alpha-1] = geoF.volume_force_c1(geo,theta_vol[i])[alpha]["V"]/1e3
           
       V_d1_vect[i] = geoF.volume_force_d1(geo,theta_vol[i])["V"]/1e3
       V_dd_vect[i] = geoF.volume_force_dd(geo,theta_vol[i])["V"]/1e3
       
       if (geo.theta_d >= geo.theta_m and (theta_vol[i] >= geo.theta_m and theta_vol[i] <= geo.theta_m)) or (geo.theta_m >= geo.theta_d and (theta_vol[i] >= geo.theta_m or theta_vol[i] <= geo.theta_d)): 
           V_ddd_vect[i] = V_d1_vect[i] + V_d1_vect[i] + V_dd_vect[i]
           V_d1_vect[i] = np.nan
         
           V_dd_vect[i] = np.nan
          
          
       else:
           V_ddd_vect[i] = np.nan
           
   
   prop = 'V'
       
   plt.figure(figsize=(6,4.5),constrained_layout=True)
   plt.xlim(vect_s1['theta'].min(), vect_s1['theta'].max())
   plt.xticks([0, pi/2, pi, 3*pi/2, 2*pi],
      [r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
    
   plt.ylim(np.nanmin(vect_dd[prop]*1e6)/1.1, np.nanmax(2*vect_s1[prop]*1e6)*1.1)
    #plt.yticks([-1, 0, +1],
    #   [r'$-1$', r'$0$', r'$+1$'])
   plt.grid()
    
    
   plt.plot(vect_s1['theta'], 2*vect_s1[prop]*1e6 , linewidth = 2, label ='s\'')
   plt.plot(theta_vol, 2*V_s1_vect ,u'#1f77b4',linestyle = '--', linewidth = 2, label ='s')
   alpha = 1
   while alpha <= geo.N_c_max:
        if alpha ==  geo.N_c_max:
            plt.plot(vect_s1['theta'], 2*vect_c1[alpha][prop]*1e6 , 'r', linewidth = 2, label ='c\'' )    
            plt.plot(theta_vol, 2*V_c1_vect[:,0] ,'r',linestyle = '--', linewidth = 2, label ='c')
        else :
            plt.plot(vect_s1['theta'], 2*vect_c1[alpha][prop]*1e6 , 'r', linewidth = 2 )   
            plt.plot(theta_vol, 2*V_c1_vect[:,0] ,'r',linestyle = '--', linewidth = 2)
        alpha += 1
       
    
   plt.plot(vect_s1['theta'], 2*vect_d1[prop]*1e6 , 'g',  linewidth = 2, label ='d\'' )
   plt.plot(theta_vol, 2*V_d1_vect ,'g',linestyle = '--', linewidth = 2, label ='d')
    
   plt.plot(vect_s1['theta'], vect_dd[prop]*1e6 , 'k', linewidth = 2, label ='dd\'' )
   plt.plot(theta_vol, V_dd_vect ,'k',linestyle = '--', linewidth = 2, label ='dd')
    
   plt.plot(vect_s1['theta'], vect_ddd[prop]*1e6 ,'y', linewidth = 2, label ='ddd\'' )
   plt.plot(theta_vol, V_ddd_vect ,'y',linestyle = '--', linewidth = 2, label ='ddd')
    
    
   plt.xlabel('Orbiting angle [rad]',fontsize=16,  fontname="Times New Roman")
   plt.ylabel('Volume [cm$^3$]',fontsize=16, fontname="Times New Roman")            
    
   plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
               mode="expand", borderaxespad=0, ncol=3)     



   
                
             