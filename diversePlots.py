# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 16:44:14 2021

@author: nicolas
"""


import GeometryFunctions as geoF
from math import pi
import numpy as np
import matplotlib.animation as animation

import matplotlib.pyplot as plt
import solverFunctions as solverF
from scipy.interpolate import interp1d
from scipy import integrate as intg



def scrollAnim(geo, mode = "comp", save = 0):
    
    """ 
    Function allowing to animate the movement of the scrolls. 
    - geo is the class containing the geometry
    - mode is a string containing either "comp" in compressor mode or "exp" in expoander mode
    - save is integer, if == 0 no save of the animation, if == 1, saving the animation
    Warning: specify the path to save the animation in the funtion!
    
    """
    
    theta = 0
    dtheta = 5*2*pi/360
    geo.getCoordFullGeometry(theta)
    
    fig = plt.figure()
    ax = plt.axes(xlim = (-70,70),ylim = (-70,70))
    (xport, yport) = geo.discharge_port()
    plt.plot(xport, yport,'k', zorder=0)
    plt.plot(geo.x_fscroll,geo.y_fscroll,'k')
    plt.fill(geo.x_fscroll, geo.y_fscroll, color='orange')
    #patch = patches.Polygon(pol,closed=True, fc='r', ec='r')
    plt.plot(geo.x_wall,geo.y_wall,'k')
    
    point = geo.vect_oxy
    
    patch = plt.Polygon(point, fc='grey')
    line, = plt.plot([],[],color='k')
     
    def init():
        
        ax.add_patch(patch)
        line.set_data([],[])
        return patch, line,
    
    def animate(i):
        if mode == "exp": 
            theta = 2*pi - i* dtheta
        else:
            theta = i* dtheta
        geo.getCoordFullGeometry(theta)
        patch.set_xy(geo.vect_oxy)
        line.set_data(geo.x_oscroll,geo.y_oscroll)
    
        return patch, line,
    
    ani = animation.FuncAnimation(fig, animate, init_func=init, frames=72, blit=True, interval=20, repeat=True)

    if save == 1:
        f = r"C:\Users\nicol\OneDrive\Documents\Doctorat\Regen_by_two\Coding\Images\scroll_SANDEN.gif" 
        writergif = animation.PillowWriter(fps=30) 
        ani.save(f, writer=writergif)
    else:
        pass
    
    
    return ani 
    

    

    
def scrollPlot(geo,theta):
    """
    

    Functions plotting the scroll set for a given orbiting angle.
   

    """
    
    
    geo.getCoordFullGeometry(theta)
    
    plt.figure()
    ax = plt.axes(xlim = (-(geo.D_wall/2+1),geo.D_wall/2+1),ylim = (-(geo.D_wall/2+1),geo.D_wall/2+1))
    
    (xport, yport) = geo.discharge_port()
    plt.plot(xport, yport,'k', zorder=1)
    plt.plot(geo.x_fscroll,geo.y_fscroll,'k',geo.x_oscroll, geo.y_oscroll, 'k')
    plt.fill(geo.x_fscroll, geo.y_fscroll, color='orange')
    plt.fill(geo.x_oscroll, geo.y_oscroll, color='grey', zorder=2)
    plt.plot(geo.x_wall,geo.y_wall,'k')

    # Plot forces 

    # scale = 10
    # f_s1 = volume_force_s1(geo,theta)["f_vect"]
    # f_s2 = volume_force_s2(geo,theta)["f_vect"]
    # f_sa = volume_force_sa(geo,theta)["f_vect"]
    
    # c1 = volume_force_c1(geo,theta)
    # c2 = volume_force_c2(geo,theta)
    # f_d1 = volume_force_d1(geo,theta)["f_vect"]
    # f_d2 = volume_force_d2(geo,theta)["f_vect"]
    # f_dd = volume_force_dd(geo,theta)["f_vect"]    
    
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
    # xy_polyS1 = volume_force_s1(geo,theta)["xy_poly"]
    # xy_polyS2 = volume_force_s2(geo,theta)["xy_poly"]
    # xy_polySA1 = volume_force_SA(geo,theta)["xy_polySA1"]
    # xy_polySA2 = volume_force_SA(geo,theta)["xy_polySA2"]
    
    # xy_polyC11 = volume_force_c1(geo,theta,geo.N_c_max)["xy_poly"]
    # xy_polyC12 = volume_force_c1(geo,theta,1)["xy_poly"]
    # xy_polyC21 = volume_force_c2(geo,theta,geo.N_c_max)["xy_poly"]
    # xy_polyC22 = volume_force_c2(geo,theta,1)["xy_poly"]
    
    #xy_polyD1 = volume_force_d1(geo,theta)["xy_poly"]
    # xy_polyD2 = volume_force_d2(geo,theta)["xy_poly"]
    # xy_polyDD = volume_force_dd(geo,theta)["xy_poly"]

    
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
    




def plot_area_port(geo):
    
    """
    

    Functions plotting the areas of the flow pathes
   

    """

    
    theta = np.linspace(0,2*pi,90)
    
    area_port =  np.zeros(len(theta))
    area_suction =  np.zeros(len(theta))
    area_discharge = np.zeros(len(theta))
    for i in range(len(theta)):
        (_,_,area_port[i]) = geoF.polyShadedPortArea(geo, theta[i])
        (_,_,area_suction[i]) = geoF.suctionFlowArea(geo, theta[i])
        (_,_,area_discharge[i]) = geoF.area_d_dd(geo,theta[i])
  
    plt.figure(figsize=(5.5,4.5))
    plt.xlim(theta.min(), theta.max())
    plt.xticks([0, pi/2, pi, 3*pi/2, 2*pi],
       [r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])

    plt.ylim(area_port.min()/1.4,area_port.max()*1.1)
    #plt.yticks([-1, 0, +1],
    #   [r'$-1$', r'$0$', r'$+1$'])
    plt.grid()
    plt.plot(theta, area_port, linewidth = 2)
    plt.xlabel('Orbiting angle [rad]',fontsize=16,  fontname="Times New Roman")
    plt.ylabel('Discharge Port Area [mm$^2$]',fontsize=16, fontname="Times New Roman")
    
    #################################
    plt.figure(figsize=(5.5,4.5))
    plt.xlim(theta.min(), theta.max())
    plt.xticks([0, pi/2, pi, 3*pi/2, 2*pi],
       [r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])

    plt.ylim(area_suction.min()/1.4,area_suction.max()*1.1)
    #plt.yticks([-1, 0, +1],
    #   [r'$-1$', r'$0$', r'$+1$'])
    plt.grid()
    plt.plot(theta, area_suction, linewidth = 2)
    plt.xlabel('Orbiting angle [rad]',fontsize=16,  fontname="Times New Roman")
    plt.ylabel('Suction Area [mm$^2$]',fontsize=16, fontname="Times New Roman")
    
    ##########################
    plt.figure(figsize=(5.5,4.5))
    plt.xlim(theta.min(), theta.max())
    plt.xticks([0, pi/2, pi, 3*pi/2, 2*pi],
       [r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])

    #plt.ylim(area_discharge.min()/,area_discharge.max()*1.1)
    #plt.yticks([-1, 0, +1],
    #   [r'$-1$', r'$0$', r'$+1$'])
    plt.grid()
    plt.plot(theta, area_discharge, linewidth = 2)
    plt.xlabel('Orbiting angle [rad]',fontsize=16,  fontname="Times New Roman")
    plt.ylabel('D-DD Area [mm$^2$]',fontsize=16, fontname="Times New Roman")

def plot_V_dV(geo):
     
     theta = np.linspace(0,2*pi,90)
     
     V_s1_vect = np.zeros(len(theta))
     V_s2_vect = np.zeros(len(theta))
     V_c1_vect = np.zeros((len(theta),geo.N_c_max))
     V_c2_vect = np.zeros((len(theta),geo.N_c_max))
     V_d1_vect = np.zeros(len(theta))
     V_d2_vect = np.zeros(len(theta))
     V_dd_vect = np.zeros(len(theta))
     V_dd_vect = np.zeros(len(theta))
     V_ddd_vect = np.zeros(len(theta))
     
     dV_s1_vect = np.zeros(len(theta))
     dV_s2_vect = np.zeros(len(theta))
     dV_c1_vect = np.zeros((len(theta),geo.N_c_max))
     dV_c2_vect = np.zeros((len(theta),geo.N_c_max))
     dV_d1_vect = np.zeros(len(theta))
     dV_d2_vect = np.zeros(len(theta))
     dV_dd_vect = np.zeros(len(theta))
     dV_dd_vect = np.zeros(len(theta))
     dV_ddd_vect = np.zeros(len(theta))
  
     for i in range(len(theta)): 
         V_s1_vect[i] = geoF.volume_force_s1(geo,theta[i])["V"]
         V_s2_vect[i] = geoF.volume_force_s2(geo,theta[i])["V"]
         
         for alpha in range(1,geo.N_c_max + 1,1):
             V_c1_vect[i,alpha-1] = geoF.volume_force_c1(geo,theta[i])[alpha]["V"]
             V_c2_vect[i,alpha-1] = geoF.volume_force_c2(geo,theta[i])[alpha]["V"]
         V_d1_vect[i] = geoF.volume_force_d1(geo,theta[i])["V"]
         V_d2_vect[i] = geoF.volume_force_d2(geo,theta[i])["V"]
         V_dd_vect[i] = geoF.volume_force_dd(geo,theta[i])["V"]
         
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
     
     
     #Volume plot
     plt.figure(figsize=(6,4.5),constrained_layout=True)
     plt.xlim(theta.min(), theta.max())
     plt.xticks([0, pi/2, pi, 3*pi/2, 2*pi],
       [r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])

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
     
     plt.xlabel('Orbiting angle [rad]',fontsize=16,  fontname="Times New Roman")
     plt.ylabel('Volume [cm$^3$]',fontsize=16, fontname="Times New Roman")   
         
      #Derivative plot
     plt.figure(figsize=(6,4.5))
     plt.xlim(theta.min(), theta.max())
     plt.xticks([0, pi/2, pi, 3*pi/2, 2*pi],
       [r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])

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

     plt.xlabel('Orbiting angle [rad]',fontsize=16,  fontname="Times New Roman")
     plt.ylabel('Volume Derivative [cm$^3$/rad]',fontsize=16, fontname="Times New Roman")       
         
     
def plot_evolution(geo, results, P_out, P_in):
    
     # Recorded results related to chambers are pulled out of the 'results' dictionary
     CV_s1  = results['CV']['s1']
     CV_c1  = results['CV']['c1']
     CV_d1  = results['CV']['d1']
     CV_dd  = results['CV']['dd']
     CV_ddd  = results['CV']['ddd']
     
     # Data are converted in a way that make this function work
     (vect_s1, vect_c1, vect_d1, vect_dd, vect_ddd) = solverF.convert_recorded_values(geo, CV_s1, CV_c1, CV_d1, CV_dd, CV_ddd)
     
     xdim = 5.5
     ydim = 4.5
     
     # PRESSURE
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
     
     # TEMPERATURE
     
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
           
     # ENTROPY
     
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
     
     # VAPOR QUALITY
     
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

     # MASS
     
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

     # OIL FLOW RATE
     
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

    # DENSITY
     
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

     # SUCTION MASS FLOW 
     
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

     # DISCHARGE MASS FLOW
     
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
     
     # MASS FLOW 'D-DD'
     
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
     
     # FORCES
     
     plt.figure(dpi=100,figsize=(xdim,ydim),constrained_layout=True)
     plt.xlim(vect_s1['theta'].min(), vect_s1['theta'].max())
     plt.xticks([0, pi/2, pi, 3*pi/2, 2*pi],
       [r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])

     max_limit = np.nanmax([np.nanmax(results['radial_force']),np.nanmax(results['tangential_force']),np.nanmax(results['axial_force'])]) 
     min_limit = np.nanmin([np.nanmin(results['radial_force']),np.nanmin(results['tangential_force']),np.nanmin(results['axial_force'])])      
     plt.ylim(min_limit*1.1, max_limit*1.1)

     plt.grid()
     
     plt.plot(vect_s1['theta'], results['radial_force'], linewidth = 2, label = 'radial force' )
     plt.plot(vect_s1['theta'], results['tangential_force'], linewidth = 2, label = 'tangential force' )
     plt.plot(vect_s1['theta'], results['axial_force'], linewidth = 2, label = 'axial force' )
     
     plt.xlabel('Orbiting angle [rad]', fontsize = 16,  fontname = "Times New Roman")
     plt.ylabel('Force [N]', fontsize = 16, fontname = "Times New Roman")
     plt.legend(fontsize=16, loc='best')
     
     # TILTING MOMENT
     
     plt.figure(dpi=100,figsize=(xdim,ydim),constrained_layout=True)
     plt.xlim(vect_s1['theta'].min(), vect_s1['theta'].max())
     plt.xticks([0, pi/2, pi, 3*pi/2, 2*pi],
       [r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])

     plt.ylim(np.nanmin(results['tilting_moment'])*0.9, np.nanmax(results['tilting_moment'])*1.1)

     plt.grid()
     
     plt.plot(vect_s1['theta'], results['tilting_moment'], linewidth = 2, label = 'tilting moment' )

     plt.xlabel('Orbiting angle [rad]', fontsize=16,  fontname="Times New Roman")
     plt.ylabel('Moment [Nm]', fontsize = 16, fontname="Times New Roman")
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
    
    
def PV_diagramm(geo, CV, P_out, P_in):  
    
     CV_s1   = CV['s1']
     CV_c1   = CV['c1']
     CV_d1   = CV['d1']
     CV_dd   = CV['dd']
     CV_ddd  = CV['ddd']
     
     (vect_s1, vect_c1, vect_d1, vect_dd, vect_ddd) = solverF.convert_recorded_values(geo, CV_s1, CV_c1, CV_d1, CV_dd, CV_ddd)
     
     
     theta = np.linspace(0,2*pi,181)
     V_dd = np.zeros(len(theta))
     V_ddd = np.zeros(len(theta))
     for i in range(len(theta)):
         V_dd[i] = geoF.volume_force_dd(geo,theta[i],False)["V"]/1e9
         V_ddd[i] = geoF.volume_force_dd(geo,theta[i],False)["V"]/1e9
     
     coefs_V_dd = np.array(np.polyfit(theta, V_dd, 5)) 
     V_dd_vect = np.polyval(coefs_V_dd, CV_s1.theta_vect_final)
     
     # plt.figure(figsize=(7,4.5),constrained_layout=True)
     # plt.plot(theta,V_dd)
     # plt.plot(CV_s1.theta_vect_final,V_dd_vect)
     
     (V_s, P_s) = (2*vect_s1['V'], vect_s1['P'])
     
     V_c = 2*vect_c1[1]['V']
     P_c = vect_c1[1]['P']
     for alpha in range(2,geo.N_c_max):
          (V_c, P_c) = (np.concatenate((V_c, 2*vect_c1[alpha]['V']), axis=None), np.concatenate((P_c, vect_c1[alpha]['P']), axis=None))


     V_ddd_vect_1 =  vect_ddd['V'].copy() - V_dd_vect
     V_ddd_vect_2 =  vect_ddd['V'].copy() - V_dd_vect
     P_ddd_vect_1 = vect_ddd['P'].copy()
     P_ddd_vect_2 = vect_ddd['P'].copy()
     
     (V_c_last, P_c_last) = (2*np.nan_to_num(vect_c1[geo.N_c_max]['V']), np.nan_to_num(vect_c1[geo.N_c_max]['P']))
     (V_d_last, P_d_last) = (2*np.nan_to_num(vect_d1['V']), np.nan_to_num(vect_d1['P']))
     V_dis = V_d_last[V_d_last != 0][0]
     
     for i in range(len(CV_s1.theta_vect_final)):
         if  CV_s1.theta_vect_final[i] <= geo.theta_d:
             V_ddd_vect_1[i] = 0
             P_ddd_vect_1[i] = 0
     V_ddd_vect_1 = np.nan_to_num(V_ddd_vect_1)
     P_ddd_vect_1 = np.nan_to_num(P_ddd_vect_1)
     (V_ddd_1, P_ddd_1) = (np.nan_to_num(CV_ddd.V_vect_final), np.nan_to_num(CV_ddd.P_vect_final))
     
     (V_dis1, P_dis1) = (V_c_last + V_d_last + V_ddd_vect_1, P_c_last + P_d_last + P_ddd_vect_1)
     
     for i in range(len(CV_s1.theta_vect_final)):
         if  CV_s1.theta_vect_final[i] > geo.theta_d:
             V_ddd_vect_2[i] = np.nan
             P_ddd_vect_2[i] = np.nan
     
     V_vect = np.concatenate((V_s, V_c , V_dis1, V_ddd_vect_2), axis=None)*1e6
     P_vect = np.concatenate((P_s, P_c, P_dis1, P_ddd_vect_2), axis=None)/1e5
     

     plt.figure(figsize=(4,3),constrained_layout=True)
     plt.xlim(0, np.nanmax(V_vect)*1.1)
     
     plt.ylim(P_in*0.9*1e-5, P_out*1.1*1e-5)
     # plt.ylim(np.nanmin(P_vect)/1.1, np.nanmax(P_vect)*1.1)
     #plt.yticks([-1, 0, +1],
     #   [r'$-1$', r'$0$', r'$+1$'])
     plt.grid()
     
     
     V_plot = np.linspace(0,np.nanmax(V_vect)*1.1, 10)
     P_plot = np.linspace(1,9, 10)
     plt.plot(V_vect, P_vect, color = u'#ff7f0e', linewidth = 2 )
     plt.plot(  V_plot  ,  P_out/1e5*np.ones(len( V_plot)), 'k--' , label = r'$P_{out}$')
     plt.plot(  V_plot ,  P_in/1e5*np.ones(len( V_plot)), 'k--' , label = r'$P_{in}$')
     
     plt.plot(  V_dis*1e6*np.ones(len( P_plot)) ,  P_plot, 'k-.' , label = r'$P_{in}$')
     
     plt.xlabel('Volume [cm$^3$]', fontsize=16,  fontname="Times New Roman")
     plt.ylabel('Pressure [bar]', fontsize=16, fontname="Times New Roman")       
     
     V_vect_plus = np.concatenate((V_s, V_c , V_dis1, V_ddd_vect_2, V_s[0]), axis=None)*1e6
     P_vect_plus = np.concatenate((P_s, P_c, P_dis1, P_ddd_vect_2, P_s[0] ), axis=None)/1e5
      
     W = intg.trapz(np.nan_to_num(V_vect)/1e6, np.nan_to_num(P_vect)*1e5 )
     
     # W_new = polygonArea(np.nan_to_num(V_vect_plus)/1e6, np.nan_to_num(P_vect_plus)*1e5)

     return W
 

 
    
def plot_volumes(geo, CV):
    
   CV_s1   = CV['s1']
   CV_c1   = CV['c1']
   CV_d1   = CV['d1']
   CV_dd   = CV['dd']
   CV_ddd  = CV['ddd']
    
   (vect_s1, vect_c1, vect_d1, vect_dd, vect_ddd) = geoF.convert_recorded_values(geo, CV_s1, CV_c1, CV_d1, CV_dd, CV_ddd)
        
    
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



   
                
             