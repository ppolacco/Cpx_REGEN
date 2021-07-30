# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 11:46:14 2021

@author: nicol
"""

import solverFunctions as solverF
import matplotlib as plt
import diversePlots as divPlot

#%%
    
point = solverF.open_results('test_2')

geo = point['geo_inputs']['geo']
P_out = point['inputs']['P_out']
CV = point['CV'][0]


CV_s1   = CV['s1']
CV_c1   = CV['c1']
CV_d1   = CV['d1']
CV_dd   = CV['dd']
CV_ddd   = CV['ddd']



plt.close('all')

#plot_V_dV(geo)
# plot_area_port(geo)

divPlot.plot_evolution(geo, CV, point['inputs']['P_out'], point['inputs']['P_in'])
# basic_plot(np.array(point['inputs']['N']), np.array( point['eta_is']), 'N [RPM]', '$\eta_{is}$ [-]' )

# theta = 1.2
# scrollPlot(geo,theta)


#%%

divPlot.plot_volumes(geo, CV)
