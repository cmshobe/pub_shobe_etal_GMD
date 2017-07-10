# -*- coding: utf-8 -*-
"""
*This script is a driver for the SPACE Landlab component.

*It accompanies Shobe et al., manuscript submitted to Geosci. Mod. Dev.

*This script replicates the test case shown in section 5.3 of Shobe et al,
*in which the SPACE model is compared to analytical solutions for channel
*slope, sediment thickness, and sediment flux in a case where both sediment 
entrainment and bedrock erosion are happening at once.

*Derivations of analytical solutions used here for comparison may be found
*in Shobe et al., manuscript submitted to Geosci. Mod. Dev.

@author: Charles M. Shobe, University of Colorado, charles.shobe@colorado.edu
"""
##Import necessary libraries and components#################################
import numpy as np
import matplotlib.pyplot as plt
from landlab import RasterModelGrid
from landlab.components.flow_routing import FlowRouter
from landlab.components import DepressionFinderAndRouter
from landlab.components import Space

#Instantiate model grid and add initial condition###########################
#Instantiate model grid
num_rows = 20
num_cols = 20
node_spacing = 100.0 #m
mg = RasterModelGrid((num_rows, num_cols), node_spacing)

#Create initial topography
np.random.seed(seed = 5000) #constant seed for constant random roughness
mg.add_zeros('node', 'topographic__elevation')
mg['node']['topographic__elevation'] += mg.node_y / 100000 + mg.node_x / \
    100000 + np.random.rand(len(mg.node_y)) / 10000
    
#Create a grid field for soil depth
mg.add_zeros('node', 'soil__depth')

#Give an initial soil depth, in this case 0 m
mg.at_node['soil__depth'][:] = 0.0 #meters

#Create a field for bedrock elevation
mg.add_zeros('node', 'bedrock__elevation')

#Make bedrock elevation equal to topographic elevation
mg.at_node['bedrock__elevation'][:] = mg.at_node['topographic__elevation']

#Increase topographic elevation by soil depth, in case there is any soil
mg.at_node['topographic__elevation'][:] += mg.at_node['soil__depth']  
    
    
##SET GRID BOUNDARY CONDITIONS################################################
#Close all domain borders
mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True,
                                       left_is_closed=True,
                                       right_is_closed=True,
                                       top_is_closed=True)
 
#Open lower-left (southwest) corner                                      
mg.set_watershed_boundary_condition_outlet_id(0,
    mg['node']['topographic__elevation'], -9999.)
    
##Instantiate components###################################################
#Instantiate flow router
fr = FlowRouter(mg) 

#Instantiate depression finder and router; optional
df = DepressionFinderAndRouter(mg)

#Define parameters
K_sed = 0.01
K_br = 0.005
F_f = 0.
phi = 0.
H_star = 1.
v_s = 5.0
m_sp = 0.5
n_sp = 1.0
sp_crit_sed = 0.
sp_crit_br = 0.

#Instantiate SPACE model with chosen parameters
sp = Space(mg, K_sed=K_sed, K_br=K_br, 
                    F_f=F_f, phi=phi, H_star=H_star, v_s=v_s, m_sp=m_sp, 
                    n_sp=n_sp, sp_crit_sed=sp_crit_sed,
                    sp_crit_br=sp_crit_br, method='simple_stream_power')

##Run time loop###########################################################
#Set model timestep
timestep = 1.0 #years

#Set elapsed time to zero
elapsed_time = 0

#Set model run time
run_time = 100000 #years

#Set rock uplift rate
rock_uplift_rate=1e-4 #m/yr


while elapsed_time < run_time: 
    #Run the flow router
    fr.run_one_step()
    
    #Run the depression finder and router; optional
    df.map_depressions()
    
    #Get list of nodes in depressions; only 
    #used if using DepressionFinderAndRouter
    flooded = np.where(df.flood_status==3)[0]
    
    #Run the SPACE model for one timestep
    sp.run_one_step(dt = timestep, flooded_nodes=flooded)
    
    #Move bedrock elevation of core nodes upwards relative to baselevel
    #at the rock uplift rate
    mg.at_node['bedrock__elevation'][mg.core_nodes] += rock_uplift_rate * timestep
    
    #Strip any soil from basin outlet so that all topographic change is due to 
    #rock uplift
    mg.at_node['soil__depth'][0] = 0. 

    #Recalculate topographic elevation to account for rock uplift
    mg.at_node['topographic__elevation'][:] = \
        mg.at_node['bedrock__elevation'][:] + mg.at_node['soil__depth'][:]
    
    if elapsed_time % 1000 == 0:
        print elapsed_time
    elapsed_time += timestep

##Save data for testing against analytical solutions#########################

#Save the soil__depth field
sed_depth = np.zeros(len(mg.at_node['topographic__elevation'][mg.core_nodes]))
sed_depth[:] = mg.at_node['soil__depth'][mg.core_nodes]
np.save('sed_depth.npy', sed_depth)

#Save the drainage_area field
area = np.zeros(len(mg.at_node['topographic__elevation'][mg.core_nodes]))
area[:] = mg.at_node['drainage_area'][mg.core_nodes]
np.save('area.npy', area)

#Save the slope field
slope = np.zeros(len(mg.at_node['topographic__elevation'][mg.core_nodes]))
slope[:] = mg.at_node['topographic__steepest_slope'][mg.core_nodes]
np.save('slope.npy', slope)

#save the sediment flux field
qs = np.zeros(len(mg.at_node['topographic__elevation'][mg.core_nodes]))
qs[:] = mg.at_node['sediment__flux'][mg.core_nodes]
np.save('qs.npy', qs)

##Compare slope, sediment depth, and sediment flux against predictions#####

#Define runoff parameter r, where Q=Ar
r = 1.0 #m/yr

#First, sediment flux. Steady state Qs at every point is expect to be Qs=U*A
#as long as F_f=0 and porosity = 0

#Instantiate figure
sedflux_fig = plt.figure()

#Instantiate subplot
sedflux_plot = plt.subplot()

#Plot drainage area on the x-axis and Qs on the y-axis.
sedflux_plot.scatter(area, qs, marker = 'o', color = 'k', label = 'Model Qs')

#set x-axis to log scale
sedflux_plot.set_xscale('log')

#Label axes
sedflux_plot.set_xlabel(r'Drainage area [m$^2$]')
sedflux_plot.set_ylabel(r'Sediment flux [m$^3$/s]')

#Plot analytical solution
analytical_domain = np.arange(min(area),max(area), 1)
sedflux_plot.plot(analytical_domain, rock_uplift_rate * analytical_domain, linewidth = 2, 
                  color = 'grey', linestyle = '-', label = 'Analytical Qs')

#Make a legend
sedflux_plot.legend(loc='upper left')

#Now compare slope-are relationship against analytical prediction

#Instantiate figure
slope_area_fig = plt.figure()

#Instantiate subplot
slope_area_plot = plt.subplot()

#Plot drainage area on the x-axis and slope on the y-axis.
slope_area_plot.scatter(area, slope, marker = 'o', color = 'k', 
                        label = 'Model slope')

#set x-axis to log scale
slope_area_plot.set_xscale('log')
slope_area_plot.set_yscale('log')

#Label axes
slope_area_plot.set_xlabel(r'Drainage area [m$^2$]')
slope_area_plot.set_ylabel('Slope [-]')

#Calculate and plot analytical solution
slope_analytical_soln = np.power((rock_uplift_rate * v_s) / \
    (K_sed * np.power(analytical_domain, m_sp) * r) + rock_uplift_rate / \
    (K_br * np.power(analytical_domain, m_sp)), 1 / n_sp)
    
slope_area_plot.plot(analytical_domain, slope_analytical_soln, linewidth = 2, 
                  color = 'grey', linestyle = '-', label = 'Analytical slope')

#Make a legend
slope_area_plot.legend(loc='upper right')
               
#Finally, check sediment depth against the analytical solution.
#Because sediment thickness should be constant, we expect a flat line

#Instantiate figure
sed_depth_fig = plt.figure()

#Instantiate subplot
sed_depth_plot = plt.subplot()

#Plot drainage area on the x-axis and slope on the y-axis.
sed_depth_plot.scatter(area, sed_depth, marker = 'o', color = 'k', 
                        label = 'Model depth')

#set x-axis to log scale
sed_depth_plot.set_xscale('log')

#Label axes
sed_depth_plot.set_xlabel(r'Drainage area [m$^2$]')
sed_depth_plot.set_ylabel('Sediment depth [m]')

#Calculate and plot analytical solution
#Calculate analytical sediment depth
sed_depth_analytical_soln = -H_star * np.log(1 - (v_s / ((K_sed * r / K_br) + v_s)))
sed_depth_analytical_array = np.repeat(sed_depth_analytical_soln, 
                                       len(analytical_domain))
    
sed_depth_plot.plot(analytical_domain, sed_depth_analytical_array, linewidth = 2, 
                  color = 'grey', linestyle = '-', label = 'Analytical depth')

#Make a legend
sed_depth_plot.legend(loc='upper right')

#Save all three plots
sedflux_fig.savefig('sediment_depth.eps')
slope_area_fig.savefig('slope_area.eps')
sed_depth_fig.savefig('sediment_depth.eps')