#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 15:35:46 2020

@author: ali
"""

import simulation as sim
import balls as ba
import numpy as np

# note: all numerical values retrieved where using 
# np.random.seed(100) enabled in simulation.py for reproducability
# this only ensures reproducability for fresh consoles (where the cell
# in question hasn't been run before)



#%% test animation to check if simulation is working properly
# also used commented out savefigs in run method to see what's happening
# frame by frame


sim1 = sim.Simulation(nballs = 5, mass = 1, radius = 0.5, container_radius = 10)


sim1.run(50, animate = True) 

# looks correct: balls not escaping container, jumping
# balls colliding with each other properly 
# very that balls are colliding with each other by looking in console
# and seeing frames that don't have id = -1 (id of the container)

# to verfiy that balls do not start overlapping
# uncomment savefig 's in run method and look at frame 0
# they are arranged to never overlap in frame 0 (or any other frame)

#%% test simulation time and momt imparted to container

print("Simulation has run for %.3g s" % (sim1.sim_time()))
print("A total of %.3g kgm/s of momentum has been imparted to the container" \
      % (sim1.momt_imparted()))
    
# values seem reasonable 

#%% task 9 histograms PRODUCES TASK 9 HISTOGRAMS
10
# run simulation with many frames and many balls to get accurate behaviour of 
# system over time
# should reset console after each time when plotting, to get system behavior from start

# run with many balls over a long time (many frames) to get good statistics
simh = sim.Simulation(nballs = 100, mass = 1, radius = 0.1, container_radius = 10)
# must run sim before plotting histograms
simh.run(300, animate = False)

simh.centre_distance_hist(15) 
# makes sense that skewed to the right as most balls will be at edges


simh.dist_between_balls_hist()
# makes sense that peaks in middle as most probable separation is the avg separation
#%% task 10

# momentum is also conserved, not just KE
# this will be shown in plots in task  11
# in 2D, 2/2 T = avg KE --> T = avg KE

# test temperature calculation 
sim10 = sim.Simulation(50)

# store simulation after being run as variable to avoid running twice
run_sim10 = sim10.run(100, animate = False, values = True)

# temp 
print(run_sim10[0]) # constant temp at circa 8e24 K (running on a fresh console)
# astronomical temp is because we are using m = 1kg

# pressure
print(run_sim10[1]) # pressure stabilises at the end to about 6 Pa (the time I ran it)
# initially no pressure as takes time for first container collision to happen
# as in this sim, many balls so start close to each other
#%%
# if velocities are doubled, temperature should quadruple 
# particle motion should also randomise faster 

sim101 = sim.Simulation(50, vfactor = 2)

sim101_run = sim101.run(100, values = True)

# temp 
print(sim101_run[0]) # constant temp circa 3.2e25
# so about quadruple as expected, 

# pressure should double

print(sim101_run[1]) # get it settling at circa 20 Pa (varies on how many times 
# it's been run)


#%% try sim101 with less frames to check how quickly system settles

sim1012 = sim.Simulation(50, vfactor = 2)

sim1012_run = sim1012.run(25, values = True)

# temp 
print(sim1012_run[0]) # constant temp circa 2.8e25
# so about quadruple as expected. Diff to last cel
# due to randomness in vels balls are initialised with

# pressure should double from 2 cells ago

print(sim1012_run[1]) # settles at circa 20 Pa 


#%% how animation would change
sim10.run(100, animate = True)
sim101.run(100, animate = True)

# there's a shorter time between collisions as all the balls are moving
# twice as fast so are more likely to collide with themselves/ container

#%% test pressure, time evolution

simp = sim.Simulation(100)

simp.run(300, pressure_graph= True ) # still climbing as many balls, so 
# not many collisions with container to reach equilibrium temp

#%% try with less balls ###PRODUCES PRESSURE VS TIME PLOT

simp2 = sim.Simulation(25)

simp2.run(500, pressure_graph = True) # eventually evens out

# sometimes you get an initial spike because a ball immediately hit the container
#%% task 11 PRODCUES CONNSERVATION PLOTS

# check consv of E

sime = sim.Simulation(100) # 100 balls

### RUN THIS TWICE TO GET PLOT SHOWN IN REPORT
# there seems to be a numerical error running it the first time
# that causes the axis scale for the com plot to interfere w title
sime.run(50, consv_graph = True) # horizontal line showing constant system KE and constant momt
# energy conserved 
# momt conserved

#uncertainy on gradient is 2.46e-13 J/s (coe)
#uncertainy on gradient is 1.14e-14 kgm/s^2 respectively, so v low (com)

# rest of tasks in physics_investigations




#%% test for large number of balls

simtest = sim.Simulation(nballs = 1000, mass = 1, radius = 0.1, container_radius = 50)

#simtest.run(100) # very, very slow
# in future, should try to reimplement order n calculation for next_collision
# had previously implemented this with dictionaries, but it was taking too long to debug
# (it would work for approx 10 frames and then balls would start escaping container)
# so rewrote an n^2 complexity algorithim due to time concerns

#%% try with realistic ball sizes (scale of atoms)

#gas = pressure_v_temp(vfactors, ball_radius = 1e-10, ball_mass = 1e-27 )

# THIS WILL NOT RUN
# reason why is ball radius is so small that the get_positions method
# tries to inscribe too many concentric circles in the container
# that the computer freezes
# it can't handle ball radii of order 0.001m, this could be improved in future


#%% test repr and str

simz = sim.Simulation(100)
print(repr(simz))
# Simulation(nballs = 100, mass = 1, radius = 0.5, container_radius = 10 
    #vfactor = 1, vel_sd = 10) # works
print(simz)
# Simulation of 100 balls with mass 1 kg, radius 0.5 m in a container of radius 10 m # works
