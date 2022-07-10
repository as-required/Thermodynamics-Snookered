#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 10:42:44 2020

@author: ali
"""

""" tests two ball next collision. Trying to debug n calculations 
next_collision method, but was abandoned and I just created an n^2 
calculation method due to time """

import testsimulation as sim
import balls as ba
import numpy as np
#%% test 



sim = sim.Simulation(nballs = 2, mass =1, radius = 2, container_radius = 15)



#%% run simulation 

sim.run(100, animate = True)

# the two balls are still not colliding, says their collision is in the past
# so time_to_collision must be wrong

# I found the error in time_to_collision and fixed it
# I was checking for np.dot(rel_pos, rel_vel) > 0 
# when it should have been < 0
# However I then had a logic error in next_collision, which I got close 
# to fixing however, I ran out of time