#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 15:46:25 2020

@author: ali
"""
# this script was testing the simulation class prior to adjusting for many balls
import simulation_singleball as sim
import balls as ba
import numpy as np




#%% task 4: testing next_collision

ball = ba.Ball(rad = 1.0, pos = [-5.0, 0.0], vel = [1.0, 0.0])

# give container effectively infinite mass and needs negative radius
container = ba.Ball(m = 1e100, rad = -10) 


simul = sim.Simulation(ball, container)

print(ball.time_to_collision(container)) # should be 14 s
# 14.0 as expected


simul.next_collision()

print(ball.pos()) # should be container radius - ball radius = 9 for x pos
# returns [9. 0.] as expected

#%% test __repr__ and __str__

print(repr(simul))
print(simul) # gives __str__

# works as expected

#%% check velocity after collision

print(ball.vel()) 

# gives [-1.  0.] which is what we expect as elastic collision 

print(container.vel())

# gives order e-100 so effectly 0

print(ball.time_to_collision(container))

# gives 18 > 0 which is what we expect as the centre of the ball travels 18 m

#%% perform next collision and check again

simul.next_collision()

print(ball.vel()) 

# gives [1.  0.] which is what we expect 

print(ball.time_to_collision(container))

# gives 18 s again as going to ither side of container 

#%% test for ball with 2d vel (task 4)

a = ba.Ball(rad = 1.0, pos = [-5.0, 0.0], vel = [1.0, 1.0])

simuly = sim.Simulation(a, container)

simuly.next_collision()

#%% check velocity after collision

print(a.vel()) 

# gives [ 0.03115432 -1.41387036] which is valid

print(container.vel())

# gives order 1e-100 so effectly 0

print(a.time_to_collision(container))

# gives 11.704699910719622 > 0 which is valid

#%% Task 5 animation

# first test for stationary ball
testball = ba.Ball(pos = [0,0], rad = 1)

testsimul = sim.Simulation(testball, container)

#testsimul.run(10) # commented out to avoud raising the error

#shows ball for first frame and then ball disappears as
#there is no next collision and we get a divide by 0 error. Extension would be to modify
#the run method to stay on the same frame if the ball is 
#stationary 

#%% test for ball with only x vel

simul.run(10) # runs as expected

#%% test for ball with x and y vel

simuly.run(10) # looks plausible 



#%% Task 6 test CoE and conv of momt as well as pressure calc

print(simul.coe(ball, container)) # returns True so KE conserved 

print(simul.com(ball, container)) # returnd False
# this is because container has v high mass and v low vel
# giving momt of order 1, this will be fixed for many balls

# pressure was attempted (commented out in simulation_singleball),
# however I decided to move on to the many balls simulation and implement it
# properly there due to time constraintd 
# if I was to implement it, I would have it append time to collision
# and momt imparted to container for each collision in run and then use 
# that in the pressure method
