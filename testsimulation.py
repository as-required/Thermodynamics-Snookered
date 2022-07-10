#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 10:33:39 2020

@author: ali
"""

import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import balls as ba
from itertools import combinations 

# for debugging to get same simulation each time 
np.random.seed(100)

# define constants 

k = 1.38064852e-23

def freedman_diaconis(data, n):
        """ Returns optimised number of bins for histograms with n samples
        uses bin width = 2 IQR / cube root no. of samples"""
        
        bin_width = 2 * np.subtract(*np.percentile(data, [75, 25]))/ np.cbrt(n)
        no_bins = int((np.max(data) - np.min(data)) / bin_width)
        
        return no_bins
    
def maxwell(v, m, T):
    """ returns maxwell-boltzman distribution. v should be a numpy array, m is
    the mass of one ball and T is the temperature of the system"""
    probabilities = []
    for i in range(len(v)):
        expo = np.exp((-m * np.dot(v[i], v[i]))/(2 * k * T))
        # use v not v**2 as working in 2D
        coeff = (4 * np.pi * np.linalg.norm(v[i])) * ((m / (2 * np.pi * k * T)) ** (3/2))
        prob = coeff * expo
        probabilities.append(prob)
    
    return probabilities
    


class Simulation:
    """ 
    Class to create simulation instances for one ball in a circular container
    
    Provides methods for getting the constituent Ball objects
    and for advancing the simulation to the next collision
    """
    
    def __init__(self, nballs, mass = 1.0, radius = 0.5, container_radius = 10.0, temperature = 273.0):
        # container should be instance of Ball
        # nballs is number of balls to be initialised in simulation 
        # mass and radius refer to the balls
        # temperature in K
        self._nballs = nballs
        self._mass = float(mass)
        self._radius = float(radius) 
        # initialise container, so users don't have to worry about mass
        # also don't have to worry about making container radius -ve
        self._container_radius = float(container_radius) 
        self._container = ba.Ball(m = 1e100, rad = -self._container_radius)
        
        self._balls = []
        
        # arrange balls in randomly selected positions 
        positions = self.get_available_positions(self._container.rad(), radius, nballs)
        
        # initialising balls with position and velocity
        for i, pos in enumerate(positions):
            vel = np.random.normal(scale = 10, size = 2) # gaussian distt of vels w/ sd 10 m/s
            self._balls.append(ba.Ball(mass, radius, pos, vel, ball_id=i+1))
            
        # test for two balls to see why not colliding 
        #self._balls.append(ba.Ball(mass, radius, [0,0], [0,0], ball_id=10))
        #self._balls.append(ba.Ball(mass, radius, [-5,0], [1,0], ball_id=11))

        self._t_collision_per_pair = dict()
        # total time since start of simulation 
        self._time = 0
        # total momentum imparted to container since start of simulation 
        self._momt = 0
        
        
    def balls(self):
        return self._balls
    
    def container(self):
        return self._container
    
    def system_ke(self):
        """ Returns total KE for the system """
        total = self._container.ke()
        
        for b in self._balls:
            total += b.ke()
        
        return total
    
    
    
    
    def centre_distance_hist(self, bins = 10):
        """ plots a histogram of distance the balls extend from the centre"""
        dist = []
        for b in self._balls:
            dist.append(np.linalg.norm(b.pos()))
        
            
        fig, ax = plt.subplots()
        ax.set_title("Distance balls extend from centre")
        ax.set_xlabel("Distance from centre")
        ax.set_ylabel("Number of balls")
        ax.hist(dist, bins) #freedman - diaconis doesn't seem to work well here
        fig.savefig("distance_balls_from_centre")
        plt.show()
            
            
    def dist_between_balls_hist(self):
        
        """ plots a histogram of distance between pairs of balls"""
        
        pair_balls = combinations(self._balls, 2)
        list_pair_balls = list(pair_balls)
        dist = []
        
        for pair in list_pair_balls:
            dist.append(np.linalg.norm(pair[0].pos() - pair[1].pos()))
        
        
        fig, ax = plt.subplots()
        ax.set_title("Distance between pairs of balls")
        ax.set_xlabel("Distance")
        ax.set_ylabel("Number of pairs")
        ax.hist(dist, freedman_diaconis(dist, len(list_pair_balls)))
        fig.savefig("distance_between_balls")
        plt.show()
        
        
        
    
    def get_available_positions(self, r_container, r_ball, nballs):
         
        """ uses concentric circles in conatiner, with predetermined positions 
        for balls and which of these positions is filled is decided randomly"""
        
        
        r_container = abs(r_container)  # need +ve in subsequent calculations
        positions = []
        # ensures balls on diff circles don't overlap
        # container radius/ ball diameter gives balls just touching
        # so multiply by 0.75 to give space
        num_circles = int(np.floor(r_container/(2*r_ball)) * 0.75)
        
        #create concentric circles by incrementing their radii
        for i in range(num_circles):
            # first r will be 0, so one position will always exist in the centre
            r = r_container/num_circles * i
            
            # calc how many positions allowed on a given circle
            # divide circle circumference by ball diameter
            # gives theoretical max number of balls allowed on a given circle 
            # determined only by ball radius 
            # reduce this by 0.75 to increase spacing between balls
            # for r = 0 case, only one potential ball allowed
            if r == 0:
                num_balls = 1
            else:
                num_balls = int(np.floor(2 * np.pi * r/(2 * r_ball)) * 0.75)
            
            # place potential balls evenly around circle
            for n in range(num_balls):
                # get an angle that depends on number of balls
                theta = 2 * np.pi * n/num_balls
                # convert to cartesian form
                x_pos = r * np.cos(theta)
                y_pos = r * np.sin(theta)
                
                # store ball positions 
                positions.append((x_pos, y_pos))
                
        if nballs > len(positions):
            raise Exception(""""Too many balls of this size for this conatiner size, 
                  max no. of balls allowed for these parameters is %s """ % (len(positions)))
        
        # randomly select which positions get filled by the number of balls specified
        # replace is False to prevent same ball being placed twice
    
        # returns 1D numpy array with nballs no. of randomly selected integers
        # in range len(positions)                   
        random_indices = np.random.choice(len(positions), size = nballs, replace = False)
        
        positions_out = []
        # use the random indices to pick a random position and store in another list
        for random_index in random_indices:
            # access selected position from random index
            x_pos, y_pos = positions[random_index]
            positions_out.append((x_pos, y_pos))
            
        return positions_out

        
    def next_collision(self): #test with only 2 balls
        # list containing times to collision for each pair of balls, incl container
        objects = [*self._balls, self._container] # only take 2 balls
        
        # only possible pairs for 2 balls and container
        pairs = [(objects[0], objects[1]), (objects[0], objects[2]), (objects[1], objects[2])]
        
        for i in range(len(pairs)):
            ball1, ball2 = pairs[i][0], pairs[i][1]
            dt = ball1.time_to_collision(ball2)
            if dt is not None:
                ball1.move(dt)
                ball2.move(dt)
                print("Colliding %s with %s" % (ball1, ball2))
                ball1.collide(ball2)
            
        
        

        

        

    def run(self, num_frames, animate = False, graph = False):
        
        # physical observables 
        T = []
        P = []
        v = []
        t = []
        
        if animate:
            f = pl.figure()
            ax = pl.axes(xlim=(-self._container_radius, self._container_radius), ylim=(-self._container_radius, self._container_radius), aspect = 1)
            ax.add_artist(self._container.get_patch())
            for b in self._balls:
                ax.add_patch(b.get_patch())
        # saving figs to debug frame by frame        
        #pl.savefig('frame_0')
        for frame in range(num_frames):
            print("Frame: %s"%(frame+1))
            # the (n-1)th frame is the nth collision
            # each frame is a collision
            self.next_collision()
            
            # calc physical quantities and add to list to plot later
            #t.append(self._time)
            #T.append(2 * self.system_ke()/ 3 * k * len(self._balls)) # T = avg KE
            # used coeff from ideal gas eqn
            #P.append((self._momt/ self._time)/(np.pi * self._container_radius ** 2)) # avg pressure 
            
            #for b in self._balls:
                #v.append(np.linalg.norm(b.vel()))
            
            
            
            
            if animate: 
                pl.pause(0.0001)
                # saving figs to debug frame by frame 
                #pl.savefig('frame_{}'.format(frame+1))
                pl.show() # showing in for loop to display a continuous animation
            
            
    #def vel_dist(self, no_collisions):
        
        #for collision in no_collisions:
            
            
    
    def coe(self):
        # tests for conservation of energy 
        # kinetic only in this case as ideal gas
        ke_i = 0
        
        for b in self._balls:
            ke_i += 0.5 * b.m() * np.dot(b.vel(), b.vel())
        
        self.next_collision() # perform next collisoon to calculate final KE
        
        ke_f = 0
        
        for b in self._balls:
            ke_f += 0.5 * b.m() * np.dot(b.vel(), b.vel())
        # same calc as for ke_i, but now ball._vel is updated and ke_f is stored
        # in different memory so they are ke_f and ke_i are not the same
        # returns True if KE is conserved (elastic collision)
        # we must account for floats being slightly different
        return abs(ke_f - ke_i) < np.finfo(float).eps 
    
    def com(self):
         # tests for conservation of momentum
         p_i = 0
         
         for b in self._balls:
             p_i += b.m() * np.linalg.norm(b.vel())
             # print(np.linalg.norm(b.vel())) for checking vels
        
         self.next_collision()
         p_f = 0
         
         for b in self._balls:
             p_f += b.m() * np.linalg.norm(b.vel())
             # print(np.linalg.norm(b.vel())) for checking vels
         # same calc as for p_i, but now ball._vel is updated and p_f is stored
         # in different memory so they are p_f and p_i are not the same
         # returns True if momentum is conserved
         # we must account for floats being slightly different
            
         return abs(p_f - p_i) < np.finfo(float).eps
     
        
     

    def __repr__(self):
        return """%s(ball = %s, container = %s)""" % ("Simulation", repr(self._balls), repr(self._container))
        
    def __str__(self):
            return """Simulation of ball with mass %g kg, radius %g m at position %s m
with velcoity %s m/s in a container of radius %g m""" \
        % (self._balls.m(), self._balls.rad(), self._balls.pos(), self._balls.vel(),\
        -self._container.rad())
            #use access methods here as now accesssing hidden attributes from another class (Ball)
        