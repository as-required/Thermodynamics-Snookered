#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 15:31:09 2020

@author: ali
"""
import numpy as np
import pylab as pl


class Simulation:
    """ 
    Class to create simulation instances for one ball in a circular container
    
    Provides methods for getting the constituent Ball objects
    and for advancing the simulation to the next collision
    """
    
    def __init__(self, ball, container): # task 3
        # ball and container should both be instances of Ball
        self._ball = ball
        self._container = container
        
    def ball(self):
        return self._ball
    
    def container(self):
        return self._conatiner
        
    def next_collision(self):
        # find time to next collison
        dt = self._ball.time_to_collision(self._container)
        # move system to time of next collision
        self._ball.move(dt)
        
        # perform the collision
        self._ball.collide(self._container)
        
    def __repr__(self):
            return """%s(ball = %s,
container = %s)""" \
        % ("Simulation", repr(self._ball), repr(self._container))
        
    def __str__(self):
            return """Simulation of ball with mass %g kg, radius %g m at position %s m
with velcoity %s m/s in a container of radius %g m""" \
        % (self._ball.m(), self._ball.rad(), self._ball.pos(), self._ball.vel(),\
        -self._container.rad())
            #use access methods here as now accesssing hidden attributes from another class (Ball)
        

    def run(self, num_frames, animate = True): # task 5
        if animate:
            f = pl.figure()
            ax = pl.axes(xlim=(-10, 10), ylim=(-10, 10))
            ax.add_artist(self._container.get_patch())
            ax.add_patch(self._ball.get_patch())
            
            for frame in range(num_frames):
                # the (n-1)th frame is the nth collision
                # each frame is a collision
                self.next_collision()
                
                # showing in for loop to display a continuous animation
                # in the same window
                pl.pause(0.001)
                pl.show()
            
    
    ### task 6
    def coe(self, *balls):
        # tests for conservation of energy 
        # used *balls to pass in an arbitrary number of balls and the container to be general
        # kinetic only in this case as ideal gas
        ke_i = 0
        
        for b in balls:
            ke_i += 0.5 * b.m() * np.dot(b.vel(), b.vel())
        
        self.next_collision() # perform next collisoon to calculate final KE
        
        ke_f = 0
        
        for b in balls:
            ke_f += 0.5 * b.m() * np.dot(b.vel(), b.vel())
        # same calc as for ke_i, but now ball._vel is updated and ke_f is stored
        # in different memory so they are ke_f and ke_i are not the same
        # returns True if KE is conserved (elastic collision)
        # we must account for floats being slightly different
        return abs(ke_f - ke_i) < np.finfo(float).eps 
    
    def com(self, *balls):
         # tests for conservation of energy
         # used *balls to pass in an arbitrary number of balls and the container to be general
         p_i = 0
         
         for b in balls:
             p_i += b.m() * np.linalg.norm(b.vel())
             # print(np.linalg.norm(b.vel())) for checking vels
        
         self.next_collision()
         p_f = 0
         
         for b in balls:
             p_f += b.m() * np.linalg.norm(b.vel())
             # print(np.linalg.norm(b.vel())) for checking vels
         # same calc as for p_i, but now ball._vel is updated and p_f is stored
         # in different memory so they are p_f and p_i are not the same
         # returns True if momentum is conserved
         # we must account for floats being slightly different
            
         return abs(p_f - p_i) < np.finfo(float).eps
     
    """def pressure(self, *balls, container, time):
        # calculates (avg) pressure ON container 
        # total change of motm/ area of container/ time simulation has run for
        # should be const after long enough time 
        # can pass in as many balls as needed because of *balls
        
        P = 0
        A = np.pi * self._container.rad() ** 2
        # total time 
        t = 
        
        for b in balls:
            rel_rad = b.rad() + container.rad()
            # pressure only created from ball-container collisions
            if rel_rad < 0:
                P += (2 * b.m() * b.vel()) / (t * A) """
        
        
        
        
    
        
        
        
        
