#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 09:24:10 2020

@author: ali
"""
import numpy as np
import pylab as pl

np.random.seed(100)

class Ball: # task 2
    """
    Class for creating the balls and  the container in the simulation.
    For the container, radius must be negative.
    
    ---------------------------------------------------------------------------
    Methods:
        - a plethora of access methods for variables hidden in the initialisation
        - ke: calculates the kinetic energy of a ball
        - momt: calculates the momentum of a ball
        - move: moves a ball to a given time
        - time_to_collision: calculates the time to the collision of two balls
        - collide: collides two balls and updates their velocities 
        - get_patch: creates and returns a patch for a ball to be 
        displayed in the simulation animation
    """
    def __init__(self, m = 1.0, rad = 0.1, \
                 pos = np.array([0.0, 0.0], 'float64'), \
                 vel = np.array([0.0, 0.0], 'float64'),
                 ball_id = -1): # in SI units
        """
        Parameters (all in SI units where appropriate)
        ----------
        m : float, optional
            mass of ball. The default is 1.0 
        rad : float, optional
            radius of ball. The default is 0.1
        pos : numpy array, optional
            2D position of ball. The default is np.array([0.0, 0.0])
        vel : numpy array, optional
            2D velocity of ball. The default is np.array([0.0, 0.0])
        ball_id: int, optional
            unique ball identifier. The default is -1 (should only be used for
                                                  container)

        Returns
        -------
        None.
        
        Raises
        ------
        Exception
            pos only takes 2D vectors
        Exception
            vel only takes 2D vectors
        """
        
        # id is there for debugging purposes
        # default is -1 for a container, balls should have id >= 0
        
        # initialise as np array automatically so vector operators always work
        self._ball_id = ball_id
        self._m = float(m)
        self._rad = float(rad)
        # allows user to input a list
        self._pos = np.array(pos, 'float64')
        self._vel = np.array(vel, 'float64')
        # all args will automatically be made into floats
        
        if len(self._pos) > 2:
            raise Exception("pos only takes 2D vectors")
        if len(self._vel) > 2:
            raise Exception("vel only takes 2D vectors")
            
    def ball_id(self):
        return self._ball_id
    def pos(self):
        return self._pos
    def vel(self):
        return self._vel
    def m(self):
        return self._m
    def rad(self):
        return self._rad
    
    def ke(self):
        """ returns KE for a given ball"""
        return 0.5 * self._m * np.dot(self._vel, self._vel)
    def momt(self):
        """ returns vector momt for a given ball"""
        return self._m * self._vel
    
    def move(self, dt):
        """
        Moves ball forward in time by increment dt

        Parameters
        ----------
        dt : float
            time increment to move ball forward by

        Returns
        -------
        None.

        """
        # check that dt is a valid float to avoid an error in next_collision()
        # method in Simulation class
        if dt != None:
            self._pos += self._vel * dt # updating position
        
    def time_to_collision(self, other): # from task 1
        """
        Calculates time to collision with other ball by solving quadratic 
        equation for two elastically colliding spheres.
        
        Uses the discriminant and dot product of the relative positions and 
        relative velocities of the balls to discard imaginary solutions (balls 
        travelling parallel to each other) and negative solutions (balls collided
        in the past), respectively.
                                                                   
        Accounts for ball-ball collisions, as well as ball-container collisions.

        Parameters
        ----------
        other : Ball object 
            Other ball involved in collision

        Returns
        -------
        root_one or root_two: float
            time to collision
        None.
        
        
        Raises
        ------
        Exception
            A relative velocity of 0 will yield infinite dt as 
            they will never collide

        """
        #define relative variables 
        # round floats to 7dp to avoid float problems
        rel_pos = np.round(self._pos - other._pos, 7)
        rel_vel = np.round(self._vel - other._vel, 7)
        rel_rad = np.round(self._rad + other._rad, 7)
        # if other is container, it will be a subtraction
        # as radius of container is -ve and rel_rad will be < 0
        
        # avoid dividing by 0 when finding roots
        comparison = rel_vel == np.array([0,0])
        if comparison.all(): # returns True if all elements evaluated as True
            raise Exception(""""A relative velocity of 0 will yield infinite dt as 
they will never collide""")
        
        # define values to find roots here to make more readable 
        a = np.dot(rel_vel, rel_vel)
        b = 2 * np.dot(rel_pos, rel_vel)
        c = np.dot(rel_pos, rel_pos) - rel_rad ** 2


        discrim = b ** 2 - 4 * a * c
            
        
        # discarding imaginary sols (balls parallel)
        if discrim < 0:
            return None
        # ball-conatiner collisions
        elif rel_rad < 0 :
            root_one = (-b + np.sqrt(discrim))\
                       / (2 * a) # max +ve sol
                       
            """ don't need to check if "balls" moving towards each other
            as if the ball is in the container, it will always have a future
            collision, so we know the max root in this case will be +ve"""
            
            return root_one
        # discarding balls that collided in the past
        else: 
            if np.dot(rel_pos, rel_vel) < 0:
                # dot product must be -ve if balls are moving towards each other
                root_two =(-b - np.sqrt(discrim))\
                       / (2 * a) # min sol
                """ root_two is the min +ve sol, as only calculated if balls are moving
                       towards each other """
                return root_two
            else:
                return None
            
        
    def collide(self, other):
        """
        Collides two Ball objects, updating their velocities

        Parameters
        ----------
        other : Ball object
            Other ball involvved in collision

        Returns
        -------
        None.

        """
             
        # makes changes to velocities of both colliding balls
        
        # define relative variables 
        rel_pos = self._pos - other._pos
        rel_vel = self._vel - other._vel
        
        # update velocities 
        self._vel -= ((2 * other._m) * np.dot(rel_pos, rel_vel) * rel_pos) / \
        (np.dot(rel_pos, rel_pos) * (self._m + other._m))
        
        other._vel -= ((2 * self._m) * np.dot(rel_pos, rel_vel) * (-1 * rel_pos)) / \
        (np.dot(rel_pos, rel_pos) * (self._m + other._m))
        
        
        
    def get_patch(self):
        """
        Generates a patch for a ball to be visible in an animation.

        Returns
        -------
        patch : matplotlib patch object
            Circle denoting the ball. No fill for container Ball objects

        """
        # method for getting the circle primitive of the instance
        # used for animation
        
        # give different colours for debugging purposes 
        # does not indicate different species of particle
        noise = np.random.random(size=3) * 0.3
        color = (0.5 + noise[0], 0.5 + noise[1], 0.5 + noise[2])
        
        patch = pl.Circle(self._pos, self._rad, fc=color, fill=(self._rad > 0))
        # only fills if instance is a ball, container is outline                  
        return patch
        
        
        
    def __repr__(self):
        """
        Representation of a given Ball object

        Returns
        -------
        string
            Description of Ball object that can be used to reproduce the object

        """
        return "%s(m = %g, rad = %g, pos = %s, vel = %s, id = %s)" \
        % ("Ball",self._m, self._rad, self._pos, self._vel, self._ball_id)
    
    def __str__(self):
        """
        User-friendly description of Ball object

        Returns
        -------
        string
            Description of Ball object for user.
        
        """
        return "Ball id %s with mass %g kg, radius %g m at position %s m with velocity %s m/s" \
        % (self._ball_id, self._m, self._rad, self._pos, self._vel)
        
    
        
        
        
