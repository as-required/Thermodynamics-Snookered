#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 15:31:09 2020

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
        """
    Calculates optimised number of bins for histograms with n samples
    using the Freedman-Diaconis Rule:
        
        bin width = 2 IQR / cube root no. of samples.

    Parameters
    ----------
    data : list
        Histogram data.
    n : int
        Number of samples in data.

    Returns
    -------
    no_bins : int
        Optimised number of bins for the data

    """
        
        bin_width = 2 * np.subtract(*np.percentile(data, [75, 25]))/ np.cbrt(n)
        no_bins = int((np.max(data) - np.min(data)) / bin_width)
        
        return no_bins
    
    


class Simulation: # task 3
    """ 
    Class to create simulation instances for many balls in a circular container
    
    Provides methods for getting the constituent Ball objects
    and for advancing the simulation to the next collision, as well as keeping
    track of physical observables 
    
    ---------------------------------------------------------------------------
    Methods:
        
        - a plethora of access methods for variables hidden in the initialisation 
        - system_ke: calculates the total kinetic energy of the system
        - system_momt: calculates the total momentum of the system
        - centre_distance_hist: plots a histogram of the distance that the balls
        extend from the centre of the container
        - distance_between_balls_hist: plots a histogram of the inter-ball 
        separation 
        - get_available_positions: psuedo-randomly selects positions for balls on 
        concentric circles inscribed into the container. The algorithm ensures
        that it is impossible for balls to be initialised overlapping
        - next_collision: advances the simulation to the next collision and 
        performs it
        - run: runs the simulation for a given number of frames. Also provides
        optional animation, physical values and graphs
    """
    
    def __init__(self, nballs, mass = 1.0, radius = 0.5, container_radius = 10.0, \
                 vfactor = 1.0, vel_sd = 10.0):
        """
        Constructor for Simulation class.
        
        Uses get_available_positions to initialise balls with pseudo-random positions
        and uses numpy.random.normal to initialise balls with velocities pseudo-randomly
        drawn from a Gaussian distribution with default mean and standard deviation of
        0 m/s and 10 m/s respectively.
        
        Also keeps track of the time the simulation has gone on for in real time and
        the total momentum imparted by the balls to the container in that time
        

        Parameters (all in SI units where necessary)
        ----------
        nballs : int
            Number of balls for simulation
        mass : float, optional
            Mass of the balls. The default is 1.0.
        radius : float, optional
            Radius of the balls. The default is 0.5.
        container_radius : float, optional
            Radius of the container. The default is 10.0.
        vfactor : float, optional
            Multiplies all starting particle velocities by this factor. The default is 1.0.
        vel_sd : float, optional
            Sets the initial standard deviation of the particle velocities. The default is 10.

        Returns
        -------
        None.

        """
        # container should be instance of Ball
        # container radius should be +ve, as the init makes it -ve
        self._nballs = nballs
        self._mass = float(mass)
        self._radius = float(radius) 
        self._vfactor = float(vfactor)
        self._vel_sd = float(vel_sd)
        # initialise container, so users don't have to worry about mass
        # also don't have to worry about making container radius -ve
        self._container_radius = float(container_radius) 
        self._container = ba.Ball(m = 1e100, rad = -self._container_radius)
        
        self._balls = []
        
        # keeps track of how long sim has gone on for 
        self._simtime = 0 # s
        
        # keeps track of total momentum imparted to container
        # to calc pressure later
        
        self._momt_imparted = 0 # kgm/s
        
        
        # arrange balls in randomly selected positions 
        positions = self.get_available_positions(self._container.rad(), radius, nballs)
        
        # initialising balls with position and velocity (task 8)
        for i, pos in enumerate(positions):
            # gaussian dist of vels w/ mean 0 m/s and sd vel_sd m/s
            vel = vfactor * np.random.normal(scale = vel_sd, size = 2)
            self._balls.append(ba.Ball(mass, radius, pos, vel, ball_id = i+1))
            # ball id's start at 1, while container id is -1
            
        # test for two balls to see why not colliding 
        #self._balls.append(ba.Ball(mass, radius, [0,0], [0,0], ball_id=10))
        #self._balls.append(ba.Ball(mass, radius, [-5,0], [1,0], ball_id=11))

        
    # access methods    
    def num_balls(self):
        return self._nballs
    def balls(self):
        return self._balls
    def ball_rad(self):
        return self._radius
    def ball_mass(self):
        return self._mass
    def vfactor(self):
        return self._vfactor
    def container(self):
        return self._container
    def container_rad(self):
        return self._container_radius
    def sim_time(self):
        return self._simtime
    def momt_imparted(self):
        return self._momt_imparted
    
    
    def system_ke(self):
        """
        Calculates total kinetic energy (in SI units) for the system by summing
        the kinetic energies of all the Ball objects obtained from
        Ball's native method

        Returns
        -------
        total : float
            Total system kinetic energy.

        """
        total = self._container.ke()
        # needed for temperature calculation 
        for b in self._balls:
            total += b.ke()
        
        # even when this is used in temperature where only the ke
        # of the balls is needed, the ke of the container is negligible
        
        return total
    
    def system_momt(self):
        """
        Calculates total momentum (in SI units) of the system by summing
        the momenta of all the Ball objects obtained from
        Ball's native method

        Returns
        -------
        total : float
            Total system momentum.

        """
        total = self._container.momt()
        
        for b in self._balls:
            total += b.momt()
        
        return total
    
    
    
    ### task 9
    def centre_distance_hist(self, bins = 10):
        """
        Plots a histogram of distance the balls extend from the centre
        should be run after calling run method with many frames
        (> 100) to show stabilised system.

        Parameters
        ----------
        bins : int, optional
            Desired number of bins. The default is 10.

        Returns
        -------
        None.

        """
        
        
        dist = []
        for b in self._balls:
            dist.append(np.linalg.norm(b.pos()))
        
            
        fig, ax = plt.subplots()
        ax.set_title("Distance balls extend from centre")
        ax.set_xlabel("Distance of ball from centre (m)")
        ax.set_ylabel("Number of balls")
        # freedman - diaconis doesn't seem to work well here
        # doesn't give enough bins so set manually
        ax.hist(dist, bins, label = "Container radius 10m") 
        # using 10m as that's the radius I use for the report
        ax.legend()
        fig.savefig("distance_balls_from_centre", dpi = 500)
        plt.show()
            
            
    def dist_between_balls_hist(self):
        """
        Plots a histogram of distance between pairs of balls
        should be run after calling run method with many
        frames (>100) to show stabilised system.
        
        The Freedman-Diaconis Rule is used to calculate the optimum
        number of bins required to accuratel represent the data

        Returns
        -------
        None.

        """
        
        
        # create list of tuples of all possible combinations of objects
        # combinations returns all possible combinations of balls with no repeats
        # in tuples
        pair_balls = combinations(self._balls, 2)
        list_pair_balls = list(pair_balls)
        dist = []
        
        for pair in list_pair_balls:
            dist.append(np.linalg.norm(pair[0].pos() - pair[1].pos()))
        
        
        fig, ax = plt.subplots()
        ax.set_title("Inter-ball Separation")
        ax.set_xlabel("Distance between balls (m)")
        ax.set_ylabel("Number of pairs")
        ax.hist(dist, freedman_diaconis(dist, len(list_pair_balls)),\
                label = "Container radius 10m")
        # using 10m as that's the radius I use for the report
        ax.legend()
        fig.savefig("distance_between_balls", dpi = 500)
        plt.show()
        
        
        
    
    def get_available_positions(self, r_container, r_ball, nballs):
        """
        Obtains psuedo-random, non-overlapping, positions for the constructor
        to initialise balls in.
        
        Uses concentric circles in conatiner, with predetermined positions 
        for balls and which of these positions is filled is decided 
        pseudo-randomly using numpy.random.choice.

        Parameters (all in SI units where necessary)
        ----------
        r_container : float
            Container radius
        r_ball : float
            Ball radius
        nballs : int
            Desired number of balls for the simulation 

        Raises
        ------
        Exception
            Too many balls of this size for this conatiner size, 
            max no. of balls allowed for these parameters is len(positions)
            

        Returns
        -------
        positions_out : list
            Pseudo-randomly selected positions for the ball centres to be 
            used in the constructor when initialising the balls

        """
        
        # task 8
         
        
        r_container = abs(r_container)  # need +ve in subsequent calculations
        positions = []
        # ensures balls on diff circles don't overlap
        # container radius/ ball diameter gives balls just touching
        # so multiply by 0.75 to give space
        num_circles = int(np.floor(r_container/(2 * r_ball)) * 0.75)
        
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

        
    def next_collision(self):
        """
        Advances the system to the time of the next collision and collides the 
        two balls involved in the next collision.
        
        Does this by using itertools.combinations to create a list of all possible
        combinations of 2 Ball objects involved in the simulation. This list then 
        has a corresponding list of the times to collision for each pair, calculated
        from Ball's native method, from which the lowest time is taken. This time
        denotes the time of the next collision, so all balls are advanced to this time
        and the pair the time corresponds to are collided.
        
        Returns
        -------
        None.

        """
        
        # task 7
    
        # list balls and container
        objects = [*self._balls, self._container]
        
        # create list of tuples of all possible combinations of objects
        # combinations returns all possible combinations of balls with no repeats
        # in tuples
        list_pair_balls = list(combinations(objects, 2))
        times = []

        # find time to next collision for each pair of balls and append to times
        for ball1, ball2 in list_pair_balls:
            dt = ball1.time_to_collision(ball2)
            if dt is not None:
                 times.append(dt)
            else:
                # having None's in times for balls that don't collide creates errors
                # instead, append large number that won't be taken as the min time
                times.append(1e5)
        
        # find next time to collision (min time)
        next_t = min(times)
        
        # index of next time to collision (the minimum time)        
        index_next_t = times.index(min(times))
        
        # move all balls to next collision
        for ball in objects: # includes container too to avoid errors
            ball.move(next_t)
        
        # update simulation time     
        self._simtime += next_t
        
        next_ball1 = list_pair_balls[index_next_t][0]
        next_ball2 = list_pair_balls[index_next_t][1]
        # make balls corresponding to that time collide
        print("Colliding %s with %s" % (next_ball1, next_ball2))
        next_ball1.collide(next_ball2)
        
        # update momentum imparted to container
        # multiply by 2 as elastic collision so will leave 
        # with same mag of vel in opposite direction
        
        # only update if container was involved in collision
        if next_ball1.ball_id() == self._container.ball_id():
            self._momt_imparted += 2 * next_ball2.m() * \
                np.linalg.norm(next_ball2.vel())
        elif next_ball2.ball_id() == self._container.ball_id():
            self._momt_imparted += 2 * next_ball1.m() * \
                np.linalg.norm(next_ball1.vel())
        else:
            pass
        

        

    def run(self, num_frames, animate = False, values = False, consv_graph = False, \
            pressure_graph = False):
        """
        Runs the simulation for a specified number of frames, calling next_collision
        each time.
        
        Also keeps track of temperature, average pressure on the 
        container, system kinetic energy, magnitude 
        of system momentum and times of each collision

        Parameters
        ----------
        num_frames : int
            Desired number of frames (collisions) that the simulation should run for
        animate : Boolean, optional
            Displays simulation animation if True. The default is False.
        values : Boolean, optional
            Returns temperature and average pressure on the 
            container at each frame if True. The default is False.
        consv_graph : Boolean, optional
            Plots graphs showing system kinetic energy and 
            system momentum magnitude against time to check that
            the system conserves energy and momentum if True. The default is False.
        pressure_graph : Boolean, optional
            Plots a graph of average pressure exerted on the container
            against time if True. The default is False.

        Returns
        -------
        T : list
            Temperature at each frame
        P : list
            Average pressure exerted on the container at each frame.

        """
        
        # task 5
        
        # task 10
        # append temp, pressure, system ke,
        # magnitude of system momentum to a list after every collision
        
        T = []
        P = []
        sys_ke = []
        sys_momt_mag = []
        
        # list of simulation times at each collision for graphs
        t = []
        
        # calculate container area once, rather than every frame
        A = np.pi * self._container_radius ** 2
        
        if animate:
            f = pl.figure()
            ax = pl.axes(xlim=(-self._container_radius, self._container_radius), ylim=(-self._container_radius, self._container_radius), aspect = 1)
            ax.add_artist(self._container.get_patch())
            for b in self._balls:
                ax.add_patch(b.get_patch())
        # saving fig to debug frame by frame        
        #pl.savefig('frame_0')
        for frame in range(num_frames):
            print("Frame: %s"%(frame+1))
            # the (n-1)th frame is the nth collision
            # each frame is a collision
            self.next_collision()
            
            # retrieve avg system KE after each collision
            # temperature comes from balls, so only average over balls
            # take system KE excl container
            
            T.append(self.system_ke()/(len(self._balls) * k)) # in Kelvin
            
            # calculate pressure after each collision
            # P = total momt imparted to container / sim time / container area
            
            P.append(self._momt_imparted / (self._simtime * A))
            
            # appends system KE and momt transferred to 
            # container after each collision
            sys_ke.append(self.system_ke())
            sys_momt_mag.append(np.linalg.norm(self.system_momt()))
            # appends time in simulation after each collision
            t.append(self._simtime)
            
            
            if animate: 
                pl.pause(0.01)
                # saving fig to debug frame by frame
                #pl.savefig('frame_{}'.format(frame+1))
                pl.show() # showing in for loop to display a continuous animation
            
        if consv_graph: #Pressure v time, sys_ke v time, momt v time
        
            # conservation plots
            fig, (axe, axm) = plt.subplots(1,2)
            plt.tight_layout()
            
            fitke, covke = np.polyfit(t, sys_ke, 0, cov = True)
            lineke = np.poly1d(fitke)
            axe.set_title("System Kinetic Energy vs Time")
            axe.set_xlabel("Time (s)")
            axe.set_ylabel("System KE (J)")
            axe.set_aspect('auto')
            axe.plot(t, sys_ke, 'bx')
            axe.plot(t, lineke(t), 'r')
            print("uncertainy on gradient is \
                  %.3g J/s" % (np.sqrt(covke[0][0])))
            
            fitm, covm = np.polyfit(t, sys_momt_mag, 0, cov = True)
            linem = np.poly1d(fitm)
            axm.set_title("Magnitude of System Momentum vs Time")
            axm.set_xlabel("Time (s)")
            axm.set_ylabel("Magnitude of System momentum (kgm/s)")
            axm.set_aspect('auto')
            axm.plot(t, sys_momt_mag , 'rx')
            axm.plot(t, linem(t), 'b')
            print("uncertainy on gradient is \
                  %.3g kgm/s^2" % (np.sqrt(covm[0][0])))
            fig.savefig("conservation_laws", dpi = 500, bbox_inches = "tight")
            plt.show()
            
        if pressure_graph:
            
            # check pressure v time
            
            fig, axp = plt.subplots()
            axp.set_title("Pressure vs Time")
            axp.set_xlabel("Time (s) ")
            axp.set_ylabel("Pressure (Pa)")
            axp.plot(t, P)
            
            fig.savefig("pressure_v_time", dpi = 500, bbox_inches = "tight")
            plt.show()
            
        if values:
            
            return T, P
        
            

    def __repr__(self):
        """
        Representation of a given Simulation object

        Returns
        -------
        string
            Description of Simulation object that can be used to reproduce the object

        """
        return """%s(nballs = %g, mass = %g, radius = %g, container_radius = %g 
    vfactor = %g, vel_sd = %g)""" % ("Simulation", self._nballs, self._mass, self._radius, \
    self._container_radius, self._vfactor, self._vel_sd)
        
    def __str__(self):
        """
        User-friendly description of Simulation object

        Returns
        -------
        string
            Description of Simulation object for user.
        
        """
        return """Simulation of %g balls with mass %g kg, radius %g m in a container of radius %g m""" \
        % (self._nballs, self._mass, self._radius, self._container_radius)
            
        
        
        
    
        
        
        
        
