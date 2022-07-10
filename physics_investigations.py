#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 15:28:38 2020

@author: ali
"""
import simulation as sim
import balls as ba
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# note: all numerical values retrieved where using 
# np.random.seed(100) enabled in simulation.py for reproducability
# this only ensures reproducability for fresh consoles (where the cell
# in question hasn't been run before)

# define constants 

k = 1.38064852e-23
#%% task 11 pressure v temp

# needs to be external funct, not method as must initialise many simulations
def pressure_v_temp(vfactors, markers, line_colours, \
                    frames = 300, nballs = 100 \
                    , ball_radii = [], ball_mass = 1.0, container_radius = 10.0):
    """
    Function utilising the Simulation class.
    
    Plots average pressure exerted on container vs temperature for as many ball radii 
    as needed.

    Parameters
    ----------
    vfactors : list
        List of vfactors to be passed into the Simulation object.
        Produces a different temperature for each iteration
    markers : list
        matplotlib.pyplot markers for data plotting.
    line_colours : list
        matplotlib.pyplot colours for best fit lines.
    frames : int, optional
        Number of frames each iteration should run for. 
        Should be at least 100 to obtain good pressure data. The default is 300.
    nballs : int, optional
        Number of balls that should be initialised in each simulation
        iteration. The default is 100.
    ball_radii : list
        Ball radii to iterate through for each simulation. Default is []
    ball_mass : float, optional
        Mass of balls in simulations. The default is 1.0.
    container_radius : float, optional
        Radius of the container in the simulations. The default is 10.0.

    Returns
    -------
    grads : list
        Gradients of each fit for each ball radius.
    grad_uncs : list
        Uncertainties in gradients of each fit for each ball radius.

    """
    temps = [[], [], []]
    pressures = [[], [], []]
    
    #create simulations with different vels (and thereby temps)
    
    for i in range(len(vfactors)):
        # loop through with simulations for diff ball radii
        for j in range(len(ball_radii)):
            simul = sim.Simulation(nballs, ball_mass, ball_radii[j], \
                           container_radius, vfactors[i])
            T, P = simul.run(frames, values = True)
            # take last temp as temp is constant for a given vfactor 
            # take last pressure value as pressure eventually settles
            # can see this in a pressure v time plot
            temps[j].append(T[-1])
            pressures[j].append(np.mean(P[-10:-1])) # take mean of pressures
            # at end as fluctuates
    
    
    fig, ax = plt.subplots()
    
    
    ax.set_title(r"Pressure v Temperature with Container Radius $%s \, m$" % (container_radius))
    ax.set_xlabel(r"Temperature ($K$)")
    ax.set_ylabel(r"Pressure (${Pa}$)")
    ax.set_aspect('auto')
    
    
    grads = []
    grad_uncs = []

    
    for i in range(len(ball_radii)):
        fit, cov = np.polyfit(temps[i], pressures[i], 1, cov = True)
        
        line = np.poly1d(fit)
    
        
        grad = fit[0]
        grad_unc = np.sqrt(cov[0][0])
        
        grads.append(grad)
        grad_uncs.append(grad_unc)
        
        ax.plot(temps[i], pressures[i], markers[i])
        # couldn't break this line into multiple lines as label would mess up
        ax.plot(temps[i], line(temps[i]), line_colours[i], label = r"Ball radius $%s \, m$, gradient = %.3g $JK^{-1}m^{-2}$" % (ball_radii[i], grad))
        print(r"Gradient for ball radius %s $m$ is %.2g +- \
              %.1g $PaK^-1$" % (ball_radii[i], grad ,grad_unc))
        
    ax.legend()
    fig.savefig("pressure_v_temp", dpi = 500, bbox_inches = "tight")
    plt.show()
    
    return grads, grad_uncs

# plot PRODUCES TASK 11, 12 PLOT

vfactors = np.linspace(1, 20, 8)

markers = ["rx", "bx", "yx"]
line_colours = ["r", "b", "y"]

pvt = pressure_v_temp(vfactors, markers, line_colours, \
                      ball_radii = [1, 0.1, 0.01], frames = 300, nballs = 23)
    
# max 23 balls allowed for ball radius of 1m to not overlap

# task 11 (assuming all other variables are held constant in all statements below)
# increasing number of balls would increase pressure, but temperature would
# decrease

# increasing container area would reduce pressure, but temperature would
# increase

# these hold for an ideal gas


#%% task 12

# using ideal gas law for container radius 10m and 23 balls
# area is 2D volume

area = np.pi * 10.0 ** 2

# gradient for ideal gas 
# here N = 23
theoretical_grad = (23 * k)/ area 

# use gradient of ball radius 0.01m  as that has the most negligible ball volume
grad_error = abs(pvt[0][2] - theoretical_grad)
grad_unc = pvt[1][2]

# returns true if our gradient is within experimental accuracy for an ideal gas

agrees = pvt[1][2] > grad_error # false

amount_outside_experimental_range = grad_error - grad_unc # around 7e-25 difference
# this is probably due to low N, as ideal gas law assumes large N
# but otherwise, nearly ideal gas

#%% task 13 curve fit

def maxwell(v, sigma):
    """
    Returns probability densities for a 2D 
    Maxwell-Boltzmann distribution.

    Parameters
    ----------
    v : float
        Velocity magnitude.
    sigma : float
        Standard deviation of the Maxwell-Boltzmann distribution.

    Returns
    -------
    prob : float
        Probability density for a given pair of parameters.

    """
    expo = np.exp(-v ** 2/ (2 * sigma ** 2))
    # use v not v**2 as working in 2D
    coeff = (4 * np.pi * v) * (1/ (np.sqrt(2 * np.pi) * sigma)) ** 3
    prob = coeff * expo
    
    return prob

def vel_dist(frames = 300, nballs = 100, ball_radius = 0.1, \
                    ball_mass = 1.0, container_radius = 10.0, vfactor = 1.0,\
                    vel_sd = 1.0, nbins = 10):
    """
    Function utilising the Simulation class.
    
    Plots the velocity distribution pdf of the Ball objects in a Simulation
    object after the simulation has run for a given number of frames, with
    a fitted 2D Maxwell-Boltzmann distribution overlaid.

    Parameters (all in SI units where appropriate)
    ----------
    frames : int, optional
        Number of frames that the simulation should run
        for before taking the velocity distribution. Should
        be a high number to allow the system to reach thermodynamic
        equilibrium for the Maxwell-Boltzmann distribution to 
        apply. The default is 300.
    nballs : int, optional
        Number of balls in the simulation, should be a high
        number ideally so the system reaches thermodynamic
        equilibrium faster. The default is 100.
    ball_radius : float, optional
        Radius of the balls in the simulation. Should be 
        smalled compared to the container radius to better 
        approximate an ideal gas. The default is 0.1.
    ball_mass : float, optional
        Mass of the balls in the simulation. The default is 1.0.
    container_radius : float, optional
        Radius of the container in the simulation. The default is 10.0.
    vfactor : float, optional
        Factor to multiply the initial velocities of the balls 
        by. The default is 1.0.
    vel_sd : float, optional
        The standard deviation of the Gaussian distribution
        that the ball velocities are initially pseudo-randomly
        drawn from. Smaller standard deviations give better
        Maxwell-Boltzmann fits. The default is 1.0.
    nbins : int, optional
        Desired number of bins for the histogram. The default is 10.

    Returns
    -------
    sigma_op : float
        The optomised value of the standard deviation of the fitted
        Maxwell-Boltzmann distribution.
    cov : numpy array
        Covariance matrix container the variance of the sigma_op
    theoretical_variance : float
        The theoretical variance of the velocity distribution of the 
        balls, assuming it to be a 2D Maxwell-Boltzmann distribution.

    """
    # use small ball radius by default to approximate ideal gas
    # use small sd by default so that you can achieve a mawellian spread
    
    k = 1.38064852e-23
    
    simul = sim.Simulation(nballs, ball_mass, ball_radius, \
                           container_radius, vfactor, vel_sd)
    # run for many frames for ball velocities to randomise    
    values = simul.run(frames, values = True)
    
    # need temp for maxwell
    T = values[0][-1] # take last temp as temp constant 
    vels = []
    
    theoretical_variance = (k * T)/ ball_mass
    
    for ball in simul.balls():
        vels.append(np.linalg.norm(ball.vel()))
    
    

    fig, ax = plt.subplots()
    ax.set_title(" Normalised velocity distribution of balls")
    ax.set_xlabel("Velocity (m/s)")
    ax.set_ylabel("Probability Density")
    
    
    # normalise with density = 1
    # not using freedman diaconis here as gives too few bins
    # more bins allows us to see the shape of the graph more clearly
    p, bins, patches = ax.hist(vels, nbins, density = 1)
    
    vel_test_fit = np.linspace(0, 5, nbins) # For maxwell fit
    
    sigma_guess = np.sqrt((k * T)/ ball_mass)
    
    sigma_op, cov = curve_fit(maxwell, vel_test_fit, p, p0 = sigma_guess )
    
    vel_test = np.linspace(0, 5, 100) # more samples to give smooth curve
    
    maxwell_fit = maxwell(vel_test, sigma_op)
    

    ax.plot(vel_test, maxwell_fit, label = "Maxwell- Boltzmann distribution")
    ax.legend()
        
    fig.savefig("vel_dist", dpi = 500, bbox_inches = "tight")
    plt.show()
    
    return sigma_op, cov, theoretical_variance
    
mb = vel_dist(frames = 300, nballs = 300, ball_radius = 0.01, nbins = 20) # PRODUCES TASK 13 PLOT

#%% task 13 compare variances 

theoretical_variance = mb[2]
empirical_variance = mb[0]
empirical_variance_unc = np.sqrt(mb[1][0][0])

print(theoretical_variance, empirical_variance, empirical_variance_unc)



#%% Task 14 Van der Waal's law
## PRODUCES TASK 14 PLOT

# a = 0 because this term is responsible for attraction between 
# balls which we have not coded for 

# b corrects for particles having non-negligible volume compared to container
# b is the volume of one particle in the gas (Nb is the volume of the whole gas),
# which is subtracted from V in the eqn to account for the volume that the gas 
# particles take up

def vdw(container_radii, nballs, frames = 600, ball_radius = 0.1, \
                    ball_mass = 1.0, vfactor = 1.0):
    """
    Function utilising the Simulation class.
    
    Van Der Waal's Law for non-interacting particles.
    Plots kT/P against V (container area as working in 2D)
    which gives -b as y intercept.

    Parameters
    ----------
    container_radii : list
        Values of container radii to iterate through to produce simulations
        with different V's for x-axis data.
    nballs : int
        Desired number of balls in each simulation.
    frames : int, optional
        Desired number of frames that each simulation should
        run for. The default is 600.
    ball_radius : float, optional
        Radii of the balls in the simulation. The default is 0.1.
    ball_mass : float, optional
        Masses of the balls in the simulation. The default is 1.0.
    vfactor : float, optional
        Factor to multiply the initial velocities of the balls 
        by. The default is 1.0.

    Returns
    -------
    b_fitted : float
        Value of b obtained from the linear fit.
    b_fitted_unc : float
        Uncertainty in the value of b obtained from the linear fit.
    fit : numpy array
        Array containing the gradient and y intercept of the linear fit.
    cov : numpy array
        Array containing the covariance matrix of the gradient and y-intercept
        from the linear fit.

    """
    
    k = 1.38064852e-23 # define in func so can be run in any script
    pressures = []
    temps = []
    y_values = []
    volumes = [] #2D volumes are areas
    
    # loop through container volumes, gathering pressure and temp
    # data for each one
    for i in range(len(container_radii)):
        volumes.append(np.pi * container_radii[i] ** 2)
        
        simul = sim.Simulation(nballs, ball_mass, ball_radius, container_radii[i], \
                               vfactor)
        T, P = simul.run(frames, values = True)
        
        temps.append(T[-1]) # take last value as temperature constant 
        # as avg ke is constant due to elastic collisions
        pressures.append(np.mean(P[-100:-1])) # avg of last 100 pressure values
    
    # get list of values for y axis
    for i in range(len(pressures)):
        y_value = (k * temps[i])/ pressures[i]
        y_values.append(y_value)
    
    fit, cov = np.polyfit(volumes, y_values, 1, cov = True)
        
    line = np.poly1d(fit)
    
        
    b_fitted = -fit[1] # as b must be +ve but y intercept -ve
    b_fitted_unc = np.sqrt(cov[1][1])
        
    fig, ax = plt.subplots()    
    ax.set_title("Determining Van Der Waal's Parameter b")
    ax.set_xlabel(r"Container area $(m^2)$")
    ax.set_ylabel(r"$kT/P$ ($J/{Pa}$)")
    
    
    ax.plot(volumes, y_values, 'rx')
    ax.plot(volumes, line(volumes), \
            label = r"Ball radius $%s m$" % (ball_radius))
    ax.legend()
    
    plt.savefig("vdw", dpi = 500, bbox_inches = "tight")
    print("b = %.3g +- %.3g m^2" % (b_fitted, b_fitted_unc))
    # return fit and cov for debugging to check grad is 1/nballs
    
    return b_fitted, b_fitted_unc, fit, cov 

# test

container_radii = np.linspace(1, 1.1, 10) # keep container radii
# range small so that vdw still applies (as if too big, 
# gas starts to become more ideal)

# choose ball radius 0.1 so non-ideal gas
van_der_waal = vdw(container_radii, nballs = 23, ball_radius = 0.1)

# if you run with ball_radius = 1 and container radii between 10 and 10.1
# b is around 20x larger than ball area, as opposed to about 1.3x larger
# for the case above. This is because the bigger the balls you pack together
# the more empty area there is between themn

# b = 0.0413 +- 0.0132 m^2 
# this is 1.27x higher than ball area as expected due to area between balls
grad = van_der_waal[2][0]
print(grad) # 0.02 approx which is very close to 1/nballs
# so vdw is the eqn of state in this container radius range


#%% b gets larger for larger containers and balls

container_radii = np.linspace(10, 10.1, 10)

# choose ball radius 0.1 so non-ideal gas
van_der_waal = vdw(container_radii, nballs = 23, ball_radius = 1)



