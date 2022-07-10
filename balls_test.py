#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 11:44:00 2020

@author: ali
"""

import balls as ba

#%% test initialisation 

b = ba.Ball()
c = ba.Ball(1, 3, [7,8])

# test representation and string

print(repr(b))
print(b) #works

print(repr(c))
print(c) #works

#%% test getter methods

print(b.pos())
print(b.vel()) #works


#%% test setter methods
a = ba.Ball(pos = [1,2], vel = [1,2])

# move
a.move(1)
print(a) # works 

#%% test time to collision exceptions

# divide by 0 exception 
d = ba.Ball(pos = [1, 2], vel = [1, 2])
e = ba.Ball(pos = [1, 2], vel = [1, 2])

print(d.time_to_collision(e))

#raises dividing by 0 exception as expected
# Exception: "A relative velocity of 0 will yield infinite dt as 
#they will never collide

# input expection 
b1 = ba.Ball(pos = [1,2], vel = [1,2,3])
b2 = ba.Ball(pos = [1,2,3], vel = [1,2]) #works

# correct exception raised
# Exception: vel only takes 2D vectors

#%% test time to collision 

# where they will collide
d = ba.Ball(pos = [1, 0], vel = [1, 0])
e = ba.Ball(pos = [2, 0], vel = [0, 0])

print(d.time_to_collision(e))
# gives effectively 0.8 s which is expected instead of 1 s as when they collide
# their relative radius is 0.2 m, so that's 0.2m they don't have to travel

#%% where they won't collide

x = ba.Ball(pos = [1, 0], vel = [0, 0])
y = ba.Ball(pos = [2, 0], vel = [0, 1])

print(x.time_to_collision(y))
# returns None as expected

#%% where they've collided in past (moving away from each other)

w = ba.Ball(pos = [1, 0], vel = [1, 0])
v = ba.Ball(pos = [-1, 0], vel = [-1, 0])

print(w.time_to_collision(v))
# returns None as expected

#%% test collide 

# d and e will collide
# original vels [1. 0.] [0. 0.]

d.collide(e)

print(d.vel(), e.vel())
# [0. 0.] [1. 0.]
# their velocities have switched, 
#exactly what we expect from a head on elastic collision



