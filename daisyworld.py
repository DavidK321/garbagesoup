#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 19:19:16 2019

@author: tammas

This is an implementation of daisy world: a theoretical 
model of the Gaia hypothesis. This theory was initially developed by 
Lovelock (1983) to demonstrate the plausibility of living things 
interacting with, and regulating, their environment.

    Wood, A. J., G. J. Ackland, J. G. Dyke, H. T. P. Williams, and T. M. 
        Lenton, 2015: Daisywolrd: a Review. Rev. Geophys., 46, RG1001, 
        https://doi.org/10.1029/2006RG000217.
"""


import numpy as np
import matplotlib.pyplot as plt


from daisyworld_definitions import albedo, daisy_replicator, beta, euler, local_temp, planetary_temp

# Define constants and variables
alphaw = 0.01   # Cover fraction of white daisies
alphab = 0.01   # Cover fraction of black daisies
p = 1           # The fraction of habitable surface
alphag = p-alphaw-alphab # Cover fraction of bare ground
aw = 0.75       # Albedo of white daisies
ab = 0.25       # Albedo of black daisies
ag = 0.5        # Albedo of bare ground
gamma = 0.3     # The death rate 
S = 1000        # Solar constant (W/m^2)
maxn = 1000     # Maximum number of iterations
tol = 0.000001  # Tollerance of solution
luminosities = np.arange(0.5,1.6, 0.002) # Stelar luminosities
alphaw_out = np.ones(len(luminosities))*np.nan # Output variable for white
alphab_out = np.ones(len(luminosities))*np.nan # Output variable for black
temp_out = np.ones(len(luminosities))*np.nan   # Output variable for temp

# Main loop for changing luminosity
for i,L in enumerate(luminosities):
    # Set a minimum for cover fractions
    if alphaw<0.01: alphaw = 0.01
    if alphab<0.01: alphab = 0.01
    alphag = p-alphaw-alphab
    # Reset counters
    n = 0
    changew, changeb = 1,1
    # Run loop for daisy earth.
    while (n<maxn) and (changew>tol) and (changeb>tol):
        # Store the initial cover fractions
        sw,sb = alphaw, alphab
        # Planetary albedo
        planet_albedo = albedo(alphaw,alphab,alphag,aw,ab,ag)
        # Planetary temperature
        T = planetary_temp(S,planet_albedo, L=L)
        # Local temperature
        Tw = local_temp(planet_albedo,aw,T)
        Tb = local_temp(planet_albedo,ab,T)
        # Birth rate
        betaw = beta(Tw, optimum=280)
        betab = beta(Tb, optimum=310)
        # Change in daisies
        dawdt = daisy_replicator(alphaw, alphag, betaw, gamma)
        dabdt = daisy_replicator(alphab, alphag, betab, gamma)
        # Integrate
        alphaw = euler(alphaw, dawdt)
        alphab = euler(alphab, dabdt)
        alphag = p-alphaw-alphab
        n += 1
    # Store the output
    alphaw_out[i] = alphaw
    alphab_out[i] = alphab
    temp_out[i] = T

# Plot the results
# Cover fractions
white = plt.plot(luminosities,alphaw_out*100,'b', label='White')
black = plt.plot(luminosities,alphab_out*100,'k', label='Black')
plt.legend(loc='upper right')
plt.xlabel('Luminosity')
plt.ylabel('Surface cover %')
plt.title('Cover fractions')
plt.show()

# Planetary temperature
plt.figure()
plt.plot(luminosities,temp_out-273.15,'r')
plt.xlabel('Luminosity')
plt.ylabel('Temperature (°C)')
plt.title('Planetary temperature')
plt.show()

#""" 1)  Change the minimum and initial area fractions or each daisy to create three scenarios:
    # A scenario with only white daisies, 2) only black daisies and 3) another where both coexist.
    #What are the ranges of luminosity where daisies can exist for each scenario? Why? 
    
    #Both: Daisies exit fron luminosity 0.62 until 1.45, white daisies from 0.7 until 1.45, black dasies from 0.62 until 1.2
    #only white: Dasies exist from 0.8 until 1.45
    #only black: Dasies exist from 0.62 until 0.99
    # If both kind of dasies exist each kind of dasiy exist within a broader range of luminosity, 
    #because black dasies help white ones to reach their needed tempreture on the planet earlier and white dasies help
    # black ones to cool down the planet through the higher albedo, even within higher luminosity
    
    #2) Increase/decrease the albedos of the daisies, together and separately.
    #How does this affect the survival of the daisies and their ability to moderate the climate?
    
    # add 10% to both albedos: black dasies resist 0.2 luminosity more, white dasies cover more and more surface with higher luminosity
    # climate is similar until 1.4 luminosity and much more balenced after 1.4 aroundd 25 Degrees (old was around 60)
    # substract 10% from both albedos: black dasies only exists until 1.0 luminosity but from 0.62, white ones start existing at 0.62 luminosity already but only until 1.2
    #so all interfaces are moved 0.2 to lower luminosity
    # similar climate in the beginning, but tempreture gets exttremely hot earlier  
    
    #3) Create a parameter for each daisy type that describes the optimum reproduction temperature. 
    #(note: all calculations are done in degrees K). Pass these to the beta function and try creating scenarios 
    #where the daisies have different optimum temperatures. Does having different optimum temperatures improve 
    #survival for the daisies?
    
    #Set derivatives to zero
    
    #4)  Modify the death rate gamma. What does this parameter do and why? 
    
     #when you decrease the death rate dasies will resist higher luminosity
     #when you increase the death rate dasies will resist lower only unnder lower luminosity
     # but startingpoint of covering will remain always the same
     # it's because you have more or lower capasity of dasies, which can lower down the effects of luminosity with their albedo

    #5)  Copy the Daisyworld while loop to create a single simulation with a luminosity of 1, and initial cover fraction 
    #of alphaw=0.2 and alphab=0.5 . Run for as many generations needed until the cover fractions reach equilibrium. 
    #Compare the cover fractions and the planetary temperature. How do the cover fraction change? Why?
    
    # How would you plot the graph with only one stable Luminosity?