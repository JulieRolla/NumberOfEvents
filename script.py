# -*- coding: utf-8 -*-
import numpy as np #import numpy
import matplotlib.pyplot as plt #import matplotlib
from matplotlib.pyplot import *
rcParams['mathtext.default'] = 'regular'
import math
from scipy.interpolate import interp1d, splrep, splev
from pylab import setp

SecPerDay = 86400. #second in a year
KM2toCM2 = 1.e10 #multiply by this to convert from km^2 to cm^2 (divide to go the other way)

#plots the projected limit from a combination of ARA, ARIANNA, and PA stations
#plots limits and effective areas side by side
#also, correctly accounts for relative livetimes between an autonomous and cabled station

acceptanceEnergies=np.array([1.000E+16,3.162E+16,1.000E+17,3.162E+17,1.000E+18,3.162E+18,1.000E+19,3.162E+19,1.000E+20,3.162E+20])
analysisAcceptanceARA2=np.array([6.77E-06,1.61E-06,1.64E-06,2.59E-06,5.98E-06,1.85E-05,6.05E-05,2.71E-04,1.34E-03,7.68E-03])
triggerAcceptanceARA2=np.array([1.28E-06,5.37E-07,6.69E-07,1.28E-06,3.46E-06,1.12E-05,4.28E-05,2.02E-04,1.07E-03,6.21E-03])
Bin_edges=np.array([  3.16227766e+15,   1.00000000e+16,   3.16227766e+16,
         1.00000000e+17,   3.16227766e+17,   1.00000000e+18,
         3.16227766e+18,   1.00000000e+19,   3.16227766e+19,
         1.00000000e+20,   3.16227766e+20])


'''
		What we want to do is extract an analysis level Effective Area * Steradians
		We are going to do this by backing it out from the published limit curves
		Here are the steps
		1) First, import the data from the limit curve
		2) Then, I rescale the limit curve so that it would represent 1 year
		3) Then, we hvae to interpolate the limit curve to every half decade, to go from logE = 15.5 to 20.5
		4) Then, we can convert from the limit curve to an effective area by multiplying by factors of 2.3, seconds, and converting km^2 to cm^2
'''


'''
		Below I comment a line-by-line example of how this works
'''

#import data
ara2_data = np.genfromtxt("limits/ara2_limit.csv",delimiter=',',skip_header=1,names=['energy','limit'])
#extract the 'limit' column, and rescale to 1 year
ara2_limit=ara2_data['limit'] * (445./365.) #take out livetime factor (includes 2 station correction)
#convert to an effective area in units of km^2
ara2_aeff = 2.3/ara2_limit/(365.*SecPerDay)/KM2toCM2
#these are the energies for the limit curve in logE
ara2_energy_logev = np.array([16.,16.5,17.,17.5,18.,18.5,19.,19.5,20,20.5])
#which we can unlog to just get the energy
ara2_energy = ara2_data['energy'] 
#
#	#now, we can set up an interpolator to interpolate across the curve
#	#we do this in log-log space, because it's easier to do that interpolation
#	ara2_interpolator = splrep(ara2_energy_logev, np.log10(ara2_limit),k=2)
#	#we want to define the log(energies) we want have interpolated limits at
#	ara2_energy_logev_interp = np.array([15.5,16.,16.5,17.,17.5,18.,18.5,19.,19.5,20,20.5])
#	#convert from logE to E
#	ara2_energy_interp = np.power(10.,ara2_energy_logev_interp)
#	#actually evaluate the spline interpolator at logE
#	ara2_limit_interp = np.power(10.,splev(ara2_energy_logev_interp, ara2_interpolator))
#	#this will be the actual ara2_aeff at the interpolated energies
#	ara2_aeff_interp = 2.3/ara2_limit_interp/(365.*SecPerDay)/KM2toCM2

ara2_energy_full=np.concatenate([[3.16000000e+15],ara2_energy])
ara2_energy_logev = np.array([15.5,16.,16.5,17.,17.5,18.,18.5,19.,19.5,20,20.5])
ara2_limit_full=np.concatenate([[0],ara2_limit])
ara2_aeff_full=np.concatenate([[0],ara2_aeff])

'''
		To estimate the sensitivity of the phased array, we are going to cheat just a touch
		We are going to cheat in the sense that we're just going to slide the ARA2 curve to the left
		on an EF(E) plot by a factor of two. Then we repeat the whole business above of interpolating and converting.
'''
#
#	pa_limit=ara2_limit
#	pa_energy = ara2_energy/2 #scoot the curve over
#	pa_energy_logev = np.log10(pa_energy)
#	pa_aeff = 2.3/pa_limit/(365.*SecPerDay)/KM2toCM2
#
#	pa_interpolator = splrep(pa_energy_logev, np.log10(pa_limit),k=1)
#	pa_energy_logev_interp = np.array([15.5,16.,16.5,17.,17.5,18.,18.5,19.,19.5,20,20.5])
#	pa_energy_interp = np.power(10.,pa_energy_logev_interp)
#	pa_limit_interp = np.power(10.,splev(pa_energy_logev_interp, pa_interpolator))
#pa_aeff_interp = 2.3/pa_limit_interp/(365.*SecPerDay)/KM2toCM2

pa_limit=ara2_limit
pa_energy = ara2_energy/2 #scoot the curve over
pa_energy_logev = np.log10(pa_energy)
pa_aeff = 2.3/pa_limit/(365.*SecPerDay)/KM2toCM2
	
pa_interpolator = splrep(pa_energy_logev, np.log10(pa_limit),k=1)
pa_energy_logev_interp = np.array([16.,16.5,17.,17.5,18.,18.5,19.,19.5,20,20.5])
pa_energy_interp = np.power(10.,pa_energy_logev_interp)
pa_limit_interp = np.power(10.,splev(pa_energy_logev_interp, pa_interpolator))
pa_aeff_interp = 2.3/pa_limit_interp/(365.*SecPerDay)/KM2toCM2

pa_energy_full=np.concatenate([[3.16000000e+15],pa_energy_interp])
pa_energy_logev = np.array([15.5,16.,16.5,17.,17.5,18.,18.5,19.,19.5,20,20.5])
pa_limit_full=np.concatenate([[0],pa_limit_interp])
pa_aeff_full=np.concatenate([[0],pa_aeff_interp])


testbed_data = np.genfromtxt("limits/testbed_limit.csv",delimiter=',',skip_header=1,names=['energy','limit'])
testbed_limit=testbed_data['limit'][0:-1] * (415./365.) #strip off last point to match ARIANNA, and take out livetime
testbed_aeff = 2.3/testbed_limit/(365.*SecPerDay)/KM2toCM2
testbed_energy_logev = np.array([17,17.5,18.,18.5,19,19.5,20,20.5])
testbed_energy = np.power(10.,testbed_energy_logev)

#testbed_interpolator = splrep(testbed_energy_logev[1:], np.log10(testbed_limit[1:]),k=3)
testbed_energy_logev_full = np.array([15.5,16.,16.5,17.,17.5,18.,18.5,19.,19.5,20,20.5])
testbed_energy_full = np.power(10.,testbed_energy_logev_full)
testbed_limit_full = np.concatenate([[0,0,0],testbed_limit])
testbed_aeff_interp = np.concatenate([[0,0,0],testbed_aeff])
	
	
hra3_data = np.genfromtxt("limits/hra3_limit.csv",delimiter=',',skip_header=1,names=['energy','limit'])
hra3_limit=hra3_data['limit']/(hra3_data['energy']/1e9) * (170./365.)#need to divide out the GeV, and take out livetime (includes 3 station correction)
hra3_aeff = 2.3/hra3_limit/(365.*SecPerDay)/KM2toCM2
hra3_energy_logev = np.array([17,17.5,18.,18.5,19,19.5,20,20.5])
hra3_energy =  np.power(10.,hra3_energy_logev)

#arianna_interpolator = splrep(hra3_energy_logev[1:], np.log10(hra3_limit[1:]),k=3)
arianna_energy_logev_full = np.array([15.5,16.,16.5,17.,17.5,18.,18.5,19.,19.5,20,20.5])
arianna_energy_full= np.power(10.,arianna_energy_logev_full)
arianna_limit_full = np.concatenate([[0,0,0],hra3_limit])
arianna_aeff_full = np.concatenate([[0,0,0],hra3_aeff])

	
'''
		Okay, now we need to import the IceCube measurement
		They have two flux measurements, one from their "thru-mu" and one from "combind-fit"
'''

icecube_energy_logev=np.array([15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,20.5])
icecube_energy = np.power(10.,icecube_energy_logev)

	#E is in eV
	#returns number/eV/cm^2/s/sr

#we can first have the IceCube thru-mu fluxes
#we can have the nomina, and upper and lower 1 sigma possibilities
def icecube_thrumu_function(E):
	return 3.03 * ((E/1.e14)**-2.19) * 1e-27
def icecube_thrumu_upper_function(E):
	return 3.81 * ((E/1.e14)**-2.09) * 1e-27
def icecube_thrumu_lower_function(E):
	return 2.34 * ((E/1.e14)**-2.29) * 1e-27
	
def icecube_combined_function(E):
	return 6.7 * ((E/1.e14)**-2.50) * 1e-27
def icecube_combined_upper_function(E):
	return 7.8 * ((E/1.e14)**-2.41) * 1e-27
def icecube_combined_lower_function(E):
	return 5.5 * ((E/1.e14)**-2.59) * 1e-27
#to compute to a curve on EF(E) we have to multiply by energy
icecube_thrumu_efe = icecube_thrumu_function(icecube_energy) * icecube_energy
icecube_thrumu_upper_efe = icecube_thrumu_upper_function(icecube_energy) * icecube_energy
icecube_thrumu_lower_efe = icecube_thrumu_lower_function(icecube_energy) * icecube_energy
icecube_combined_efe = icecube_combined_function(icecube_energy) * icecube_energy
icecube_combined_upper_efe = icecube_combined_upper_function(icecube_energy) * icecube_energy
icecube_combined_lower_efe = icecube_combined_lower_function(icecube_energy) * icecube_energy


'''
		Okay, now we need to pull in some theoretical limits
		We'll do ahlers and kotera
'''
ahlers_data = np.genfromtxt("limits/ahlers_2012.csv",delimiter=',',skip_header=1,names=['energy','flux'])
ahlers_energy_logev = ahlers_data['energy']
ahlers_energy = np.power(10.,ahlers_energy_logev)
ahlers_limit_log = ahlers_data['flux']
ahlers_limit = np.power(10.,ahlers_limit_log) 

ahlers_interpolator = splrep(ahlers_energy_logev, np.log10(ahlers_limit),k=1)
ahlers_energy_logev_interp = np.array([15.5,16.,16.5,17.,17.5,18.,18.5,19.,19.5,20,20.5])
ahlers_energy_interp = np.power(10.,ahlers_energy_logev_interp)
ahlers_limit_interp = np.power(10.,splev(ahlers_energy_logev_interp, ahlers_interpolator))

	
	
	
	



kotera_max_data = np.genfromtxt("limits/kotera_max.csv",delimiter=',',skip_header=1,names=['energy','flux'])
kotera_max_energy_logev = kotera_max_data['energy']
kotera_max_energy = np.power(10.,kotera_max_energy_logev)
kotera_max_limit_log = kotera_max_data['flux']
kotera_max_limit = np.power(10.,kotera_max_limit_log) 
kotera_min_data = np.genfromtxt("limits/kotera_min.csv",delimiter=',',skip_header=1,names=['energy','flux'])
kotera_min_energy_logev = kotera_min_data['energy']
kotera_min_energy = np.power(10.,kotera_min_energy_logev)
kotera_min_limit_log = kotera_min_data['flux']
kotera_min_limit = np.power(10.,kotera_min_limit_log)
#to make a fill between for the kotera, we need to interpolate the min data to the max data energy points
kotera_min_interpolator = splrep( kotera_min_energy_logev , kotera_min_limit_log,k=1)
kotera_min_limit_interp = np.power(10.,splev(kotera_max_energy_logev,kotera_min_interpolator))

'''
	Okay, now we can design our new "ideal" detectors
	We will want it to be some combination of ARA stations "num_ara"
	phased array enhances ARA stations "num_pa"
	arianna stations "num_arianna"
	and number of years we want them to run for "num_years"
'''

num_ara=1.
num_pa=0
num_arianna=1.
num_years=5.

hybrid_energy_logev = np.array([15.5,16.,16.5,17.,17.5,18.,18.5,19.,19.5,20,20.5])
hybrid_energy = np.power(10.,hybrid_energy_logev)
'''
We want to break up our calculation into two pieces
An "autonomous" piece, which is composed of just ARIANNA stations
An "cabled" piece, which is composed of ARA + PA stations
'''

#so these are the effective areas in km^2 * sr
#this is the part that sets the total cabled aeff
hybrid_aeff_cabled = (num_ara*ara2_aeff_full) + (num_pa*pa_aeff_full)	
#this is the part that sets the total autonomous aeff
hybrid_aeff_autonomous = (num_arianna*arianna_aeff_full)
#this is the part that sets the total aeff	
hybrid_aeff = hybrid_aeff_cabled + hybrid_aeff_autonomous

#to compute the expected limit, we now need to take into accountt the *livetimes*, convert to seconds, and to cm^2
hybrid_aeff_cabled_with_livetime = hybrid_aeff_cabled*(num_years*365.*SecPerDay)*KM2toCM2
hybrid_aeff_autonomous_with_livetime = hybrid_aeff_autonomous*(num_years*365.*SecPerDay*0.6)*KM2toCM2
hybrid_aeff_with_livetime = hybrid_aeff_cabled_with_livetime + hybrid_aeff_autonomous_with_livetime
	
#now we can actually compute the limit setting power
#hybrid_cabled_limit = 1./hybrid_aeff_cabled_with_livetime
#hybrid_cabled_limit[hybrid_cabled_limit == np.inf] = 0
#
#hybird_autonomous_limit = 1./hybrid_aeff_autonomous_with_livetime
#hybird_autonomous_limit[hybird_autonomous_limit == np.inf] = 0
#
#hybrid_limit = 1./hybrid_aeff_with_livetime
#hybrid_limit[hybrid_limit == np.inf] = 0
