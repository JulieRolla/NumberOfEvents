import numpy as np #import numpy

from NumberOfEventsFunction import NumberOfEvents
import sys, imp
imp.reload(sys.modules['NumberOfEventsFunction'])


SecPerDay = 86400. #second in a day
KM2toCM2 = 1.e10

execfile('/Users/JulieRolla/Downloads/Proposal2018-master/comparison_plots/script.py')

'''
Energy Bins
'''
# Sets the energy bin edges. If this is changed, the acceptances and fluxes below must also be changed
Bin_edges=np.array([3.16227766e+15,   1.00000000e+16,   3.16227766e+16,
         1.00000000e+17,   3.16227766e+17,   1.00000000e+18,
         3.16227766e+18,   1.00000000e+19,   3.16227766e+19,
         1.00000000e+20,   3.16227766e+20])



'''
Fluxes
in units of cm^-2 sr^-1 eV^-1 s^-1
'''
## IceCube  (Brian found in his code, confirmed seperately)
IC_thuMu=icecube_thrumu_efe/Bin_edges
IC_combined=icecube_combined_efe/Bin_edges

## GZK (Brian found in his code)
#Ahlers
GZK_Ahlers=ahlers_limit_interp/Bin_edges

# Kotera would need to interpolate and then find flux at each energy bin. Can do if needed


'''
Acceptances
all converted to units of cm^2 sr
'''
### Brian's findings (took from limits: aeff = 2.3/limit/(365.*SecPerDay))
# this doesn't include an ln(10), 4*pi or an efficiency term (athough not sure if that matters)

#ARA
ARAsingleAcc_B = ara2_aeff_full * KM2toCM2
ARA37Acc_B= ARAsingleAcc_B*37

#ARIANNA
ARIANNAsingleAcc_B = arianna_aeff_full * KM2toCM2
ARIANNA1296Acc_B = 1296*ARIANNAsingleAcc_B

# Phased Array Detectors
PAsingleAcc_B=pa_aeff_full* KM2toCM2

#### Hybrid Combinations
# Change these to get numbers of each type
num_ara=1.
num_pa=0.
num_arianna=1.

hybridAcc_B=(num_ara*ARAsingleAcc_B)+(num_pa*PAsingleAcc_B)+(num_arianna*ARIANNAsingleAcc_B)


### Julie's acceptances from candidacy findings (From figure four https://arxiv.org/pdf/1410.7352.pdf then get Acceptance using thin target approximation)

#ARA
ARA37Acc_J = np.array([0, 1.08E+06,2.18E+07,1.84E+08,7.42E+08,2.33E+09,5.15E+09,1.03E+10,1.76E+10,3.02E+10,5.22E+10])

#ARIANNA
ARIANNA1296Acc_J = np.array([0,0,0, 2.84E+08,1.54E+09,5.35E+09,1.38E+10,2.99E+10,5.92E+10,1.11E+11,1.91E+11])


'''
Examples of function working
'''
#NumberOfEvents(Bin_edges, DetectorAcceptance, Flux, Days)
N= NumberOfEvents(Bin_edges, ARIANNA1296Acc_J, GZK_Ahlers, 365)
print(np.sum(N))


