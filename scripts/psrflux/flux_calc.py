# Script developed by Christopher Lee (paper : https://ui.adsabs.harvard.edu/abs/2022PASA...39...42L/abstract)
# e.g. python flux_calc.py <sig_obs> <sig_sim> <no of bins> <no of non-clipped bins> <pulsar period> <profile sum>

import sys
import math
import prof_utils

bins = int(sys.argv[3])
sig_obs = float(sys.argv[1])
sig_sim = float(sys.argv[2])*math.sqrt(bins)
n = int(sys.argv[4])
P0 = float(sys.argv[5])
Delta_t = P0 / bins
sig_sim_1 = float(sys.argv[2])
prof_sum = float(sys.argv[6])

u_sig_obs = 0 # uncertainty in rms noise from observation
u_sig_sim = 0.3 * sig_sim # uncertainty in simulated noise
u_Phi = math.sqrt(n) * sig_obs * Delta_t # uncertainty in unscaled sum

zeta = sig_sim / sig_obs # scaling factor
zeta_msok_test = float(sys.argv[2])/float(sys.argv[1])
print("DEBUG : zeta = %.8f" % (zeta_msok_test))

# profile, length = prof_utils.get_from_ascii('ascii_archive.txt')

Phi = prof_sum * Delta_t # unscaled sum
Phi_Jys = Phi * zeta # scaled sum

S = Phi_Jys / P0 # mean flux density
u_S = Phi_Jys / P0 * math.sqrt( (u_sig_obs/sig_obs)**2 + (u_sig_sim/sig_sim)**2 + (u_Phi/Phi)**2 ) # uncertainty in mean flux density

print("Standard deviation of profile noise: %s" % sig_obs)
print("Simulated noise: %s Jy" % sig_sim_1)
print("Number of phase bins: %s" % bins)
print("Number of non-clipped phase bins: %s" % n)
print("Period: %s s" % P0)
print("Mean flux density: %s Jy" % S)
print("Uncertainty in mean flux density: %s Jy" % u_S)

# double mean_flux_data = (sum_profile_data*gSigmaSimulated)/(n_bins*gOriginalRMS);
print("DEBUG : %.8f * %.8f / %d = " % (prof_sum,zeta,bins))
mean_flux_msok_test = prof_sum*zeta/bins
print("DEBUG : mean_flux_msok_test = %.8f" % (mean_flux_msok_test))