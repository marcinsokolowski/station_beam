# Script developed by Christopher Lee (paper : https://ui.adsabs.harvard.edu/abs/2022PASA...39...42L/abstract)
# get_sigma_clip_manual.py <on pulse left> <on pulse right>

import sys
import math
import prof_utils
import matplotlib.pyplot as plt
import numpy as np

on_pulse_left = int(sys.argv[1])
on_pulse_right = int(sys.argv[2])

#filename = sys.argv[1]
filename = 'ascii_archive.txt'
profile, length = prof_utils.get_from_ascii(filename)
sigma, data_clipped = prof_utils.sigmaClip(profile)

if on_pulse_right > on_pulse_left:
	on_pulse_bins = on_pulse_right - on_pulse_left
	sum_on_pulse = sum(profile[on_pulse_left:on_pulse_right])
else:
	on_pulse_bins = on_pulse_right - on_pulse_left + length
	sum_on_pulse = sum(profile[0:on_pulse_right]) + sum(profile[on_pulse_left:])

off_pulse_bins = length - on_pulse_bins

max_flux_index = profile.index(max(profile))
plt.plot(np.linspace(0, len(profile), len(profile)), profile)
plt.axvline(x = on_pulse_right, ls='--', c='r')
plt.axvline(x = on_pulse_left, ls='--', c='r')
plt.savefig(filename.replace('.txt', '_on_pulse.png'), dpi=400, bbox_inches='tight')

print("%s %s %s" % (sigma, off_pulse_bins, sum_on_pulse))
