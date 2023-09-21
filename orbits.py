import spiceypy as sp
import os 

in_fn = os.path.expanduser('~/Documents/Papers/SAPS_mission_concept/orbit/sc0_300.0alt_65.0inc_0.8sep_20000.bsp')
sp.furnsh(in_fn)
sp.spkpos
