'''
EGEN KODE
Kode for å regne ut areal av solcellepanel
'''

import scipy.constants as const
import ast2000tools.utils as utils
from ast2000tools.solar_system import SolarSystem
from ast2000tools.space_mission import SpaceMission

seed = utils.get_seed('antonabr')
system = SolarSystem(seed)

sigma = const.sigma
AU = const.au
R_star = system.star_radius * 1e3
T_star = system.star_temperature
eff = 0.12
watt = 40
r = 6e11
r_planet_3 = 3.635 * AU

def F(r):
    return sigma*T_star**4*(R_star/r)**2

def area(r):
    return watt/(eff*F(r))

# print(f'At a distance {r:.0e} km from the star, we need an area of {area(r):.3} m^2')
print(f'For the lander to function at planet 6, we need an area of {area(r_planet_3):.3} m^2')


'''
At a distance 6e+11 km from the star, we need an area of 0.464 m^2
'''
