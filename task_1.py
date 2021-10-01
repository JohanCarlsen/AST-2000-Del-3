'''
EGEN KODE
'''

import numpy as np
import matplotlib.pyplot as plt
from numba import njit
import ast2000tools.constants as const
import ast2000tools.utils as utils
from ast2000tools.solar_system import SolarSystem
from ast2000tools.space_mission import SpaceMission

seed = utils.get_seed('antonabr')
system = SolarSystem(seed)
mission = SpaceMission(seed)

km_to_AU = 6.68458712e-9                        # constant to scale units
AU = const.AU
G_sol = const.G_sol
star_mass = system.star_mass
T_star = system.star_temperature
star_radius = system.star_radius * km_to_AU

def temp(r_mean, R_p):
    '''
    This function calculates the temperature from
    the mean distance from the planet and
    its radius
    '''
    return (T_star**4 / 4 * (star_radius / r_mean)**2)**(1/4)

# @njit
def orbit(planet_number):
    '''
    Same code as in part 2. Implementet to calculate
    one orbit of each planet
    '''
    R_p = system.radii[planet_number] * km_to_AU                                        # radius of planet in AU
    a = system.semi_major_axes[planet_number]                                           # semi major axis of planet in AU
    P = np.sqrt(4*np.pi**2*a**3 / (G_sol*(star_mass + system.masses[planet_number])))   # calculates one orbit of the planet
    total_time = P
    time_steps = int(total_time * 11000)
    dt = total_time / time_steps
    r = np.zeros((time_steps, 2))
    v = np.zeros((time_steps, 2))
    v_initial = system.initial_velocities[:, planet_number]                             # initial values for r and v are imported from ast2000tools
    r_initial = system.initial_positions[:, planet_number]
    v[0,:] = v_initial
    r[0,:] = r_initial
    r_norm = np.linalg.norm(r[0])                                                       # initial norm vector
    r_norm_new = np.zeros(time_steps)                                                   # array to contain new norms of vector r
    r_norm_new[0] = r_norm
    a = -G_sol*star_mass / r_norm**3 * r[0,:]                                           # initial acceleration given from gravitational formula
    for i in range(time_steps-1):
        '''
        This for-loop solves the diff
        equations using the leapfrog method
        '''
        r[i+1,:] = r[i,:] + v[i,:]*dt + 0.5*a*dt**2
        r_norm_new[i+1] = np.linalg.norm(r[i+1,:])
        a_new = -G_sol*star_mass / r_norm_new[i]**3 * r[i+1,:]
        v[i+1,:] = v[i,:] + 0.5*(a + a_new)*dt
        a = a_new
        v_norm = np.linalg.norm(v[i+1])
    r_mean = np.mean(r_norm_new)                                                        # finds the mean distance from the star for each planet
    print(f'Planet {planet_number} has a mean distance from the star at {r_mean:.3f} AU, and has a mean temperature of {int(temp(r_mean, R_p))} Kelvin')
    return r_norm_new

# for i in range(8):
#     orbit(i)

def hab_zone(T_p):
    '''
    This function calculates the boundries for
    the habitable zone
    '''
    return T_star**2*star_radius / (2*T_p**2)

lower_limit = 260   # lower and upper limit for the habitable zone
upper_limit = 390

# print(f'The boundires for the habitable zone ranges from {hab_zone(upper_limit):.3f} AU to {hab_zone(lower_limit):.3f} AU')
# print(f'With plus/minus 15 K, the boundries become {hab_zone(upper_limit+15):.3f} AU and {hab_zone(lower_limit-15):.2f} AU')

'''
Planet 0 has a mean distance from the star at 1.849 AU, and has a mean temperature of 349 Kelvin
Planet 1 has a mean distance from the star at 2.799 AU, and has a mean temperature of 284 Kelvin
Planet 2 has a mean distance from the star at 7.760 AU, and has a mean temperature of 170 Kelvin
Planet 3 has a mean distance from the star at 13.609 AU, and has a mean temperature of 128 Kelvin
Planet 4 has a mean distance from the star at 21.133 AU, and has a mean temperature of 103 Kelvin
Planet 5 has a mean distance from the star at 10.010 AU, and has a mean temperature of 150 Kelvin
Planet 6 has a mean distance from the star at 3.635 AU, and has a mean temperature of 249 Kelvin
Planet 7 has a mean distance from the star at 4.776 AU, and has a mean temperature of 217 Kelvin

The boundires for the habitable zone ranges from 1.484 AU to 3.339 AU
With plus/minus 15 K, the boundries become 1.376 AU and 3.76 AU
'''
