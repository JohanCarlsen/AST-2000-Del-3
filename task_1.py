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

km_to_AU = 6.68458712e-9
AU = const.AU
G_sol = const.G_sol
star_mass = system.star_mass
T_star = system.star_temperature
star_radius = system.star_radius * km_to_AU

def temp(r_mean, R_p):
    return (T_star**4 / 4 * (star_radius / r_mean)**2)**(1/4)

# @njit
def orbit(planet_number):
    R_p = system.radii[planet_number] * km_to_AU
    a = system.semi_major_axes[planet_number]
    P = np.sqrt(4*np.pi**2*a**3 / (G_sol*(star_mass + system.masses[planet_number])))
    total_time = P
    time_steps = int(total_time * 11000)
    dt = total_time / time_steps
    r = np.zeros((time_steps, 2))
    v = np.zeros((time_steps, 2))
    v_initial = system.initial_velocities[:, planet_number]
    r_initial = system.initial_positions[:, planet_number]
    v[0,:] = v_initial
    r[0,:] = r_initial
    r_norm = np.linalg.norm(r[0])
    r_norm_new = np.zeros(time_steps)
    r_norm_new[0] = r_norm
    a = -G_sol*star_mass / r_norm**3 * r[0,:]
    for i in range(time_steps-1):
        r[i+1,:] = r[i,:] + v[i,:]*dt + 0.5*a*dt**2
        r_norm_new[i+1] = np.linalg.norm(r[i+1,:])
        a_new = -G_sol*star_mass / r_norm_new[i]**3 * r[i+1,:]
        v[i+1,:] = v[i,:] + 0.5*(a + a_new)*dt
        a = a_new
        v_norm = np.linalg.norm(v[i+1])
    r_mean = np.mean(r_norm_new)
    print(f'Planet {planet_number} has a mean temperature of {int(temp(r_mean, R_p))} Kelvin')
    return r_norm_new

# for i in range(8):
#     orbit(i)

'''
Planet 0 has a mean temperature of 349 Kelvin
Planet 1 has a mean temperature of 284 Kelvin
Planet 2 has a mean temperature of 170 Kelvin
Planet 3 has a mean temperature of 128 Kelvin
Planet 4 has a mean temperature of 103 Kelvin
Planet 5 has a mean temperature of 150 Kelvin
Planet 6 has a mean temperature of 249 Kelvin
Planet 7 has a mean temperature of 217 Kelvin
'''
