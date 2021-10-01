import numpy as np
import matplotlib.pyplot as plt
from numba import njit
import ast2000tools.constants as const
import ast2000tools.utils as utils
from ast2000tools.solar_system import SolarSystem
from ast2000tools.space_mission import SpaceMission
from ast2000tools.shortcuts import SpaceMissionShortcuts
import time

"""
SNARVEI!!, Men skrevet skjelettkode i bunnen for å vise
hvordan det ville vært implementert.
"""

# Snarvei
seed = utils.get_seed('antonabr')
system = SolarSystem(seed)
mission = SpaceMission(seed)
shcut = SpaceMissionShortcuts(mission, [14143])
total_thrust_force = 1.6e13*5.140997440620973e-09
mass_loss_pr_s = 1.6e13*1.007317606364996e-12

shcut.place_spacecraft_on_escape_trajectory(total_thrust_force, mass_loss_pr_s, 0.3, 100000, np.pi/2, 100)

# Skjelettkode
"""
EGEN SKJELETTKODE!!! Dette er hvordan vu har tenkt at koden ville vært implemetert. Denne kjører
altså ikke, (fordi koden fra 1F ikke fungerte og vi ikke får til å fikse den),
men dette er mer en tanke om hvordan vi hadde løst det.
"""
# Verdier som vi trenger å ha med videre
G = const.G
rotational_time = system.rotational_periods[0]*24*60*60     # Omgjort far dager til sekunder for hjemplaneten
planet_radius = system.radii[0]*1000        # Omgjort fra km til meter for hjemplaneten
spacecraft_mass = mission.spacecraft_mass   # Allerede oppgitt i kg
planet_mass = system.masses[0]*1.989e30     # Gjort om fra solmasse til kg for hjemplaneten
v_escape = np.sqrt(2*G*planet_mass / planet_radius)     # Unnslipningshastighet for hjemplaneten i m/s
from_m_to_AU = 1 / 149597870700      # Konverterer fra meter til AU
from_s_to_yr = 1 / (365*24*60*60)

planet_position_data = ...                # Data hentet for planetens posisjon fra task 2 fra Del 2
planet_velocity_data = ...                # Data hentet for planetens hastighet fra task 2 fra Del 2

@njit
def rocket_launch(speed, fuel, planet_launch_angle, excpeted_launch_time, time_of_launch, gravity=True):
    """
    Vi regner først ut rakettens posisjon og hastighet i forhold til planeten, og transformerer
    posisjonen og hastigheteten i forhold til stjernas posisjon senere gitt planetens posisjon og hastighet i en tid
    time_of_launch gitt i år (fra task 2 fra Del 2).
    Speed, fuel, planet_launch_angle, excpeted_launch_time oppgis i SI-enheter (m/s, kg, rad, s).
    """
    # Kan skru av og på gravitasjon for testing
    if gravity == True:
        gamma = G
        print('gravity is on')
    if gravity == False:
        gamma = 0
        print('gravity is off')

    """
    Gjør noen kule triks for å finne ut hvilken posisjon planeten har i tiden time_of_launch,
    kaller denne indeksen time_of_launch i indeks senere i skjelett-koden her:
    np.where(something) ... ...
    ... ... programming words ... ...
    """
    dt = 1e-4
    time_steps = int(time/dt)
    total_mass = fuel + spacecraft_mass
    fuel_consumption = 0

    t = np.linspace(0, time, time_steps)
    v = np.zeros((time_steps, 2))
    r = np.zeros((time_steps, 2))

    v0 = 2*np.pi*planet_radius / rotational_time

    # Arrays som bestemmer retning for initalbetingelser gitt vinkelposisjonen til raketten på planeten
    set_planet_launch_initial_velocity = v0*np.array([-np.sin(planet_launch_angle), np.cos(planet_launch_angle)])
    set_planet_launch_position = planet_radius*np.array([np.cos(planet_launch_angle), np.sin(planet_launch_angle)])
    # Setter initialbetingelser
    v[0] = set_planet_launch_initial_velocity
    r[0] = set_planet_launch_position
    v_norm = np.linalg.norm(v[0])
    r_norm = np.linalg.norm(r[0])
    r_hat = r[0] / r_norm
    # Vi antar at vi alltid skyter opp radielt utover
    a = thrust_force/total_mass*r_hat - gamma*planet_mass / r_norm**2*r_hat
    dt = 1e-4
    T = 0
    for i in range(time_steps-1):
        # Bruker Leap-frog-løkke
        r[i+1] = r[i] + v[i]*dt + 0.5*a*dt**2
        r_norm = np.linalg.norm(r[i+1])
        r_hat = r[i+1] / r_norm
        # Må huske å oppdatere massen for det neste tidssteget for neste akselerajson i leap-frog
        total_mass -= mass_loss_rocket*dt
        a_ipo = thrust_force/total_mass*r_hat - gamma*planet_mass / r_norm**2*r_hat     # Regner ut akselerasjonen i det neste tidssteget
        v[i+1] = v[i] + 0.5*(a + a_ipo)*dt
        a = a_ipo

        v_norm = np.linalg.norm(v[i+1])
        fuel_consumption += mass_loss_rocket*dt
        T += dt
        # Tester om vi går tom for drivstoff
        if total_mass <= spacecraft_mass:
            print('not enough fuel')
            print('time : ', T)
            print('r=', r[i+1])
            print('v=', v_norm)
            # Transformerer fra planetkoordinatsystemet til stjernekoordinatsystemet
            R = planet_position_data[time_of_launch] + r[i+1] * from_m_to_AU      # Gjør om til posisjon i stjernesystemet
            V = planet_velocity_data[time_of_launch] + v[i+1] * from_m_to_AU / from_s_to_yr # Gjør om til hastighet i stjernesystemet
            return t, fuel_consumption, total_mass - spacecraft_mass, v, r, V, R, T
        # Tester om raketten går under bakken på planeten
        if r_norm < planet_radius:
            r[i+1] = planet_radius*e_r      # Om raketten synker nedover i planeten setter jeg den tilbake til bakken
        # Tester om vi har nådd farten vi ønsker
        if v_norm >= speed:
            print('Sucessful launch in', T, 's')
            print('v=', v_norm)
            # Transformerer fra planetkoordinatsystemet til stjernekoordinatsystemet
            R = planet_position_data[time_of_launch] + r[i+1] * from_m_to_AU      # Gjør om til posisjon i stjernesystemet
            V = planet_velocity_data[time_of_launch] + v[i+1] * from_m_to_AU / from_s_to_yr # Gjør om til hastighet i stjernesystemet
            return t, fuel_consumption, total_mass - spacecraft_mass, v, r, V, R, T

    # Transformerer fra planetkoordinatsystemet til stjernekoordinatsystemet
    R = planet_position_data[time_of_launch] + r[i+1] * from_m_to_AU      # Gjør om til posisjon i stjernesystemet
    V = planet_velocity_data[time_of_launch] + v[i+1] * from_m_to_AU / from_s_to_yr # Gjør om til hastighet i stjernesystemet
    return t, fuel_consumption, total_mass - spacecraft_mass, v, r, V, R, T

t, mass_consumed, fuel_left, v, r, V, R, T = rocket_launch(speed=v_escape, fuel=initial_fuel, planet_launch_angle=np.pi/2, excpeted_launch_time=1200, time_of_launch=1)
print(f'consumption : {mass_consumed_boost}\nfuel left : {fuel_left}')
# Ville plottet data for å sjekke hvordan de ser ut
fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax1.plot(r[:,0], r[:,1])
ax2.plot(t, v)
plt.show()
