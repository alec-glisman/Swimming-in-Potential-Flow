# SECTION: Dependencies

# External dependencies
import numpy as np                 # mathematical data structures
import sys                         # access system data
import os                          # access system data
from optparse import OptionParser  # Get user input

# Internal Dependencies (NOTE: Assuming file called from project base directory)
sys.path.insert(0, os.getcwd() + '/python')
from GSDUtil import GSDUtil  # noqa: E402

# !SECTION (Dependencies)


# SECTION: Parse user input options

parser = OptionParser()
parser.add_option("--GSD-path", dest="u_gsd_path",
                  help="Path to GSD file to be created",
                  metavar="string")
parser.add_option("--dt", dest="u_dt",
                  help="(dimensionless) time-step for numerical integration",
                  metavar="double")
parser.add_option("--R-avg", dest="u_R_avg",
                  help="(dimensionless) average inter-particle pair separation during articulation period",
                  metavar="double")
parser.add_option("--Z-height", dest="u_Z_height",
                  help="(dimensionless) Initial height of collinear swimmer above wall",
                  metavar="double")
parser.add_option("--phase-angle", dest="u_phase_angle",
                  help="phase angle between oscillating pairs",
                  metavar="double")
parser.add_option("--U0", dest="u_U0",
                  help="velocity oscillation amplitude",
                  metavar="double")
parser.add_option("--omega", dest="u_omega",
                  help="velocity oscillation frequency",
                  metavar="double")

# !SECTION (Parse user input options)


# SECTION: Parameters

# Integrator parameters
t = 0.0
tf = 25.0
num_steps_output = 1000

# Material parameters
fluid_density = 1.0
particle_density = 1.0

# Potential parameters
wca_epsilon = 0.0
wca_sigma = 0.0

# Particle parameters
n = 6
types = ['Swimmer', 'Image']
typeid = [0, 0, 0, 1, 1, 1]

# !SECTION (Parameters)


# SECTION: GSD Creation

def initializeGSD():
    global gsd_class

    gsd_class = GSDUtil(gsd_path, create_gsd=True)
    gsd_class.setLogParameters(dt, t, tf, tau, num_steps_output,
                               fluid_density, particle_density, wca_epsilon, wca_sigma)
    gsd_class.setParticleParameters(n, types, typeid)


def setInitialConditions():

    # position
    pos = np.zeros((n, 3), dtype=np.double)
    pos[0] = [R_avg, 0.0, Z_height]
    pos[1] = [0.0, 0.0, Z_height]
    pos[2] = [- R_avg + (U0 / omega) * np.sin(phase_angle), 0.0, Z_height]
    pos[3] = [R_avg, 0.0, -Z_height]
    pos[4] = [0.0, 0.0, -Z_height]
    pos[5] = [- R_avg + (U0 / omega) * np.sin(phase_angle), 0.0, -Z_height]

    # velocity NOTE[epic=Assumptions]: must be calculated in simulation system (C++)
    vel = np.zeros_like(pos, dtype=np.double)

    # acceleration NOTE[epic=Assumptions]: must be calculated in simulation system (C++)
    acc = np.zeros_like(pos, dtype=np.double)

    # output data
    gsd_class.setKinematics(pos, vel, acc)


def setSystemData():

    # get snapshot
    snapshot = gsd_class.snapshot

    # convert data to GSD expected type
    l_R_avg = np.array([R_avg], dtype=np.double)
    l_Z_height = np.array([Z_height], dtype=np.double)
    l_U0 = np.array([U0], dtype=np.double)
    l_omega = np.array([omega], dtype=np.double)
    l_phase_angle = np.array([phase_angle], dtype=np.double)

    # output data
    snapshot.log['swimmer/R_avg'] = l_R_avg
    snapshot.log['swimmer/Z_height'] = l_Z_height
    snapshot.log['swimmer/U0'] = l_U0
    snapshot.log['swimmer/omega'] = l_omega
    snapshot.log['swimmer/phase_shift'] = l_phase_angle

# !SECTION: (GSD Creation)


# SECTION: For use when being called from command line

if __name__ == "__main__":
    global gsd_path, dt, R_avg, Z_height, phase_angle, U0, omega, epsilon, tau

    # Parse user input
    options, remainder = parser.parse_args(sys.argv[1:])
    gsd_path = str(options.u_gsd_path)
    dt = np.double(options.u_dt)
    R_avg = np.double(options.u_R_avg)
    Z_height = np.double(options.u_Z_height)
    phase_angle = np.double(options.u_phase_angle)
    U0 = np.double(options.u_U0)
    omega = np.double(options.u_omega)

    epsilon = (U0 / omega) / R_avg
    tau = (2.0 * np.pi) / omega

    initializeGSD()         # Create GSD, set log parameters, set particle parameters
    setInitialConditions()  # set kinematics
    setSystemData()         # set data specific to given simulation system
    gsd_class.saveSnapshot()

# !SECTION (For use when being called from command line)