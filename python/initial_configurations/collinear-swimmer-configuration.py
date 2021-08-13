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
parser.add_option("--R_avg", dest="u_R_avg",
                  help="(dimensionless) average inter-particle pair separation during articulation period",
                  metavar="double")
parser.add_option("--phase_angle", dest="u_phase_angle",
                  help="phase angle between oscillating pairs",
                  metavar="double")
parser.add_option("--epsilon", dest="u_epsilon",
                  help="(dimensionless) relative oscillation amplitude compared to average inter-particle separation",
                  metavar="double")

# !SECTION (Parse user input options)


# SECTION: Parameters

# Swimmer parameters
omega = 2.0 * np.pi
tau = (2.0 * np.pi) / omega

# Integrator parameters
t = 0.0
tf = 1.0
num_steps_output = 1000

# Material parameters
fluid_density = 1.0
particle_density = 1.0

# Potential parameters
wca_epsilon = 0.0
wca_sigma = 3.0

# Particle parameters
n = 3
types = ['Swimmer']
typeid = [0, 0, 0]

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
    pos = np.zeros((3*n), dtype=np.float32)
    pos[0] = R_avg  # particle 0
    pos[6] = - R_avg + (U0 / omega) * np.sin(phase_angle)

    # velocity NOTE[epic=Assumptions]: must be calculated in simulation system (C++)
    vel = np.zeros_like(pos)

    # acceleration NOTE[epic=Assumptions]: must be calculated in simulation system (C++)
    acc = np.zeros_like(vel)

    # output data
    gsd_class.setKinematics(pos, vel, acc)


def setSystemData():

    # get snapshot
    snapshot = gsd_class.snapshot

    # convert data to GSD expected type
    l_R_avg = np.array([R_avg], dtype=np.float32)
    l_U0 = np.array([U0], dtype=np.float32)
    l_omega = np.array([omega], dtype=np.float32)
    l_phase_angle = np.array([phase_angle], dtype=np.float32)

    # output data
    snapshot.log['swimmer/R_avg'] = l_R_avg
    snapshot.log['swimmer/U0'] = l_U0
    snapshot.log['swimmer/omega'] = l_omega
    snapshot.log['swimmer/phase_shift'] = l_phase_angle

# !SECTION: (GSD Creation)


# SECTION: For use when being called from command line

if __name__ == "__main__":
    global gsd_path, dt, R_avg, phase_angle, epsilon, U0

    # Parse user input
    options, remainder = parser.parse_args(sys.argv[1:])
    gsd_path = str(options.u_gsd_path)
    dt = np.double(options.u_dt)
    R_avg = np.double(options.u_R_avg)
    phase_angle = np.double(options.u_phase_angle)
    epsilon = np.double(options.u_epsilon)

    # Calculate parameters
    U0 = epsilon * R_avg * omega

    initializeGSD()         # Create GSD, set log parameters, set particle parameters
    setInitialConditions()  # set kinematics
    setSystemData()         # set data specific to given simulation system
    gsd_class.saveSnapshot()

# !SECTION
