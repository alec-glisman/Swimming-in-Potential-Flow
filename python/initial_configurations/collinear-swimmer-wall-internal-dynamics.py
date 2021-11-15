"""Script to create GSD file for specified system

__author__ = "Alec Glisman"

Example:
    python3 python/initial_configurations/collinear-swimmer-wall-configuration.py --GSD-path=data.gsd --dt=1e-6 --R-avg=4.0 --Z-height=6.0 --phase-angle=-1.57079632679e+0 --U0=2.0 --omega=1.0 --image-system=1
"""

# SECTION: Dependencies

# External dependencies
import numpy as np                 # mathematical data structures
import sys                         # access system data
import os                          # access system data
from optparse import OptionParser  # Get user input

# Internal Dependencies (NOTE: Assuming file called from project base directory)
sys.path.insert(0, os.getcwd() + '/python')
from GSDUtilPy import GSDUtilPy  # noqa: E402

# !SECTION (Dependencies)


# SECTION: Parse user input options

parser = OptionParser()
parser.add_option("--GSD-path", dest="u_gsd_path",
                  help="Path to GSD file to be created",
                  metavar="string")

parser.add_option("--dt", dest="u_dt",
                  help="(dimensionless) time-step for numerical integration",
                  metavar="double")

parser.add_option("--ti", dest="u_ti",
                  help="(dimensionless) initial integration time",
                  metavar="double")
parser.add_option("--tf", dest="u_tf",
                  help="(dimensionless) final integration time",
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

parser.add_option("--number-bodies", dest="u_number_bodies",
                  help="number of bodies to simulate",
                  metavar="int")
parser.add_option("--image-system", dest="u_image_system",
                  help="Is the system under consideration an image system (about z-axis)",
                  metavar="int")

# !SECTION (Parse user input options)


# SECTION: Parameters

# Integrator parameters
num_steps_output = 10000

# Material parameters
fluid_density = 1.0
particle_density = 1.0

# Potential parameters
wca_epsilon = 0.0
wca_sigma = 0.0

# Particle parameters
num_particles_per_body = 3
types = ['constrained', 'locater']

# !SECTION (Parameters)


# SECTION: GSD Creation

def initializeGSD(gsd_path,
                  N, M,
                  dt, ti, tf, tau,
                  num_steps_output,
                  fluid_density, particle_density,
                  wca_epsilon, wca_sigma,
                  image_system):
    """Create GSD, set log parameters, set particle parameters

    Args:
        gsd_path (string): path to GSD file to create
        N (int): number of particles to simulate
        M (int): number of bodies to simulate
        num_particles_per_body (int): number of particles per body
        dt (double): integrator finite delta-t
        ti (double): integrator initial time
        tf (double): integrator final time
        tau (double): characteristic simulation time
        num_steps_output (int): number of frames to write to GSD during simulation
        fluid_density (double): fluid mass density
        particle_density (double): particle mass density
        wca_epsilon (double): WCA potential energy scale
        wca_sigma (double): WCA potential length scale
        image_system (int): Integer boolean deciding if simulation is an image system

    Returns:
        GSD_class: GSD helper class
    """
    # Set particle type-ids
    typeid = [None] * N

    for i in range(M):
        typeid[3 * i] = 1
        typeid[3 * i + 1] = 0
        typeid[3 * i + 2] = 0
    # Convert typeid to type required by GSD schema
    typeid = np.array(typeid, dtype=np.uint32)

    # Create GSD
    gsd_class = GSDUtilPy(gsd_path, create_gsd=True)
    # Set GSD log parameters
    gsd_class.setLogParameters(dt, ti, tf, tau,
                               num_steps_output,
                               fluid_density, particle_density,
                               wca_epsilon, wca_sigma,
                               image_system)
    # Set GSD particle parameters
    gsd_class.setParticleParameters(N, types=types, typeid=typeid)

    return gsd_class


def setInitialConditions(gsd_class,
                         N, num_particles_per_body,
                         Z_height,
                         image_system):
    """Set initial kinematics of GSD file

    Args:
        gsd_class (GSD_class): GSD helper class
        N (int): number of particles to simulate
        num_particles_per_body (int): number of particles per body
        Z_height (double): z-axis height of bodies
        image_system (int): Integer boolean deciding if simulation is an image system
    """

    # orientation
    quat = np.zeros((N, 4), dtype=np.double)
    no_rotation_quat = np.array([1.0, 0.0, 0.0, 0.0], dtype=np.double)

    for i in range(N):

        quat[i] = no_rotation_quat

    # position (only need to specify locater position, articulation calculated in C++ simulation)
    pos = np.zeros((N, 3), dtype=np.double)

    for i in range(M):

        Ni = 3 * M

        for j in range(num_particles_per_body):

            # "Real" particles in image system
            if ((i < (M/2)) and (image_system == 1)):

                pos[Ni + j] = [0.0, 0.0, Z_height]

            # "Image" particles in image system
            elif ((i < (M/2)) and (image_system == 1)):

                pos[Ni + j] = [0.0, 0.0, -Z_height]

            # All particles in non-image system
            else:
                pos[Ni + j] = [0.0, 0.0, Z_height]

    vel = np.zeros_like(pos, dtype=np.double)  # Calculated in C++ simulation

    acc = np.zeros_like(pos, dtype=np.double)  # Calculated in C++ simulation

    # output data
    gsd_class.setKinematics(quat, pos, vel, acc)


def setSystemData(gsd_class,
                  M,
                  R_avg, Z_height,
                  U0, omega, phase_angle):
    """Set system-specific parameters to GSD class

    Args:
        gsd_class (GSD_class): GSD helper class
        M (int): number of bodies to simulate
        R_avg (double): average spacing between particle pairs over course of simulation
        Z_height (double): z-axis height of body initially
        U0 (double): internal oscillation amplitude of body
        omega (double): internal oscillation frequency of body
        phase_angle (double): internal oscillation phase angle between particle pairs
    """

    # get snapshot
    snapshot = gsd_class.snapshot

    # convert data to GSD expected type
    l_R_avg = np.array([R_avg], dtype=np.double)
    l_Z_height = np.array([Z_height], dtype=np.double)
    l_U0 = np.array([U0], dtype=np.double)
    l_omega = np.array([omega], dtype=np.double)
    l_phase_angle = np.array([phase_angle], dtype=np.double)
    l_U_swim = np.zeros(6 * M, dtype=np.double)
    l_A_swim = np.zeros(6 * M, dtype=np.double)

    # output data
    snapshot.log['swimmer/R_avg'] = l_R_avg
    snapshot.log['swimmer/Z_height'] = l_Z_height
    snapshot.log['swimmer/U0'] = l_U0
    snapshot.log['swimmer/omega'] = l_omega
    snapshot.log['swimmer/phase_shift'] = l_phase_angle
    snapshot.log['swimmer/U_swim'] = l_U_swim
    snapshot.log['swimmer/A_swim'] = l_A_swim


def saveGSDFrame(gsd_class):
    """Save GSD frame to file

    Args:
        gsd_class (string): GSD helper class
    """

    gsd_class.saveSnapshot()

# !SECTION: (GSD Creation)


# SECTION: For use when being called from command line

if __name__ == "__main__":
    """Main method
    """

    # Parse user input
    options, remainder = parser.parse_args(sys.argv[1:])

    gsd_path = str(options.u_gsd_path)

    dt = np.double(options.u_dt)
    ti = np.double(options.u_ti)
    tf = np.double(options.u_tf)

    R_avg = np.double(options.u_R_avg)
    Z_height = np.double(options.u_Z_height)
    phase_angle = np.double(options.u_phase_angle)
    U0 = np.double(options.u_U0)
    omega = np.double(options.u_omega)

    # REVIEW: {1: image system, 0: regular system}
    image_system = int(options.u_image_system)

    M = int(options.u_number_bodies)
    N = num_particles_per_body * M

    tau = (2.0 * np.pi) / omega
    epsilon = (U0 / omega) / R_avg

    gsd_class = initializeGSD(gsd_path,
                              N, M,
                              dt, ti, tf, tau,
                              num_steps_output,
                              fluid_density, particle_density,
                              wca_epsilon, wca_sigma,
                              image_system)

    setInitialConditions(gsd_class,
                         N, num_particles_per_body,
                         Z_height,
                         image_system)

    setSystemData(gsd_class,
                  M,
                  R_avg, Z_height,
                  U0, omega, phase_angle)

    saveGSDFrame(gsd_class)

# !SECTION (For use when being called from command line)
