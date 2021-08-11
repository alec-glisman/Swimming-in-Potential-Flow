# External dependencies
import numpy as np  # mathematical data structures
import gsd.hoomd  # Simulation data-structure
from optparse import OptionParser  # Get user input

# Parse user input
parser = OptionParser()
parser.add_option("--dt", dest="u_dt",
                  help="(dimensionless) time-step for numerical integration", metavar="double")
parser.add_option("--R_avg", dest="u_R_avg",
                  help="(dimensionless) average inter-particle pair separation during articulation period", metavar="double")
parser.add_option("--phase_angle", dest="u_phase_angle",
                  help="phase angle between oscillating pairs", metavar="double")
parser.add_option("--epsilon", dest="u_epsilon",
                  help="(dimensionless) relative oscillation amplitude compared to average inter-particle separation", metavar="double")


# SECTION: Parameters
# Simulation parameters
gsd_name = "test.gsd"
step = np.array([0], dtype=np.uint64)
num_dimensions = np.array([3], dtype=np.uint8)

# Integrator parameters
dt = np.array([1e-5], dtype=np.float32)
t = np.array([0.0], dtype=np.float32)
tf = np.array([1.0], dtype=np.float32)
tau = np.array([1.0], dtype=np.float32)
num_steps_output = np.array([1000], dtype=np.uint64)

# Length scales
a = 1.0  # NOTE[epic=assumptions,seq=1]: a must always be 1
wall_separation = 10 * a
box_size = 30 * a

# Material parameters
fluid_density = np.array([1.0], dtype=np.float32)
particle_density = np.array([1.0], dtype=np.float32)

# Potential parameters
wca_epsilon = np.array([0.0], dtype=np.float32)
wca_sigma = np.array([3.0*a], dtype=np.float32)

# Initial positions
R_1 = np.array([5.0, 0.0, wall_separation], dtype=np.float32)
R_2 = np.array([0.0, 0.0, wall_separation], dtype=np.float32)
R_3 = np.array([-5.0, 0.0, wall_separation], dtype=np.float32)
R_4 = np.array([5.0, 0.0, -wall_separation], dtype=np.float32)
R_5 = np.array([0.0, 0.0, -wall_separation], dtype=np.float32)
R_6 = np.array([-5.0, 0.0, -wall_separation], dtype=np.float32)

# Initial velocities
U_1 = np.array([1.0, 0.0, 0.0], dtype=np.float32)
U_2 = np.array([1.0, 0.0, 0.0], dtype=np.float32)
U_3 = np.array([1.0, 0.0, 0.0], dtype=np.float32)
U_4 = np.array([1.0, 0.0, -0.0], dtype=np.float32)
U_5 = np.array([1.0, 0.0, -0.0], dtype=np.float32)
U_6 = np.array([1.0, 0.0, -0.0], dtype=np.float32)

# Initial accelerations
A_1 = np.array([1.0, 0.0, 0.0], dtype=np.float32)
A_2 = np.array([1.0, 0.0, 0.0], dtype=np.float32)
A_3 = np.array([1.0, 0.0, 0.0], dtype=np.float32)
A_4 = np.array([1.0, 0.0, -0.0], dtype=np.float32)
A_5 = np.array([1.0, 0.0, -0.0], dtype=np.float32)
A_6 = np.array([1.0, 0.0, -0.0], dtype=np.float32)
# !SECTION (Parameters)


# SECTION: Snapshot initialization
s = gsd.hoomd.Snapshot()

s.configuration.step = step  # simulation time step
# number of spatial dimensions in simulation
s.configuration.dimensions = num_dimensions
# simulation box size (periodic)
#   box[0:3]: (𝑙𝑥,𝑙𝑦,𝑙𝑧)  the box length in each direction, in length units
#   box[3:]: (𝑥𝑦, 𝑥𝑧, 𝑦𝑧) the tilt factors, unitless values
s.configuration.box = [box_size, box_size, box_size, 0, 0, 0]

# Needed, as GSD does not have these saved by default for some reason?
s.log['configuration/step'] = step
s.log['configuration/dimensions'] = num_dimensions

s.log['integrator/dt'] = dt
s.log['integrator/t'] = t
s.log['integrator/tf'] = tf
s.log['integrator/tau'] = tau
s.log['integrator/num_steps_output'] = num_steps_output

s.log['material_parameters/fluid_density'] = fluid_density
s.log['material_parameters/particle_density'] = particle_density

s.log['wca/epsilon'] = wca_epsilon
s.log['wca/sigma'] = wca_sigma

s.particles.N = 6  # number of particles
s.particles.types = ['Squirmer', 'Image']  # number of particle types
s.particles.typeid = [0, 0, 0,  # Type-ID of each particle
                      1, 1, 1]
s.particles.diameter = [2*a, 2*a, 2*a,  # diameter
                        2*a, 2*a, 2*a]

s.particles.position = [  # initial positions
    R_1, R_2, R_3, R_4, R_5, R_6
]
s.particles.velocity = [  # initial velocities
    U_1, U_2, U_3, U_4, U_5, U_6
]
# NOTE[epic=data,seq=1]: Acceleration data is being stored in
# moment_inertia field as there is not a proper storage field for such data
s.particles.moment_inertia = [  # initial accelerations
    A_1, A_2, A_3, A_4, A_5, A_6
]
# !SECTION (Snapshot initialization)


# SECTION: Trajectory initialization
traj = gsd.hoomd.open(name=gsd_name, mode='wb')
traj.append(s)
# !SECTION: (Trajectory initialization)

# For use when being called from command line
if __name__ == "__main__":
    options, remainder = parser.parse_args(sys.argv[1:])

    RelPath = str(options.user_relPath)
    OutputDir = str(options.user_outputDir)

    aggregate_plots(RelPath, OutputDir)
