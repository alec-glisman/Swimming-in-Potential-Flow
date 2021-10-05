# @AUTHOR: alec-glisman (GitHub)

# SECTION: Dependencies
import numpy as np                 # mathematical data structures
import gsd.hoomd                   # Simulation data-structure
from pathlib import Path           # Modify path variables
# !SECTION


class GSDUtil:

    # SECTION: Static parameters

    # !SECTION (Static parameters)

    # SECTION: Initializer
    def __init__(self, gsd_path, create_gsd=True):
        self.newGSD = create_gsd

        self.gsdPath = gsd_path
        self.gsdName = Path(gsd_path).name

        self.trajectory = None

        if create_gsd:
            self.trajectory = gsd.hoomd.open(name=self.gsdPath, mode='wb')
            self.snapshot = gsd.hoomd.Snapshot()
        else:
            self.trajectory = gsd.hoomd.open(name=self.gsdPath, mode='rb')
            # load last snapshot in GSD file
            self.snapshot = self.trajectory.read_frame(
                self.trajectory.file.nframes - 1)

        # box parameters (not currently used)
        #   box[0:3]: (l_x, l_y, l_z)  the box length in each direction, in length units
        #   box[3:]: (xy, xz, yz) the tilt factors, unitless values
        self.snapshot.configuration.box = [30, 30, 30, 0, 0, 0]
    # !SECTION (Initializer)

    def setLogParameters(self,
                         dt, t, tf, tau,
                         num_steps_output=1000,
                         fluid_density=1.0, particle_density=1.0,
                         wca_epsilon=0.0, wca_sigma=0.0):
        # Convert data types to GSD expected type
        dt = np.array([dt], dtype=np.double)
        t = np.array([t], dtype=np.double)
        tf = np.array([tf], dtype=np.double)
        tau = np.array([tau], dtype=np.double)
        num_steps_output = np.array([num_steps_output], dtype=np.uint64)

        fluid_density = np.array([fluid_density], dtype=np.double)
        particle_density = np.array([particle_density], dtype=np.double)

        wca_epsilon = np.array([wca_epsilon], dtype=np.double)
        wca_sigma = np.array([wca_sigma], dtype=np.double)

        # Save data
        self.snapshot.log['integrator/dt'] = dt
        self.snapshot.log['integrator/t'] = t
        self.snapshot.log['integrator/tf'] = tf
        self.snapshot.log['integrator/tau'] = tau
        self.snapshot.log['integrator/num_steps_output'] = num_steps_output

        self.snapshot.log['material_parameters/fluid_density'] = fluid_density
        self.snapshot.log['material_parameters/particle_density'] = particle_density

        self.snapshot.log['wca/epsilon'] = wca_epsilon
        self.snapshot.log['wca/sigma'] = wca_sigma

    def setParticleParameters(self, N, types=None, typeid=None, diameter=None):
        if self.newGSD:
            step = 0
        else:
            step = len(self.trajectory)
        step = np.array([step], dtype=np.uint64)

        dimensions = np.array([3], dtype=np.uint8)

        self.snapshot.configuration.step = step  # simulation time step
        # number of spatial dimensions in simulation
        self.snapshot.configuration.dimensions = dimensions
        # simulation box size (periodic)
        # Needed, as GSD does not have these saved by default for some reason?
        self.snapshot.log['configuration/step'] = step
        self.snapshot.log['configuration/dimensions'] = dimensions

        self.snapshot.particles.N = int(N)  # number of particles

        if types is not None:
            self.snapshot.particles.types = types  # number of particle types
        if typeid is not None:
            self.snapshot.particles.typeid = typeid  # Type-ID of each particle

        if diameter is not None:
            self.snapshot.particles.diameter = diameter
        else:
            self.snapshot.particles.diameter = [2] * int(N)

    def setKinematics(self, Q, X, U, A):
        # Convert data types to personal expected type
        Q = np.array(Q, dtype=np.double)
        X = np.array(X, dtype=np.double)
        U = np.array(U, dtype=np.double)
        A = np.array(A, dtype=np.double)

        # Save kinematics
        self.snapshot.log['particles/double_orientation'] = Q
        self.snapshot.log['particles/double_position'] = X
        self.snapshot.log['particles/double_velocity'] = U
        self.snapshot.log['particles/double_moment_inertia'] = A

        # Convert data types to GSD expected type
        Q = np.array(Q, dtype=np.single)
        X = np.array(X, dtype=np.single)
        U = np.array(U, dtype=np.single)
        A = np.array(A, dtype=np.single)

        # Save kinematics
        self.snapshot.particles.orientation = Q
        self.snapshot.particles.position = X
        self.snapshot.particles.velocity = U
        self.snapshot.particles.moment_inertia = A

    def saveSnapshot(self):
        self.trajectory.append(self.snapshot)
