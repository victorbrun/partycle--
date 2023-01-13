from typing import List, Tuple

import numpy as np
import matplotlib.pyplot as plt
import csv
import quaternion as qn

########### PACKING SETTINGS ###########
DOMAIN_CSV = "build/domain.csv"
DOMAIN_X_BOUNDS = (0, 10)
DOMAIN_Y_BOUNDS = (0, 10)
DOMAIN_Z_BOUNDS = (0, 10)
PARTICLE_PLOTTING_RESOLUTION = 20
COMPONENT_ID_TO_COLOUR = {
        1: "b", # blue 
        2: "g", # green
        3: "r", # red 
        4: "c", # cyan
        5: "m", # purple
        6: "y", # yellow
        7: "k", # black
}

class Superellipsoid:
    def __init__(self, component_id: int, scale_params: Tuple[float,float,float], shape_params: Tuple[float,float], rotation: qn.quaternion, center: np.ndarray):
        self.component_id = component_id
        self.scale_params = scale_params
        self.shape_param = shape_params
        self.rot = rotation
        self.center = center

    def surface(self, eta: np.ndarray, omega: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Generate meshgrid over the superellipsiod surface
        Returns:
            tuple[np.ndarray, np.ndarray, np.ndarray]: Tuple (x,y,z) of grid coordinates for each component
        """
        eta, omega = np.meshgrid(eta, omega)
            
        #Compute the surface in local coords
        a, b, c = self.scale_params[0], self.scale_params[1], self.scale_params[2]
        n1, n2 = self.shape_param[0], self.shape_param[1]
    
        #Surface vector
        x = a*np.sign(np.cos(eta))*np.abs(np.cos(eta))**(2/n1) * np.sign(np.cos(omega))*np.abs(np.cos(omega))**(2/n2) 
        y = b*np.sign(np.cos(eta))*np.abs(np.cos(eta))**(2/n1) * np.sign(np.sin(omega))*np.abs(np.sin(omega))**(2/n2) 
        z = c*np.sign(np.sin(eta))*np.abs(np.sin(eta))**(2/n1)
        
        # Rotating the meshgrid
        R = self.rot
        for ix in range(x.shape[0]):
            for jx in range(x.shape[1]):
                vec = np.array([x[ix,jx], y[ix,jx], z[ix,jx]])
                vec = qn.rotate_vectors(R, vec)

                # Updates old coords to rotated coords
                x[ix, jx] = vec[0]
                y[ix, jx] = vec[1]
                z[ix, jx] = vec[2]

        # Translating particle into its actual position
        x = x + self.center[0]
        y = y + self.center[1]
        z = z + self.center[2]
            
        return (x, y, z)

def load_particles(file_path: str) -> List[Superellipsoid]:
    particles = []

    with open(file_path, "r") as file:
        reader = csv.reader(file)
        
        header = True
        for row in reader:
            entries = row[0].split(";")
            
            # Skipping header
            if header:
                header = False
                continue

            cid = int(entries[0])
            
            a = float(entries[1])
            b = float(entries[2])
            c = float(entries[3])
            n1 = float(entries[4])
            n2 = float(entries[5])

            x = float(entries[6])
            y = float(entries[7])
            z = float(entries[8])

            rot_w = float(entries[9])
            rot_x = float(entries[10])
            rot_y = float(entries[11])
            rot_z = float(entries[12])

            particles.append(Superellipsoid(cid, (a,b,c), (n1,n2), qn.quaternion(rot_w, rot_x, rot_y, rot_z), np.array([x,y,z])))

    return particles

def draw_particles(particles: List[Superellipsoid], resolution: int) -> None:
    """
    Generates pyplot of domain with particles in it. The resolution 
    specifies resolution of the parametrisation of each particle.
    """
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.axes.set_xlim3d(left=DOMAIN_X_BOUNDS[0], right=DOMAIN_X_BOUNDS[1]) 
    ax.axes.set_ylim3d(bottom=DOMAIN_Y_BOUNDS[0], top=DOMAIN_Y_BOUNDS[1])
    ax.axes.set_zlim3d(bottom=DOMAIN_Z_BOUNDS[0], top=DOMAIN_Z_BOUNDS[1]) 

    eta = np.linspace(-np.pi/2, np.pi/2, resolution)
    omega = np.linspace(-np.pi, np.pi, resolution)
    for particle in particles:
        x,y,z = particle.surface(eta, omega)
        ax.plot_surface(x, y, z, color=COMPONENT_ID_TO_COLOUR[particle.component_id])

    plt.show()

def main():
    domain = load_particles(DOMAIN_CSV)
    draw_particles(domain, PARTICLE_PLOTTING_RESOLUTION)

if __name__ == "__main__":
    main()
