from typing import List, Tuple
from matplotlib.cbook import contiguous_regions

import numpy as np
import matplotlib as plt
import csv
import quaternion as qn

########### PACKING SETTINGS ###########
DOMAIN_CSV = "build/domain.csv"

class Superellipsoid:
    def __init__(self, component_id: int, scale_params: Tuple[float,float,float], shape_params: Tuple[float,float], rotation: qn.quaternion, center: np.ndarray):
        self.component_id = component_id
        self.scale_params = scale_params
        self.shape_param = shape_params
        self.rot = rotation
        self.center = center


def load_domain(file_path: str) -> List[Superellipsoid]:
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

            x = float(entries[5])
            y = float(entries[7])
            z = float(entries[8])

            rot_w = float(entries[9])
            rot_x = float(entries[10])
            rot_y = float(entries[11])
            rot_z = float(entries[12])

            particles.append(Superellipsoid(cid, (a,b,c), (n1,n2), qn.quaternion(rot_w, rot_x, rot_y, rot_z), np.array([x,y,z])))

    return particles
    
def main():
    domain = load_domain(DOMAIN_CSV)
    print(domain)

if __name__ == "__main__":
    main()
