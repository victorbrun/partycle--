# partycle--
C++ program for geometrically generating particle mixtures and computing contact statistics. Since Github does not support "+" in repo names, the name is decremented instead of incremented.

## Getting started


### Components list format
In order for the program to generate a mixture it must know which components to include in the mixture, this is specified by providing a ";"-separated CSV file. The fields in this CSV file is specified in below table.

| **Column name**       | **Description** | **Example**        |
| --------------------- | --------------- | ------------------ |
| component_id          | Unique integer. Used to distinguish between components.         | 1        |
| volume_distribution   | String representing a distribution and its parameters which describes the volume of the components reference particle. Available distributions are discussed below.	| normal(10,1) |
| volume_fraction       | Fraction describing how large part of the mixture that is composed of this component. This volume fraction regards the total volume of all particles and does not include the domain volume in which they are put, i.e. porosity is not included here (or it is more precisely assumed to be 0). This means that the sum of a mixture's components respective volume fraction must add up to 1. This quantity can be computed using the weight fraction of the components together with their respective densities. | 0.789	             |
| reference_particle_a  | Scale parameter in x-axis (before possible rotation).   	  | 1.5      |
| reference_particle_b  | Scale parameter in y-axis (before possible rotation).		  | 10       |
| reference_particle_c  | Scale parameter in z-axis (before possible rotation).		  | 0.123345 |
| reference_particle_n1 | Shape parameter in xy-plane (before possible rotation). A value of 2 makes the shape circular in xy-plane while a higher value makes it more cubic, but with smooth edges. The parameter must lie in the interval [2, 8].| 2.9      |
| reference_particle_n2 | Shape parameter in z-axis (before possible rotation). A value of 2 makes the shape circular in z-plane while a higher value makes it more cubic, but with smooth edges. The parameter must lie in the interval [2, 8].      | 8        |

An example of a components file in accordance with above format follows below.
```
component_id;volume_distribution;volume_fraction;reference_particle_a;reference_particle_b;reference_particle_c;reference_particle_n1;reference_particle_n2
1;uniform(1,10);0.25;1;1;1;2;2
2;log-normal(1,0.25);0.25;1;1;1;8;8
3;normal(10,2);0.25;5;1;1;5;2
4;weibull(1,5);0.25;1;3;1;2;2
```
