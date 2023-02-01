# partycle--
A C++ program for geometrically generating particle mixtures and computing contact statistics. Since Github does not support "+" in repo names, the name is decremented instead of incremented.

The program was developed as part of the course "Project course in mathematical and statistical modelling" at Chalmers University of Technology. The project was hosted by AstraZeneca and it resulted in a report which can be found in the repo (`report.pdf`).

## Getting started
### Building the project
To use the program it first needs to be built. To do these follow the below steps:
1. Ensure that [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page) is installed where `g++` can find it.
2. Clone the repo and create a new folder at the project root named `build`.
3. `cd` into `build/` and run `cmake ..`, this will generate the necessary Make files.
4. Run `mak all` in `build/`.

### Running the program
Afte the project has been succesfully built there will be an executable in `build/bin/` named `partycle--`. This can be executed using according to below:
- `partycle-- -d [DOMAIN DEFINITION] -cf [FILE] -ct [CONTACT TOLERANCE]`
- `partycle-- --domain [DOMAIN DEFINITION] --component-file [FILE] --contact-tolerance [CONTACT TOLERANCE]`

**Argument specification:**
- The domain definition follows the format `[ax,bx]x[ay,by]x[az,bz]`, where `a*` is the lower bound in respective axis while `b*` is the upper bound. Note that no white spaces are alowed and that the first brackets define x-space, second y-space, and the third z-space.
- The provided file need to be a `csv`-file without blank spaces. More about the format can be found under [Componenst list format](#Components-list-format).
- The default contact tolerance is set to 1e-2, this means that if `-ct`/`--contact-tolerance` is left out the program will use 1e-2 as contact tolerance.

When the program is finished a `csv`-file called `mixture.csv` has been created containing the mixture. This can be utilisd through the `plot.py` script to plot the mixture in 3d space. Note that this will take a LONG time for large mixtures. To change colours, plot axes, boundaries etc., just edit the constants in `plot.py` 
  

## Components list format
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

### Volume distributions
The volume distribution of a given component describes the volume of a particle of that component. This distribution does not describe the shape of the particle, this is done by the components specified reference particle. In other words; if a realisation of the volume distribution for some particle of some component yields 2, then the component's reference particle will have its scale parameters equally scaled so that the volume becomes 2.  

To specify the volume distribution for a component, a string on the form `<distribution_name>(<parameters...>)` is utilised in the components file. The available distributions and relating parameters is presented in the below table.
| **String format**    | **Parameters** | **Example** |
| -------------------- | -------------- | ----------- |
| uniform(a,b)         | `a` and `b` defines the interval [a,b] on which is the distributions support. | uniform(1,9.9) |
| normal(mu,sigma)     | `mu` is the expectation of a random variable with this distribution. `sigma` is the variance of a random variable with this distribution. | normal(2.5,0.25) |
| log-normal(mu,sigma) | `mu` is the expectation of the logarithm of a random variable with this distribution. `sigma` is the variance of the logarithm of a random variable with this distribution. | log-normal(1,0.25) |
| weibull(k, lambda)   | `k` is the shape parameter. `lambda` is the scale parameter. | weibull(0,0.25) |
