#include "domain.hpp"
#include "superellipsoid.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <tuple>

/***** Auxilirary functions *****/

/**
 * Checks p1 and p2 are colliding.
 *
 * @param `p1`: superellipsoid object
 * @param `p2`: superellipsoid object
 * @return 2 element tuple with first element being true if `p1` and 
 * 		   `p2` are colliding, otherwise false. The second element
 * 		   is the minimum distance between two points, one on each
 * 		   superellipsoid.
 */
/*
std::tuple<bool, double> check_collision(Superellipsoid* p1, Superellipsoid* p2) {
	Eigen::Vector3d c1 = p1->get_center();
	Eigen::Vector3d c2 = p2->get_center();
	double center_dist2 = (c1[0] - c2[0])*(c1[0] - c2[0]) + 
						  (c1[1] - c2[1])*(c1[1] - c2[1]) + 
		  			      (c1[2] - c2[2])*(c1[2] - c2[2]);
	
	double r1 = p1->inscribed_sphere_radius();
	double r2 = p2->inscribed_sphere_radius();
	double R1 = p1->circumscribed_sphere_radius();
	double R2 = p2->circumscribed_sphere_radius();

	std::tuple<bool,double> ret;
	if ( center_dist2 < r1*r1 + r2*r2 ) {
		// Inside inscribed sphere radius => collision 
		ret = {true, 0.0};
	} else if ( center_dist2 < R1*R1 + R2*R2 ) {
		// Between inscribed and circumscribed sphere so 
		// checks collision using exact distance
		std::tuple<bool,double> temp = Superellipsoid::distance(p1, p2);
		ret = {!std::get<0>(temp), std::get<1>(temp)};
	} else {
		// Outside circumscribed sphere so no collision so just returning
		// false and the exact distance.
		std::tuple<bool,double> temp = Superellipsoid::distance(p1, p2);
		ret = {false, std::get<1>(temp)};
	}
	return ret;
}
*/
/***** Private *****/

void Domain::remove_particle_af(Superellipsoid* p) { 
	// Searches for the index at which p is located in advancing_front
	int remove_idx = -1;
	for (int ix = 0; ix < this->advancing_front.size(); ix++) {
		if (this->advancing_front.at(ix) == p) remove_idx = ix;
	}

	// Only remove element if p was found in af
	if (remove_idx > -1) this->advancing_front.erase(this->advancing_front.begin() + remove_idx); 
}

void Domain::add_particle_af(Superellipsoid* p) {
	this->advancing_front.push_back(p);
}

std::vector<Superellipsoid*> Domain::particles_in_subdomain(double x_range[2], double y_range[2], double z_range[2]) {
	return this->particles.particles_in_domain(x_range, y_range, z_range);
}

void Domain::initialise_outward_advancing_front(Superellipsoid* particles[4]) {
	double x_mid = (this->x_bounds[0] + this->x_bounds[1])/2;
	double y_mid = (this->x_bounds[0] + this->x_bounds[1])/2;
	double z_mid = (this->x_bounds[0] + this->x_bounds[1])/2;
	Eigen::Vector3d domain_center(x_mid, y_mid, z_mid);
	
	// Setting up variables to make it easy to set centers of particles
	double r_max = 0;
	for (int ix = 0; ix < 4; ix++) {
		double r = particles[ix] -> circumscribed_sphere_radius();
		if (r_max < r) r_max = r; 
	}
	Eigen::Vector3d p1_local_center(r_max, 0, -r_max/std::sqrt(2));
	Eigen::Vector3d p2_local_center(-r_max, 0, -r_max/std::sqrt(2));
	Eigen::Vector3d p3_local_center(0, r_max, r_max/std::sqrt(2));
	Eigen::Vector3d p4_local_center(0, -r_max, r_max/std::sqrt(2));

	// Setting centers of particles and adds them to domain and af 
	particles[0]->set_center(p1_local_center + domain_center);
	this->add_particle(particles[0]);
	this->add_particle_af(particles[0]);
	particles[1]->set_center(p2_local_center + domain_center);
	this->add_particle(particles[1]);
	this->add_particle_af(particles[1]);
	particles[2]->set_center(p3_local_center + domain_center);
	this->add_particle(particles[2]);
	this->add_particle_af(particles[2]);
	particles[3]->set_center(p4_local_center + domain_center);
	this->add_particle(particles[3]);
	this->add_particle_af(particles[3]);

	// Todo: move the particles into close contact while keeping their shared center constant
	
}

/*
int Domain::binary_approach(Superellipsoid* fixed_particles[3], Superellipsoid* mobile_particle) { 
	Superellipsoid* fp1 = fixed_particles[0];
	Superellipsoid* fp2 = fixed_particles[1];
	Superellipsoid* fp3 = fixed_particles[2];

	Eigen::Vector3d c_mp = mobile_particle->get_center();
	double R1 = mobile_particle->circumscribed_sphere_radius();
	double r1 = mobile_particle->inscribed_sphere_radius();

	Eigen::Vector3d c_fp1 = fp1->get_center();
	Eigen::Vector3d c_fp2 = fp2->get_center();
	Eigen::Vector3d c_fp3 = fp3->get_center();
	
	Eigen::Vector3d mid_point = (c_fp1 + c_fp2 + c_fp3)/3;
	Eigen::Vector3d domain_mid_point((this->x_bounds[0] + this->x_bounds[1])/2,
									 (this->y_bounds[0] + this->y_bounds[1])/2,
									 (this->z_bounds[0] + this->z_bounds[1])/2);

	// Computing the normal vector to the plane intersecting
	// all three fixed particles' center 
	Eigen::Vector3d n_vec = (c_fp2 - c_fp1).cross(c_fp3 - c_fp1);
	n_vec.normalize();
	// Ensuring that n_vec is pointed somewhat in the direction away from the domain center 
	if ( (mid_point - domain_mid_point).norm() > (mid_point - domain_mid_point + 0.1*n_vec).norm() ) n_vec = -n_vec;

	double r_max = R1;
	for (int ix = 0; ix < 3; ix++) {
		double r = fixed_particles[ix]->circumscribed_sphere_radius();
		if (r > r_max) r_max = r;
	}
	Eigen::Vector3d starting_point = mid_point + r_max * n_vec;

	int relocation_counter = 0;
	double lambda = 0; // When this is 1 new_center is equal to mid_point
	double lambda_hi = 1;
	double lambda_lo = 0;
	while ( true ) {
		Eigen::Vector3d new_center = starting_point - lambda*r_max*n_vec;
		mobile_particle->set_center(new_center);
		relocation_counter++;

		// Domain which can include particles to colide with
		// TODO: this domain can be shrunk and do not need to 
		// be updated on every iteration
	 	double d = 2 * this->larges_circumscribing_sphere_radius;
		double x_bounds[2] = {new_center[0] - d, new_center[0] + d};
		double y_bounds[2] = {new_center[1] - d, new_center[1] + d};
		double z_bounds[2] = {new_center[2] - d, new_center[2] + d};
		std::vector<Superellipsoid*> possible_collisions = this->particles_in_subdomain(x_bounds, y_bounds, z_bounds);
	
		// Loops through possible collisions and checks for collisions
		bool collision = false; // Initialises to flase will make it cover the case when possible_collisions.size() = 0.
		for (int ix = 0; ix < possible_collisions.size(); ix++) {
			std::tuple<bool,double> coll_dist = check_collision(mobile_particle, possible_collisions[ix]);
			collision = std::get<0>(coll_dist);
			
			// Exits collision check loop on first collision
			if (collision) break;
		}

		// Computes distance between mobile particle and fixed particles
		double mobile_fixed_dist[3];
		for (int ix = 0; ix < 3; ix++) {
			std::tuple<bool,double> temp = Superellipsoid::distance(mobile_particle, fixed_particles[ix]);
			mobile_fixed_dist[ix] = std::get<1>(temp);
		}

		if (collision) {
			// If we have collision we take a step away from mid_point
			std::cout << "[INFO]: collision detected between mobile particle and at least one other particle." << std::endl;
			lambda_lo = lambda;
		} else if (mobile_fixed_dist[0] <= this->contact_tol || mobile_fixed_dist[1] <= this->contact_tol || mobile_fixed_dist[2] <= this->contact_tol) {
			// If mobile particle is close enough to any of the fixed particles we terminate
			break;
		} else {
			// If no collision nor close enough to fixed particle 
			// we take a step towards mid_point
			lambda_hi = lambda;
		}
		// Updates step length and direction 
		lambda = lambda_lo + (lambda_hi - lambda_lo)/2;
	}	
	
	std::cout << "[INFO]: relocated particle " << relocation_counter << " times." << std::endl;
	return 0;
}
*/
/***** Public *****/

Domain::Domain(double x_range[2], double y_range[2], double z_range[2], double contact_tol) {
	this->contact_tol = contact_tol;
	this->x_bounds[0] = x_range[0];
	this->x_bounds[1] = x_range[1];
	this->y_bounds[0] = y_range[0];
	this->y_bounds[1] = y_range[1];
	this->z_bounds[0] = z_range[1];
	this->z_bounds[1] = z_range[1];

	this->particles = CoordinateIndexer();
}

Domain::Domain(double x_range[2], double y_range[2], double z_range[2], double contact_tol, int n_particles) {
	this->contact_tol = contact_tol;
	this->x_bounds[0] = x_range[0];
	this->x_bounds[1] = x_range[1];
	this->y_bounds[0] = y_range[0];
	this->y_bounds[1] = y_range[1];
	this->z_bounds[0] = z_range[1];
	this->z_bounds[1] = z_range[1];

	this->particles = CoordinateIndexer(n_particles);
}

Domain::Domain(double x_range[2], double y_range[2], double z_range[2], double contact_tol, std::vector<Superellipsoid*>* particles) {
	this->contact_tol = contact_tol;
	this->x_bounds[0] = x_range[0];
	this->x_bounds[1] = x_range[1];
	this->y_bounds[0] = y_range[0];
	this->y_bounds[1] = y_range[1];
	this->z_bounds[0] = z_range[1];
	this->z_bounds[1] = z_range[1];

	this->particles = CoordinateIndexer(particles);
}

int Domain::n_particles() { return this->particles.n_particles(); }

void Domain::add_particle(Superellipsoid* p) {
	this->particles.add_particle(p);
	
	double r = p->circumscribed_sphere_radius();
	if (r > this->larges_circumscribing_sphere_radius) larges_circumscribing_sphere_radius = r; 
}

double Domain::volume(void) {
	double x_len = this->x_bounds[1] - this->x_bounds[0];
	double y_len = this->y_bounds[1] - this->y_bounds[0];
	double z_len = this->z_bounds[1] - this->z_bounds[0];

	return x_len * y_len * z_len;
}


