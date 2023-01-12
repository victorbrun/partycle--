#include "domain.hpp"
#include "superellipsoid.hpp"
#include "matplotlibcpp.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <tuple>
#include <random>

bool check_collision(Superellipsoid* p1, Superellipsoid* p2) {
	Eigen::Vector3d c1 = p1->get_center();
	Eigen::Vector3d c2 = p2->get_center();
	double center_dist2 = (c1[0] - c2[0])*(c1[0] - c2[0]) + 
						  (c1[1] - c2[1])*(c1[1] - c2[1]) + 
		  			      (c1[2] - c2[2])*(c1[2] - c2[2]);
	
	double r1 = p1->inscribed_sphere_radius();
	double r2 = p2->inscribed_sphere_radius();
	double R1 = p1->circumscribed_sphere_radius();
	double R2 = p2->circumscribed_sphere_radius();

	bool ret;
	if ( center_dist2 < r1*r1 + r2*r2 ) {
		// Inside inscribed sphere radius => collision 
		ret = true;
	} else if ( center_dist2 < R1*R1 + R2*R2 ) {
		// Between inscribed and circumscribed sphere so 
		// checks collision using exact distance
		ret = Superellipsoid::distance(p1, p2);
	} else {
		// Outside circumscribed sphere so no collision so just returning
		// false and the exact distance.
		ret = false;
	}
	return ret;
}

/***** Private *****/

void Domain::remove_particle_af(Superellipsoid* p) { 
	// Searches for the index at which p is located in advancing_front
	int remove_idx = -1;
	int sz = this->advancing_front.size();
	for (int ix = 0; ix < sz; ix++) {
		if (this->advancing_front.at(ix) == p) remove_idx = ix;
	}

	// Only remove element if p was found in af
	if (remove_idx > -1) this->advancing_front.erase(this->advancing_front.begin() + remove_idx); 
}

void Domain::add_particle_af(Superellipsoid* p) {
	this->advancing_front.push_back(p);
}

std::vector<Superellipsoid*> Domain::particles_in_subdomain(double x_range[2], double y_range[2], double z_range[2]) {
	return this->particles->particles_in_domain(x_range, y_range, z_range);
}

void Domain::increment_advancing_front(Superellipsoid* p) {
	std::vector<Superellipsoid*> af = this->advancing_front;
	int sz = af.size();
	if (sz == 0) {
		return; 
	}

	// Draw a reference particle
	std::random_device rd; 
    std::mt19937 mt(rd()); 
    std::uniform_int_distribution<> uni(0, sz);
	int idx = uni(mt);
	Superellipsoid* p_ref = af[idx];
	
	
	// Get number of contacts for p_ref
	int no_nonzeros = p_ref->get_contacts().size();

	bool all_fails = true;  // Default to true, update if success
	bool cc; // Tracks collisions
	// Test all combinations of neighbour pairs, given that we have at least 2 contacts
	// Running binary approach. If successfull we do not try any other pair.
	//  If unsuccessfull we continue to next pair of particles. If no pair of particles 
	//  are successfull we deem the reference point unreachable and delete it from the 
	//  advancing front.
	if (no_nonzeros >= 2) {
		for (int ix = 0;  ix < no_nonzeros-1; ix++) {
			Superellipsoid* p1 = af.at(ix);
			for (int jx = ix + 1; jx < no_nonzeros; jx++) {
				Superellipsoid* p2 = af.at(jx);
				
				std::vector<Superellipsoid*> particles; 
				particles.push_back(p_ref);
				particles.push_back(p1);
				particles.push_back(p2);

				// Try to approach chosen particles
				cc = this->binary_approach(particles, p);
				if (cc) {
					// Binary approach succeeded:
					// Add the particle to the advancing front and domain list
					this->add_particle_af(p);
					this->add_particle(p);

					//Update contacts
					p->add_contact(particles[0]);
					p->add_contact(particles[1]);
					p->add_contact(particles[2]);
					
					particles[0]->add_contact(p);
					particles[1]->add_contact(p);
					particles[2]->add_contact(p);

					//
					all_fails = false; 
					break;
				} else {
					// Binary approach failed:  
					//pair[0].remove_contact(p_ref)
                	//pair[1].remove_contact(p_ref)

				}

			}
		}
	}
	if (all_fails) {
		// All possible particle assemblies failed: 
		this->remove_particle_af(p);

	}
	return; 
}

int Domain::binary_approach(std::vector<Superellipsoid*> fixed_particles, Superellipsoid* mobile_particle) { 
	Superellipsoid* fp1 = fixed_particles.at(0);
	Superellipsoid* fp2 = fixed_particles.at(1);
	Superellipsoid* fp3 = fixed_particles.at(2);

	double R1 = mobile_particle->circumscribed_sphere_radius();
	//double r1 = mobile_particle->inscribed_sphere_radius();

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

	// Get the largest radius of ref particles
	double r_max = R1;
	for (int ix = 0; ix < 3; ix++) {
		double r = fixed_particles[ix]->circumscribed_sphere_radius();
		if (r > r_max) r_max = r;
	}
	
	Eigen::Vector3d starting_point = mid_point + (r_max+0.1) * n_vec;

	// Return value; true indicates success
	bool ret = true; 

	// Counters
	int iterations = 0; 
	int collision_counter = 0;
	int relocation_counter = 0;

	// Step parameters
	double lambda = 0;    // lambda ranges from 0 to 1, moving the particle from the start to end point
	double lambda_hi = 1; // lambda and hi/lo are updated depending on collision check outcomes
	double lambda_lo = 0;

	// Misc
	Eigen::Vector3d p_vec = mid_point-starting_point;
	bool collision_old;
	double dlambda; 
	// Tolerances & termination parameters
	double tol1 = 1e-2;	   // Overlap tolerance
	double tol2 = 1e-2;    // Within-distance tolerance
	int N = 20; 		   // Max iteration count
	int K = 5;             // Max sequential collision count

	// Termination flags
	bool F1 = false;	  // Success: stepsize < tol1 && no cc to cc : collision with small (tolerable) overlap
	bool F2 = false;	  // Success: stepsize < tol2 && cc to no cc : no collision with small distance between particles
	bool F3 = false;	  // Failure: i > K                          : Maximum number of iterations reached
	bool F4 = false;	  // Failure: i == 0 && cc == true           : Non - viable startpoint
	bool F5 = false;	  // Failure: collision_counter > K          :  Maximum number of sequential collisions reached

	double d = 2 * this->larges_circumscribing_sphere_radius;
	while (true) {
		Eigen::Vector3d new_center = starting_point + lambda*p_vec;
		mobile_particle->set_center(new_center);
		relocation_counter++;

		// Domain which can include particles to colide with
		// TODO: this domain can be shrunk and do not need to 
		// be updated on every iteration
		double x_bounds[2] = {new_center[0] - d, new_center[0] + d};
		double y_bounds[2] = {new_center[1] - d, new_center[1] + d};
		double z_bounds[2] = {new_center[2] - d, new_center[2] + d};
		std::vector<Superellipsoid*> possible_collisions = this->particles_in_subdomain(x_bounds, y_bounds, z_bounds);
	
		// Loops through possible collisions and checks for collisions
		bool collision = false; // Initialization to false will make it cover the case when possible_collisions.size() = 0.
		for (int ix = 0; ix < int(possible_collisions.size()); ix++) {
			collision = check_collision(mobile_particle, possible_collisions[ix]);
			// Exits collision check loop on first collision
			if (collision) break;
		}

		if (collision) {
			// If we have collision we take a step away from mid_point
			std::cout << "[INFO]: collision detected between mobile particle and at least one other particle." << std::endl;
			lambda_lo = lambda;
			collision_counter ++; 
		} else {
			// If no collision nor close enough to fixed particle 
			// we take a step towards mid_point and 
			lambda_hi = lambda;
			collision_counter = 0;
		}

		// Check termination criteria
		if (iterations > 0) {
			F1 = (d < tol1) && (collision_old = false) && (collision = true);
			F2 = (d < tol2) && (collision_old = true) && (collision = false);
		}
		F3 = (iterations > N);
		F4 = (iterations == 0) && collision;
		F5 = collision_counter > K;

		if (F1) {
			// Success
			ret = true;
			break; 
		} else if(F2){
			// Success
			ret = true;
			break; 
		} else if(F3){
			// Failure
			ret = false;
			break; 
		} else if(F4){
			// Failure
			ret = false;
			break; 
		} else if(F5){
			// Failure
			ret = false;
			break; 
		}

		// Update step size
		dlambda = lambda - (lambda_hi - lambda_lo)/2;
	
		// Updates step length and direction 
		lambda = lambda_lo + (lambda_hi - lambda_lo)/2;

		// Set old collision counter for the sake of exit flags
		collision_old = collision;

		// Update iterations
		iterations++; 
	}	
	
	std::cout << "[INFO]: relocated particle " << relocation_counter << " times." << std::endl;
	return ret;
}

/***** Public *****/

Domain::Domain(double x_range[2], double y_range[2], double z_range[2], double contact_tol) {
	this->contact_tol = contact_tol;
	this->x_bounds[0] = x_range[0];
	this->x_bounds[1] = x_range[1];
	this->y_bounds[0] = y_range[0];
	this->y_bounds[1] = y_range[1];
	this->z_bounds[0] = z_range[0];
	this->z_bounds[1] = z_range[1];

	this->particles = new CoordinateIndexer();
}

Domain::Domain(double x_range[2], double y_range[2], double z_range[2], double contact_tol, int n_particles) {
	this->contact_tol = contact_tol;
	this->x_bounds[0] = x_range[0];
	this->x_bounds[1] = x_range[1];
	this->y_bounds[0] = y_range[0];
	this->y_bounds[1] = y_range[1];
	this->z_bounds[0] = z_range[0];
	this->z_bounds[1] = z_range[1];

	this->particles = new CoordinateIndexer(n_particles);
}

Domain::Domain(double x_range[2], double y_range[2], double z_range[2], double contact_tol, std::vector<Superellipsoid*>* particles) {
	this->contact_tol = contact_tol;
	this->x_bounds[0] = x_range[0];
	this->x_bounds[1] = x_range[1];
	this->y_bounds[0] = y_range[0];
	this->y_bounds[1] = y_range[1];
	this->z_bounds[0] = z_range[0];
	this->z_bounds[1] = z_range[1];

	this->particles = new CoordinateIndexer(particles);
}

Domain::~Domain() {
	std::cout << "[INFO]: destroying domain instance" << std::endl;
	delete particles;
}
	
int Domain::n_particles() { return this->particles->n_particles(); }

std::vector<Superellipsoid*>* Domain::get_particles() { return this->particles->get_particles(); }

void Domain::add_particle(Superellipsoid* p) {
	this->particles->add_particle(p);
	
	double r = p->circumscribed_sphere_radius();
	if (r > this->larges_circumscribing_sphere_radius) larges_circumscribing_sphere_radius = r; 
}

double Domain::volume(void) {
	double x_len = this->x_bounds[1] - this->x_bounds[0];
	double y_len = this->y_bounds[1] - this->y_bounds[0];
	double z_len = this->z_bounds[1] - this->z_bounds[0];

	return x_len * y_len * z_len;
}

void Domain::initialise_outward_advancing_front(Superellipsoid* particles[4]) {
	double x_mid = (this->x_bounds[0] + this->x_bounds[1])/2;
	double y_mid = (this->x_bounds[0] + this->x_bounds[1])/2;
	double z_mid = (this->x_bounds[0] + this->x_bounds[1])/2;
	Eigen::Vector3d domain_center(x_mid, y_mid, z_mid);
	
	// Setting up variables to make it easy to set centers of particles
	// so that they are placed in the corners of a tetrahedron.
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

	// Add the particles to each others contact lists
	particles[0]->add_contact(particles[1]);
	particles[0]->add_contact(particles[2]);
	particles[0]->add_contact(particles[3]);

	particles[1]->add_contact(particles[0]);
	particles[1]->add_contact(particles[2]);
	particles[1]->add_contact(particles[3]);

	particles[2]->add_contact(particles[0]);
	particles[2]->add_contact(particles[1]);
	particles[2]->add_contact(particles[3]);

	particles[3]->add_contact(particles[0]);
	particles[3]->add_contact(particles[1]);
	particles[3]->add_contact(particles[2]);
	
	// All of the particles closer and closer to shared centrum until at least one pair is 
	// touching. This could be done faster by utilising binary approach but it does not 
	// really fit the written binary approach.
	double center_distance = r_max * std::sqrt(3/2); // Distance for each partile to center of tetrahedron
	double step_size = 0.01 * r_max;
	Eigen::Vector3d c1 = particles[0]->get_center();
	Eigen::Vector3d c2 = particles[1]->get_center();
	Eigen::Vector3d c3 = particles[2]->get_center();
	Eigen::Vector3d c4 = particles[3]->get_center();
	for (double offset = 0; offset < center_distance; offset += step_size) {

		// Decrements each particles distance to the shared centrum.
		particles[0]->set_center(c1 - p1_local_center*offset);
		particles[1]->set_center(c2 - p2_local_center*offset);
		particles[2]->set_center(c3 - p3_local_center*offset);
		particles[3]->set_center(c4 - p4_local_center*offset);

		// Checks if any of the particles are colliding,
		// if so we take a step back and break the loop,
		// otherwise we keep moving the particles forward.
		bool collision = check_collision(particles[0], particles[1]) || check_collision(particles[0], particles[2]) ||
						 check_collision(particles[0], particles[3]) || check_collision(particles[1], particles[2]) ||
						 check_collision(particles[1], particles[3]) || check_collision(particles[2], particles[3]);
		if (collision) {
			// Steps the particles back one step since we know that they do not touch there
			particles[0]->set_center(c1 - p1_local_center * (offset-step_size));
			particles[1]->set_center(c2 - p2_local_center * (offset-step_size));
			particles[2]->set_center(c3 - p3_local_center * (offset-step_size));
			particles[3]->set_center(c4 - p4_local_center * (offset-step_size));
			break;
		}
	}
}

void Domain::draw(int resolution){
	//Generate linspace for eta and omega
	Eigen::VectorXd eta, omega; 
	eta.setLinSpaced(resolution, -M_PI/2, M_PI/2);
	omega.setLinSpaced(resolution, -M_PI, M_PI);

	//Generate mesh
	Eigen::MatrixXd ETA = eta.rowwise().replicate(resolution).transpose();
	Eigen::MatrixXd OMEGA = omega.rowwise().replicate(resolution);

	//Generate surface plot points
	std::vector<double> f(resolution);
	std::vector<std::vector<double>> X(resolution, f), Y(resolution, f), Z(resolution, f);
	
	//Collect all particles
	std::vector<Superellipsoid*>* particles = this->get_particles();

	//Plotting configs
	namespace plt = matplotlibcpp;
	plt::figure(1);

	//Cmap args: 'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds'
	//Todo: plot surface based on not yet implemented class for particle
	std::map<std::string, std::string> keywords = {
		{"cmap", "Greens"}
	};

	//Iterate and plot
	for (int ix = 0; ix < this->n_particles(); ix++){
		//Collect mesh points and surfaceplot		
		std::tie(X, Y, Z) = particles->at(ix)->parametric_surface(ETA, OMEGA);
		std::cout << "here" << std::endl;
		plt::plot_surface(X, Y, Z, keywords, 1);
	}
	
	plt::show();
	//plt::save("images.png");
}
