#include "domain.hpp"
#include "superellipsoid.hpp"

#include <algorithm>
#include <cmath>
#include <eigen3/Eigen/src/Geometry/Quaternion.h>
#include <fstream>
#include <iostream>
#include <random>
#include <regex>
#include <stdexcept>


/***** Private *****/

bool Domain::check_collision(Superellipsoid* p1, Superellipsoid* p2) {
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

bool Domain::check_collision(Superellipsoid* p, const std::vector<Superellipsoid*>& p_vec) {
	for (size_t ix = 0; ix < p_vec.size(); ix++) {
		if (Domain::check_collision(p, p_vec.at(ix))) return true;
	}
	return false;
}

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
	// TODO: we must choose new reference particle until we cannot chose another
	// currently we only choose one and if it does not work, i.e. should be removed from af 
	// we just terminate this function without doing anything
	std::random_device rd; 
	std::mt19937 mt(rd()); 

	while (this->advancing_front.size() > 0) {
		std::vector<Superellipsoid*> af = this->advancing_front;

		// Draw a reference particle
		std::uniform_int_distribution<> uni(0, af.size()-1);
		int idx = uni(mt);
		Superellipsoid* p_ref = af.at(idx);

		// Extracts reference particles contacts and checks so that there are enough 
		// of them to preform the binary approach algorithm.
		std::vector<Superellipsoid*> ref_contacts = p_ref->get_contacts();
		if (ref_contacts.size() < 2) {
			throw std::runtime_error("reference particle does not have enough contacts. Must be at least 2.");
		}

		int collision_result;
		for (size_t ix = 0; ix < ref_contacts.size()-1; ix++) {
			Superellipsoid* p1 = ref_contacts.at(ix);
			for (size_t jx = ix+1; jx < ref_contacts.size(); jx++) {
				Superellipsoid* p2 = ref_contacts.at(jx);
				std::vector<Superellipsoid*> fixed_ps = {p_ref, p1, p2};

				collision_result = this->binary_approach(fixed_ps, p);
				if (collision_result == 0 || collision_result == 1 || collision_result == 5) {
					// Binary approach suceeded so we add particle to domain and 
					// advancing front.
					this->add_particle(p);
					this->add_particle_af(p);

					// Update contacts
					// TODO: contact should only be added for those particles 
					// which are within contact distance 
					p->add_contact(fixed_ps.at(0));
					p->add_contact(fixed_ps.at(1));
					p->add_contact(fixed_ps.at(2));
					fixed_ps.at(0)->add_contact(p);
					fixed_ps.at(1)->add_contact(p);
					fixed_ps.at(2)->add_contact(p);

					// We will permamently place the particle at the first 
					// possible possition so since binary approach succeeded
					// we return here.
					return;
				}

				// If arrive here and it is the last iterartion of both loops
				// then the reference particle is such that no new particle can be placed 
				// in contact with it and 2 neighbours, so we remove it from the advancing front 
				if (ix == ref_contacts.size()-2 && jx == ref_contacts.size()-1) this->remove_particle_af(p_ref);
			}
		}
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
	Eigen::Vector3d domain_mid_point( this->x_bounds[0] + (this->x_bounds[1] - this->x_bounds[0])/2,
									  this->y_bounds[0] + (this->y_bounds[1] - this->y_bounds[0])/2,
									  this->z_bounds[0] + (this->z_bounds[1] - this->z_bounds[0])/2);

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
	
	// First point mobile particle is placed in and vector we move it along
	Eigen::Vector3d starting_point = mid_point + (2*r_max) * n_vec;
	Eigen::Vector3d p_vec = mid_point-starting_point;
	
	// Return value; 0 indicates success
	int exit_code = 0; 

	// Counters
	int iterations = 0; 
	int collision_counter = 0;
	int relocation_counter = 0;

	// Cache varables
	bool collision_old;
	Eigen::Vector3d old_center;

	// Tolerances & termination parameters
	double tol1 = this->contact_tol;	// Overlap tolerance
	double tol2 = this->contact_tol;	// Within-distance tolerance
	int N = 20; 		   				// Max iteration count
	int K = 10;             			// Max sequential collision count

	// Termination flags
	bool F1 = false;	  // Success: stepsize < tol1 && no cc to cc : collision with small (tolerable) overlap. 			Exit code: 0
	bool F2 = false;	  // Success: stepsize < tol2 && cc to no cc : no collision with small distance between particles.	Exit code: 1
	bool F3 = false;	  // Failure: i > N                          : Maximum number of iterations reached.				Exit code: 2					
	bool F4 = false;	  // Failure: i == 0 && cc == true           : Non - viable startpoint.								Exit code: 3
	bool F5 = false;	  // Failure: collision_counter > K          : Maximum number of sequential collisions reached.		Exit code: 4

	// Checks that we have collision in mid_point and no collision in starting_point,
	// i.e. checks that there must be a possition in close contact in (mid_point, starting_point]
	std::vector<Superellipsoid*>* all_particles = this->get_particles();
	mobile_particle->set_center(mid_point);
	if (Domain::check_collision(mobile_particle, *all_particles) == false) {
		return 5;
	}
	mobile_particle->set_center(starting_point);
	if (Domain::check_collision(mobile_particle, *all_particles) == true) {
		// Particle cannot be in conatct with other particles at starting point
		return 3;
	}

	// Step parameters. These are updated depending on collision check outcomes
	double lambda_hi = 1; // lambda and hi/lo are updated depending on collision check outcomes
	double lambda_lo = 0;
	double d = 2 * this->larges_circumscribing_sphere_radius;
	while (true) {
		double lambda = lambda_lo + (lambda_hi - lambda_lo)/2;
		
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
		//std::vector<Superellipsoid*>* possible_collisions = this->get_particles();
	
		// Updates the domain in which to search for possible placement 
		// depending on if current placement has collision or not
		bool collision = Domain::check_collision(mobile_particle, possible_collisions);
		if (collision) {
			collision_counter++;
			lambda_hi = lambda - this->contact_tol;
		} else {
			collision_counter = 0;
			lambda_lo = lambda + this->contact_tol;
		}
	
		// Check termination criteria
		if (iterations > 0) {
			double step_size = (new_center - old_center).norm();
			F1 = (step_size < tol1) && (collision_old == false) && (collision == true);
			F2 = (step_size < tol2) && (collision_old == true) && (collision == false);
		}
		F3 = (iterations > N);
		F4 = (iterations == 0) && collision;
		F5 = collision_counter > K;

		if (F1) {
			// Success
			exit_code = 0;
			break;
		} else if(F2){
			// Success
			exit_code = 1;
			break;
		} else if(F3){
			// Failure
			exit_code = 2;
			break;
		} else if(F4){
			// Failure
			exit_code = 3;
			break;
		} else if(F5){
			// Failure
			exit_code = 4;
			break;
		}

		// Updates stuff needed for next iteration
		old_center = new_center; 
		collision_old = collision;
		iterations++;
	}

	//std::cout << "[INFO]: relocated particle " << relocation_counter << " times." << std::endl;
	return exit_code;
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

Domain::Domain(double x_range[2], double y_range[2], double z_range[2], double contact_tol, size_t n_particles) {
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
	delete particles;
}
	
size_t Domain::n_particles() { return this->particles->n_particles(); }

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
		double r = particles[ix]->circumscribed_sphere_radius();
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

void Domain::write_csv(const std::string& file_name) {
	// Checking that the component-file is a csv file 
	std::regex file_format1(".*\\.csv$");
	std::regex file_format2(".*\\.CSV$");
	if ( !std::regex_match(file_name, file_format1) && !std::regex_match(file_name, file_format2) ) {
		throw std::runtime_error("invalid file format. Must be CSV");
	}

	std::ofstream file;
	std::string header = "component_id;a;b;c;n1;n2;center_x;center_y;center_z;rot_w;rot_x;rot_y;rot_z";
	std::vector<Superellipsoid*>* particles = this->get_particles();

	file.open(file_name);
	file << header << "\n";
	for (size_t ix = 0; ix < this->n_particles(); ix++) {
		Superellipsoid* p = particles->at(ix);
		int cid = p->get_component_id();
		double a = p->get_scale("a");
		double c = p->get_scale("b");
		double b = p->get_scale("c");
		double n1 = p->get_shape("n1");
		double n2 = p->get_shape("n2");
		double x = p->get_center()[0];
		double y = p->get_center()[1];
		double z = p->get_center()[2];
		
		// The quaternion is on the form w + xi + yj + zk
		double rot_w = p->get_orientation().w();
		double rot_x = p->get_orientation().x();
		double rot_y = p->get_orientation().y();
		double rot_z = p->get_orientation().z();

		file << cid << ";" << a << ";" << b << ";" << c << ";" << n1 << ";" << n2 << ";" << x << ";" << y << ";" << z << ";" << rot_w << ";" << rot_x << ";" << rot_y << ";" << rot_z << "\n"; 
	}
	file.close();
}
