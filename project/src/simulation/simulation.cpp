#include "simulation.hpp"

using namespace cgp;

// Convert a density value to a pressure
float density_to_pressure(float rho, float rho0, float stiffness)
{
	return stiffness*(rho-rho0);
}

float W_laplacian_viscosity(vec3 const& p_i, vec3 const& p_j, float h)
{
    float const r = norm(p_i-p_j);
    if (r > h) {
        return 0.0;
    }
    return 45.0 / (3.14159f * std::pow(h, 6)) * (h - r);
}

vec3 W_gradient_pressure(vec3 const& p_i, vec3 const& p_j, float h)
{
    float const r = norm(p_i-p_j);
    if (r > h) {
        return vec3(0,0,0);
    }
    return -45.0 / (3.14159f * std::pow(h, 6)) * std::pow(h - r, 2) * (p_i - p_j) / r; 
}

float W_density(vec3 const& p_i, const vec3& p_j, float h)
{
	float const r = norm(p_i-p_j);
    if (r > h) {
        return 0.0f;
    }
    // assert_cgp_no_msg(r<=h);
	return 315.0/(64.0*3.14159f*std::pow(h,9)) * std::pow(h*h-r*r, 3.0f);
}


void update_density(numarray<particle_element>& particles, float h, float m)
{
    // To do: Compute the density value (particles[i].rho) at each particle position
    //  rho_i = \sum_j m W_density(pi,pj)
    int const N = particles.size();
    for(int i=0; i<N; ++i) {
        particles[i].rho = 0.0;
        for (int j = 0; j < N; j++) {
            particles[i].rho += m * W_density(particles[i].p, particles[j].p, h);
        }
    }

}

// Convert the particle density to pressure
void update_pressure(numarray<particle_element>& particles, float rho0, float stiffness)
{
	const int N = particles.size();
    for(int i=0; i<N; ++i)
        particles[i].pressure = density_to_pressure(particles[i].rho, rho0, stiffness);
}

// Compute the forces and update the acceleration of the particles
void update_force(numarray<particle_element>& particles, float h, float m, float nu)
{
	// gravity
    const int N = particles.size();
    for(int i=0; i<N; ++i)
        particles[i].f = m * vec3{0,-9.81f,0};

    //TO Do
    // For all particles i
    //   Compute F_pressure
    //   Compute F_viscosity
    //   particles[i].f += (F_pressure + F_viscosity)
    // ...
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (j == i) continue;
            particles[i].f += -m * m / particles[i].rho / particles[j].rho / 2.0 * (particles[i].pressure + particles[j].pressure) * W_gradient_pressure(particles[i].p, particles[j].p, h);
            particles[i].f += nu * m * m / particles[j].rho * (particles[j].v - particles[i].v) * W_laplacian_viscosity(particles[i].p, particles[j].p, h);
        }
    }

}

void force_bubbles_to_liquid(numarray<particle_element>& particles, numarray<bubble_element>& bubbles) {
    for (int i = 0; i < bubbles.size(); i++) {
        bubbles[i].f = vec3(0, 1.0, 0);
        bool attached_to_particle = false;
        for (int j = 0; j < particles.size(); j++) {
            if (cgp::norm(bubbles[i].p - particles[j].p) < particles[j].d) {
                attached_to_particle = true;
            }
        }
        if (!attached_to_particle) {
            bubbles[i].f = vec3(0, -2, 0);
        }
    }
}

void simulate(float dt, numarray<particle_element>& particles, numarray<bubble_element>& bubbles, sph_parameters_structure const& sph_parameters)
{

	// Update values
    update_density(particles, sph_parameters.h, sph_parameters.m);                   // First compute updated density
    update_pressure(particles, sph_parameters.rho0, sph_parameters.stiffness);       // Compute associated pressure
    update_force(particles, sph_parameters.h, sph_parameters.m, sph_parameters.nu);  // Update forces
    force_bubbles_to_liquid(particles, bubbles);

	// Numerical integration
	float const damping = 0.005f;
    int const N = particles.size();
	float const m = sph_parameters.m;
	for(int k=0; k<N; ++k)
	{
		vec3& p = particles[k].p;
		vec3& v = particles[k].v;
		vec3& f = particles[k].f;

		v = (1-damping)*v + dt*f/m;
		p = p + dt*v;
	}
    for(int k=0; k<bubbles.size(); ++k)
	{
		vec3& p = bubbles[k].p;
		vec3& v = bubbles[k].v;
		vec3& f = bubbles[k].f;

		v = (1-damping)*v + dt*f;
		p = p + dt*v;
	}


	// Collision
    float const epsilon = 1e-3f;
    for(int k=0; k<N; ++k)
    {
        vec3& p = particles[k].p;
        vec3& v = particles[k].v;

        // small perturbation to avoid alignment
        if( p.y<-1 ) {p.y = -1+epsilon*rand_uniform();  v.y *= -0.5f;}
        if( p.x<-1 ) {p.x = -1+epsilon*rand_uniform();  v.x *= -0.5f;}
        if( p.x>1 )  {p.x =  1-epsilon*rand_uniform();  v.x *= -0.5f;}
    }


    // spawning of bubbles
    for (int i = 0; i < N; i++) {
        if (rand_uniform() < 0.005 && particles[i].bubble_count > 0) {
            particles[i].bubble_count--;
            bubble_element new_bubble;
            new_bubble.p = particles[i].p;
            new_bubble.f = vec3(0, 1.0, 0);
            bubbles.push_back(new_bubble);
        }
    }

}