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
        bubbles[i].f = vec3(0, 35.0, 0);
        int cnt_attached = 0;
        for (int j = 0; j < particles.size(); j++) {
            if (cgp::norm(bubbles[i].p - particles[j].p) < particles[j].d) {
                cnt_attached += 1;
            }
        }
        bool attached_to_particle = cnt_attached > 0;

        float max_speed = 1. + 2. * cnt_attached;
        if (!attached_to_particle) {
            bubbles[i].f = vec3(0, -9.81, 0);

        }

        float speed = std::max(cgp::norm(bubbles[i].v), 0.001f);
        bubbles[i].v = bubbles[i].v / speed * std::min(speed, max_speed);
    }
}

void reflect_bubbles(numarray<bubble_element>& bubbles) {
    for (int i = 0; i < bubbles.size(); i++) {
        int cnt_hit = 0;
        int cnt_refl = 0;
        for (int j = 0; j < bubbles.size(); j++) {
            if (i == j) continue;

            double dist = cgp::norm(bubbles[i].p - bubbles[j].p);
            if (dist < bubbles[i].d + bubbles[j].d) {
                cnt_hit++;
                double push_force = 2 * (dist - bubbles[i].d) / bubbles[j].d - 1;
                push_force = std::max(double(0), std::min(push_force, double(1)));

                auto path_to_i = bubbles[i].p - bubbles[j].p;
                auto normalized_path_to_i = path_to_i / cgp::norm(path_to_i);
                double coef_force_to_i = cgp::dot(normalized_path_to_i, bubbles[j].f);
                auto force_to_i = normalized_path_to_i * coef_force_to_i;
                auto perp_to_i = bubbles[j].f - force_to_i;

                if (coef_force_to_i > 0) {
                    cnt_refl++;
                    bubbles[j].f -= force_to_i;
                    bubbles[j].f += force_to_i * push_force;
                    // if (i < 20) 
                    //     std::cout << "\t\t " << i << ' ' << j << " = " << push_force << std::endl;
                }
            }
        }
        // if (i < 20)
        //     std::cout << i << ": " << cnt_refl << " / " << cnt_hit << std::endl;
    }
}

void reflect_bubbles_2(numarray<bubble_element>& bubbles) {
    for (int i = 0; i < bubbles.size(); i++) {
        int cnt_hit = 0;
        int cnt_refl = 0;
        for (int j = 0; j < bubbles.size(); j++) {
            if (i == j) continue;

            double dist = cgp::norm(bubbles[i].p - bubbles[j].p);
            if (dist < bubbles[i].d + bubbles[j].d) {
                cnt_hit++;
                double push_force = 1 - (dist - bubbles[i].d) / bubbles[j].d;
                push_force = std::max(double(0), std::min(push_force, double(1)));

                auto path_to_i = bubbles[i].p - bubbles[j].p;
                auto normalized_path_to_i = path_to_i / cgp::norm(path_to_i);
                bubbles[j].f -= normalized_path_to_i * push_force;
            }
        }
        // if (i < 20)
        //     std::cout << i << ": " << cnt_hit << std::endl;
    }
}

void pop_bubbles(numarray<particle_element>& particles, numarray<bubble_element>& bubbles, float dt, float pop_coef, bool more_foam) {
    std::set<int, std::greater<int>> bubbles_to_pop;
    for (int i = 0; i < bubbles.size(); i++) {
        int cnt_attached = 0;
        for (int j = 0; j < particles.size(); j++) {
            if (cgp::norm(bubbles[i].p - particles[j].p) < particles[j].d) {
                cnt_attached += 1;
            }
        }
        if (cnt_attached > 0) continue;


        int cnt_hit = 0;
        if (more_foam) {
            for (int j = 0; j < bubbles.size(); j++) {
                if (i == j) continue;

                double dist = cgp::norm(bubbles[i].p - bubbles[j].p);
                if (dist < bubbles[i].d + bubbles[j].d) {
                    cnt_hit++;
                }
            }
        }
        if (cnt_hit <= 2) {
            float prob_pop = std::min(0.5f, pop_coef * (bubbles[i].p.y + 1) / 4);
            prob_pop *= (dt * 6);
            if (rand_uniform() < prob_pop) {
                bubbles_to_pop.insert(i);
                // bubbles[i].p.x = 10 + rand_uniform();
            }
        }
    }

    for (auto i : bubbles_to_pop) {
        std::swap(bubbles[i], bubbles[bubbles.size() - 1]);
        bubbles.resize(bubbles.size() - 1);
    }
}


// TODO: set<int> active_bubbles.
// TODO: high speed particle pops a bubble.

void simulate(float dt, numarray<particle_element>& particles, numarray<bubble_element>& bubbles, sph_parameters_structure const& sph_parameters, float pop_coef, bool more_foam)
{
    // std::cout << "dt = " << dt << std::endl;
    // auto a = cgp::vec3({1, 0, 2});
    // auto b = cgp::rotation_axis_angle(a, 3.1415926);
    // cgp::rotation_transform()
    // std::cout << "a = " << a << std::endl;
    // std::cout << "b = " << b << std::endl;

	// Update values
    update_density(particles, sph_parameters.h, sph_parameters.m);                   // First compute updated density
    update_pressure(particles, sph_parameters.rho0, sph_parameters.stiffness);       // Compute associated pressure
    update_force(particles, sph_parameters.h, sph_parameters.m, sph_parameters.nu);  // Update forces
    force_bubbles_to_liquid(particles, bubbles);
    reflect_bubbles(bubbles);
    reflect_bubbles_2(bubbles);
    pop_bubbles(particles, bubbles, dt, pop_coef, more_foam);

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
    for(int k=0; k<particles.size(); ++k)
    {
        vec3& p = particles[k].p;
        vec3& v = particles[k].v;

        // small perturbation to avoid alignment
        if( p.y<-1 ) {p.y = -1+epsilon*rand_uniform();  v.y *= -0.5f;}
        if( p.x<-1 ) {p.x = -1+epsilon*rand_uniform();  v.x *= -0.5f;}
        if( p.x>1 )  {p.x =  1-epsilon*rand_uniform();  v.x *= -0.5f;}
    }
    for(int k=0; k<bubbles.size(); ++k)
    {
        vec3& p = bubbles[k].p;
        vec3& v = bubbles[k].v;

        if (p.x < 2) {
            // small perturbation to avoid alignment
            if( p.y<-1 ) {p.y = -1+epsilon*rand_uniform();  v.y *= -0.5f;}
            if( p.x<-1 ) {p.x = -1+epsilon*rand_uniform();  v.x *= -0.5f;}
            if( p.x>1 )  {p.x =  1-epsilon*rand_uniform();  v.x *= -0.5f;}
        }
    }


    // spawning of bubbles
    const float p_spawn = dt;
    if (bubbles.size() < 400) {
        for (int i = 0; i < N; i++) {
            if (particles[i].bubble_count > 0 && norm(particles[i].v) < 1 && norm(particles[i].f) < 3 && rand_uniform() < p_spawn) {
                particles[i].bubble_count--;
                bubble_element new_bubble;
                new_bubble.p = particles[i].p;
                new_bubble.f = vec3(0, 1.0, 0);
                // new_bubble.r *= 1 + 0.5 * rand_uniform();
                bubbles.push_back(new_bubble);
            }
        }
    }

}