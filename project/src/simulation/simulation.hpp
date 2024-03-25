#pragma once

#include "cgp/cgp.hpp"



// SPH Particle
struct particle_element
{
    cgp::vec3 p; // Position
    cgp::vec3 v; // Speed
    cgp::vec3 f; // Force

    float rho;      // density at this particle position
    float pressure; // pressure at this particle position
    int bubble_count; // 
    float d;


    particle_element() : p{0,0,0}, v{0,0,0}, f{0,0,0}, rho(0), pressure(0), bubble_count(2), d(0.1) {}
};

// SPH Particle
struct bubble_element
{
    cgp::vec3 p; // Position
    cgp::vec3 v; // Speed
    cgp::vec3 f; // Force

    float r = 0.03; // Display radius.
    float d = 0.04; // Private space radius.

    bubble_element() : p{0,0,0},v{0,0,0},f{0,0,0} {}
};

// SPH simulation parameters
struct sph_parameters_structure
{
    // Influence distance of a particle (size of the kernel)
    float h = 0.12f;

    // Rest density (normalized to 1 - real unit should be 1000kg/m^2)
    float rho0 = 1;

     // Total mass of a particle (consider rho0 h^2)
    float m = rho0*h*h;

    // viscosity parameter
    float nu = 0.02f;   
     
    // Stiffness converting density to pressure
    float stiffness = 8.0f;
    
};


void simulate(float dt, cgp::numarray<particle_element>& particles, cgp::numarray<bubble_element>& bubbles, sph_parameters_structure const& sph_parameters, float pop_coef, bool more_foam);