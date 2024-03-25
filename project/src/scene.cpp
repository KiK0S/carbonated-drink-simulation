#include "scene.hpp"
#include <numeric>
#include <algorithm>


using namespace cgp;

const int BIG_BUBBLES = 1;
const bool VIV = false;
const vec3 bubble_color = vec3({0.75, 0.75, 0.75});

void scene_structure::initialize()
{
	camera_projection = camera_projection_orthographic{ -1.1f, 1.1f, -1.1f, 1.1f, -10, 10, window.aspect_ratio() };
	camera_control.initialize(inputs, window); // Give access to the inputs and window global state to the camera controler
	camera_control.look_at({ 0.0f, 0.0f, 2.0f }, {0,0,0}, {0,1,0});
	global_frame.initialize_data_on_gpu(mesh_primitive_frame());

	field.resize(200, 200);
	field_quad.initialize_data_on_gpu(mesh_primitive_quadrangle({ -1,-1,0 }, { 1,-1,0 }, { 1,1,0 }, { -1,1,0 }) );
	field_quad.material.phong = { 1,0,0 };
	field_quad.texture.initialize_texture_2d_on_gpu(field);

	initialize_sph();
	auto sphere = mesh_primitive_sphere(1.0,{0,0,0},10,10);
	sphere_particle.initialize_data_on_gpu(sphere);
	sphere_particle.model.scaling = 0.01f;
	curve_visual.color = { 1,0,0 };
	curve_visual.initialize_data_on_gpu(curve_primitive_circle());

	bubble.initialize_data_on_gpu(curve_primitive_circle(BIG_BUBBLES));  // TODO: remove * 3.
	bubble.color = bubble_color;

	bubble_stick.initialize_data_on_gpu(curve_to_segments({
		{0, 0, 0}, 
		{0.1, 0.1, 0.1}, 
	}));
	bubble_stick.color = bubble_color;


}

void scene_structure::initialize_sph()
{
	// Initial particle spacing (relative to h)
	float const c = 0.5f;
	float const h = sph_parameters.h;


	// Fill a square with particles
	particles.clear();
	bubbles.clear();
	for (float x = -0.5 + h; x < 0.5f - h; x = x + c * h)
	{
		for (float y = -1.0f + h; y < 1.0f - h; y = y + c * h)
		{
			particle_element particle;
			particle.p = { x + h / 8.0 * rand_uniform(),y + h / 8.0 * rand_uniform(),0 }; // a zero value in z position will lead to a 2D simulation
			particles.push_back(particle);
		}
	}

}

void scene_structure::display_voronoi()
{
	// std::cout << "Start Voronoi" << std::endl;

	// if (bubbles.size() >= 3) {
	// 	bubble_element new_bubble;

	// 	new_bubble.p = vec3(0.46, 0, 0);
	// 	bubbles[0] = new_bubble;

	// 	new_bubble.p = vec3(0.54, 0, 0);
	// 	bubbles[1] = new_bubble;

	// 	new_bubble.p = vec3(0.5, 0.04, 0);
	// 	bubbles[2] = new_bubble;

	// 	bubbles.resize(3);	
	// }

	int n = bubbles.size();
	std::vector<std::set<int>> cand(n);  // Canditates to have a Voronoi border with us. 

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) continue;
			if (cgp::norm(bubbles[i].p - bubbles[j].p) < (bubbles[i].r + bubbles[j].r))
				cand[i].insert(j);
		}
	}

	if (VIV)
		std::cout << "Generated neighbours" << std::endl;



	for (int i = 0; i < n; i++) {
		if (VIV)
			std::cout << i << " / " << n << std::endl;
		if (cand[i].empty()) {
			bubble.model.translation = bubbles[i].p;
			bubble.model.scaling = bubbles[i].r;
			bubble.color = bubble_color;
			if (VIV)
				std::cout << "\t\t\t Drawing " << std::endl;
			draw(bubble, environment);
			if (VIV)
				std::cout << "\t\t\t Done " << std::endl;
			continue;
		}

		std::vector<std::pair<vec3, vec3>> hplanes;

		for (auto j : cand[i]) {
			auto c = (bubbles[i].p + bubbles[j].p) / 2;
			auto perp = (bubbles[i].p - bubbles[j].p);
			std::swap(perp.x, perp.y);
			perp.x = -perp.x;

			if (norm(perp) < 0.0001) continue;

			perp = normalize(perp);
			// std::cout << "\t perp = " << perp << std::endl;


			float d1 = norm(c - bubbles[i].p);
			float d2 = sqrt(bubbles[i].r * bubbles[i].r - d1 * d1);

			if (norm(c) < 0.0001 || norm(d2 * perp) < 0.0001) continue;

			hplanes.push_back({c, d2 * perp});

			// bubble_stick.initialize_data_on_gpu(curve_to_segments({
			// 	c - d2 * perp, 
			// 	c + d2 * perp, 
			// }));



			// bubble_stick.color = { 0,1,0 };
			// draw(bubble_stick, environment);
		}
		if (VIV)
			std::cout << "\t Created hplanes" << std::endl;


		std::vector<std::pair<float, float>> segs;
		for (auto [c, hp]: hplanes) {
			float l = -1, r = 1;
			for (auto [ac, ahp]: hplanes) {
				if (sum(c - ac) == 0 && sum(ahp - hp) == 0) continue;

				auto perp = hp;
				std::swap(perp.x, perp.y);
				perp.x = -perp.x;

				auto aperp = ahp;
				std::swap(aperp.x, aperp.y);
				aperp.x = -aperp.x;

				auto path_ac_c = c - ac;


				auto ahp_par_hp = normalize(hp) * dot(normalize(hp), ahp);
				auto ahp_par_perp = normalize(perp) * dot(normalize(perp), ahp);

				auto path_ac_c_par_hp = normalize(hp) * dot(normalize(hp), path_ac_c);
				auto path_ac_c_par_perp = normalize(perp) * dot(normalize(perp), path_ac_c);

				auto path_ac_p = ahp * (norm(path_ac_c_par_perp) / norm(ahp_par_perp));
				auto old_ahp = ahp;
				if (dot(path_ac_c_par_perp, ahp_par_perp) < 0) {
					// std::cout << "SHIIIT IT CAN'T BE THIS" << std::endl;

					ahp = -ahp;
					ahp_par_hp = normalize(hp) * dot(normalize(hp), ahp);
					ahp_par_perp = normalize(perp) * dot(normalize(perp), ahp);
					path_ac_p = ahp * (norm(path_ac_c_par_perp) / norm(ahp_par_perp));
					if (dot(path_ac_c_par_perp, ahp_par_perp) < 0) {
						std::cout << "The logic is broken" << std::endl;
						normalize(vec3({0, 0, 0}));
					}
				}




				auto p = ac + path_ac_p;
				auto path_c_p = p - c;


				// bubble_stick.initialize_data_on_gpu(curve_to_segments({
				// 	ac, 
				// 	c, 
				// }));
				// bubble_stick.color = { 1,0,0 };
				// draw(bubble_stick, environment);

				// bubble_stick.initialize_data_on_gpu(curve_to_segments({
				// 	ac, 
				// 	p, 
				// }));
				// bubble_stick.color = { 0,0,1 };
				// draw(bubble_stick, environment);


				// if (norm(path_c_p) > 0.0001 && norm(perp) > 0.0001) {
				// 	auto npath = normalize(path_c_p);
				// 	auto nperp = normalize(perp);
				// 	std::cout << "\t SURVIVE!!! \t" << nperp << ' ' << npath << ": " << dot(nperp, npath) << std::endl;
				// 	if (std::abs(dot(nperp, npath)) > 0.01) {
				// 		std::cout << nperp << ' ' << npath << ": " << dot(nperp, npath) << std::endl;
				// 		normalize(vec3({0, 0, 0}));
				// 	}
				// }

				// float cf_perp = dot(normalize(perp), ahp);
				float cf_perp = dot(normalize(perp), old_ahp);
				float lim = dot(hp, path_c_p) / dot(hp, hp);
				if (cf_perp < 0)
					r = std::min(r, lim);
				else
					l = std::max(l, lim);
				
			}

			if (l + 0.001 < r && true) {
				bubble_stick.initialize_data_on_gpu(curve_to_segments({
					c + l * hp, 
					c + r * hp, 
				}));
				bubble_stick.color = bubble_color;
				if (VIV)
					std::cout << "\t\t\t Drawing " << std::endl;
				draw(bubble_stick, environment);
				if (VIV)
					std::cout << "\t\t\t Done " << std::endl;
				bubble_stick.clear();

				auto s1 = (c + l * hp) - bubbles[i].p;
				auto s2 = (c + r * hp) - bubbles[i].p;

				float t1 = atan2(s1.y, s1.x);
				float t2 = atan2(s2.y, s2.x);

				if (t1 > t2) std::swap(t1, t2);
				if (t2 - t1 > M_PI) {
					t1 += 2 * M_PI;
					std::swap(t1, t2);
				}
				segs.push_back({t1, t2});
			}
		}

		if (VIV)
			std::cout << "\t Made segs" << std::endl;
		std::sort(segs.begin(), segs.end(), [](std::pair<float, float> a, std::pair<float, float> b) {
			return a.first + a.second < b.first + b.second;
		});

		// std::cout << "segs = ";
		// for (auto p : segs)
		// 	std::cout << "(" << p.first << ' ' << p.second << ") ";
		// std::cout << std::endl;

		for (int k = 0; k < segs.size(); k++) {
			if (VIV)
				std::cout << "\t seg " << k << " / " << segs.size() << std::endl;

			int j = (k + 1) % segs.size();

			float a1 = segs[k].second;
			float a2 = segs[j].first;
			if (std::abs(a2 - a1) < 0.001 || std::abs(std::abs(a2 - a1) - 2 * M_PI) < 0.001) continue;

			while (a2 < a1) a2 += 2 * M_PI;

			a1 *= (180 / M_PI);
			a2 *= (180 / M_PI);
			if (a2 - a1 < 0.1) continue;

			int cnt_seg = std::ceil((a2 - a1) / (360 / 15));

			float step = (a2 - a1) / cnt_seg;

			for (int j = 0; j < cnt_seg; j++) {
				if (VIV)
					std::cout << "\t\t subseg " << j << " / " << cnt_seg << std::endl;
				float aa1 = a1 + j * step;
				float aa2 = aa1 + step;
				aa1 *= (M_PI / 180);
				aa2 *= (M_PI / 180);

				// std::cout << "\t\t seg " << aa1 << ' ' << aa2 << std::endl;
				// vec3 pp1 = vec3({cos(aa1), sin(aa1), 0}) * 0.025 * BIG_BUBBLES;  // TODO: use radius; 
				// vec3 pp2 = vec3({cos(aa2), sin(aa2), 0}) * 0.025 * BIG_BUBBLES;

				vec3 pp1 = vec3({cos(aa1), sin(aa1), 0}) * bubbles[i].r;  // TODO: use radius; 
				vec3 pp2 = vec3({cos(aa2), sin(aa2), 0}) * bubbles[i].r;

				bubble_stick.initialize_data_on_gpu(curve_to_segments({
					bubbles[i].p + pp1,
					bubbles[i].p + pp2,
				}));
				bubble_stick.color = bubble_color;
				if (VIV)
					std::cout << "\t\t\t Drawing " << std::endl;
				draw(bubble_stick, environment);
				if (VIV)
					std::cout << "\t\t\t Done " << std::endl;
				bubble_stick.clear();
			}
		}
		if (VIV)
			std::cout << "\t finished " << i << std::endl;		
	}

	// std::cout << "Finished Voronoi" << std::endl;


	// for (int k = 0; k < bubbles.size(); k++) {
	// 	bubble.model.translation = bubbles[k].p;
	// 	bubble_stick.model.translation = bubbles[k].p;
	// 	// draw(bubble, environment);
	// 	// draw(bubble_stick, environment);
	// }
}

void scene_structure::display_frame()
{
	// Set the light to the current position of the camera
	environment.light = camera_control.camera_model.position();
	
	timer.update(); // update the timer to the current elapsed time
	float const dt = 0.01f * timer.scale;
	sph_parameters.nu = gui.nu;
	simulate(dt, particles, bubbles, sph_parameters, gui.bubble_pop_coef, gui.more_foam);


	if (gui.display_particles) {
		for (int k = 0; k < particles.size(); ++k) {
			vec3 const& p = particles[k].p;
			sphere_particle.model.translation = p;
			draw(sphere_particle, environment);
			// Citron.
		}
	}

	if (gui.display_radius) {
		curve_visual.model.scaling = sph_parameters.h;
		for (int k = 0; k < particles.size(); k += 10) {
			curve_visual.model.translation = particles[k].p;
			draw(curve_visual, environment);
		}
	}

	bubble.model.scaling = sph_parameters.h;  // TODO: what's this?
	display_voronoi();

	if (gui.display_color) {
		update_field_color(field, particles, bubbles);
		field_quad.texture.update(field);
		draw(field_quad, environment);
	}

}

void scene_structure::display_gui()
{
	ImGui::SliderFloat("Timer scale", &timer.scale, 0.01f, 4.0f, "%0.2f");

	bool const restart = ImGui::Button("Restart");
	if (restart)
		initialize_sph();

	ImGui::Checkbox("Color", &gui.display_color);
	ImGui::Checkbox("Particles", &gui.display_particles);
	ImGui::Checkbox("Radius", &gui.display_radius);
	ImGui::SliderFloat("Viscosity", &gui.nu, 0.0f, 0.5f, "%0.2f");
	ImGui::SliderInt("Spawning bubbles", &gui.bubble_count_per_particle, 1, 5);
	ImGui::SliderFloat("Pop coefficient", &gui.bubble_pop_coef, 0.7f, 4.0f, "%0.2f");
	ImGui::Checkbox("More foam", &gui.more_foam);
}

void scene_structure::update_field_color(grid_2D<vec3>& field, numarray<particle_element> const& particles, numarray<bubble_element> const& bubbles)
{
	field.fill(vec3(0.0, 0.0, 0.0));
	grid_2D<float> f_intensity(field.dimension.x, field.dimension.y);
	f_intensity.fill(0.0);
	int const Nf = int(field.dimension.x);
	for (int k = 0; k < particles.size(); k++) {
		vec3 const& pi = particles[k].p;
		int ci = (pi.x / 2.0f + 0.5f)  * (Nf - 1.0f);
		int cj = (pi.y / 2.0f + 0.5f)  * (Nf - 1.0f);

		int bb = 20;
		for (size_t kx = (size_t)std::max(0, -bb + ci); kx <= (size_t) std::min(Nf - 1, ci + bb); kx++) {
			for (size_t ky = (size_t)std::max(0, -bb + cj); ky <= (size_t) std::min(Nf - 1, cj + bb); ky++) {
				vec3 const p0 = { 2.0f * (kx / (Nf - 1.0f) - 0.5f), 2.0f * (ky / (Nf - 1.0f) - 0.5f), 0.0f };
				float const r = norm(pi - p0) / particles[k].d;
				f_intensity(kx, Nf - 1 - ky) += 0.25f * std::exp(-r * r);
			}
		}
	}

	
	for (int k = 0; k < bubbles.size(); k++) {
		vec3 const& pi = bubbles[k].p;
		
		int ci = (pi.x / 2.0f + 0.5f)  * (Nf - 1.0f);
		int cj = (pi.y / 2.0f + 0.5f)  * (Nf - 1.0f);

		int bb = 20;
		for (size_t kx = (size_t)std::max(0, -bb + ci); kx <= (size_t) std::min(Nf - 1, ci + bb); kx++) {
			for (size_t ky = (size_t)std::max(0, -bb + cj); ky <= (size_t) std::min(Nf - 1, cj + bb); ky++) {
				vec3 const p0 = { 2.0f * (kx / (Nf - 1.0f) - 0.5f), 2.0f * (ky / (Nf - 1.0f) - 0.5f), 0.0f };
				float const r = norm(pi - p0) / (bubbles[k].r / 1.5);
				f_intensity(kx, Nf - 1 - ky) -= 2 * std::exp(-r * r);
			}
		}

	}

	
	for (int kx = 0; kx < Nf; ++kx) {
		for (int ky = 0; ky < Nf; ++ky) {
			float f = f_intensity(kx, Nf - 1 - ky);
			if (f < 0)
					field(kx, Nf - 1 - ky) = environment.background_color - clamp(f, -1., 0.) * (vec3(1., 1., 1.) - environment.background_color);
				else
					field(kx, Nf - 1 - ky) = environment.background_color - clamp(f, 0., 1.) * vec3(0.05, 0.31, 0.91);
		}
	}



	for (int k = 0; k < bubbles.size(); k++) {
		vec3 const& pi = bubbles[k].p;
		
		int ci = (pi.x / 2.0f + 0.5f)  * (Nf - 1.0f);
		int cj = (pi.y / 2.0f + 0.5f)  * (Nf - 1.0f);

		int bb = 15;
		for (size_t kx = (size_t)std::max(0, -bb + ci); kx <= (size_t) std::min(Nf - 1, ci + bb); kx++) {
			for (size_t ky = (size_t)std::max(0, -bb + cj); ky <= (size_t) std::min(Nf - 1, cj + bb); ky++) {
				vec3 const p0 = { 2.0f * (kx / (Nf - 1.0f) - 0.5f), 2.0f * (ky / (Nf - 1.0f) - 0.5f), 0.0f };
				float const r = norm(pi - p0) / (bubbles[k].r / 1.5);
				if (norm(p0 - pi) < bubbles[k].r * 0.7)
					field(kx, Nf - 1 - ky) =  vec3(1, 1, 1);
			}
		}

	}
}

void scene_structure::mouse_move_event()
{
	// Do not mode the camera
	/* 
	if (!inputs.keyboard.shift)
		camera_control.action_mouse_move(environment.camera_view);
		*/
}
void scene_structure::mouse_click_event()
{
	camera_control.action_mouse_click(environment.camera_view);
}
void scene_structure::keyboard_event()
{
	camera_control.action_keyboard(environment.camera_view);
}
void scene_structure::idle_frame()
{
	camera_control.idle_frame(environment.camera_view);
}

