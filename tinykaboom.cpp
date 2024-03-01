#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "geometry.h"


const float sphere1_radius   = 0.5; // Rayon pour la petite boule
const float sphere2_radius  = 0.6; // Rayon pour la boule moyenne
const float sphere3_radius   = 0.7; // Rayon pour la grande boule

const float bouton = 0.08; // rayon pour les boutons
//const float sphere_radius   = 0.3;

const float noise_amplitude = 0.15; // texture



template <typename T> inline T lerp(const T &v0, const T &v1, float t) {
    return v0 + (v1-v0)*std::max(0.f, std::min(1.f, t));
}

float hash(const float n) {
    float x = sin(n)*43758.5453f;
    return x-floor(x);
}

float dot(const Vec3f& a, const Vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

float noise(const Vec3f &x) {
    Vec3f p(floor(x.x), floor(x.y), floor(x.z));
    Vec3f f(x.x-p.x, x.y-p.y, x.z-p.z);
    f = f*(f*(Vec3f(3.f, 3.f, 3.f)-f*2.f));
    float n = p*Vec3f(1.f, 57.f, 113.f);
    return lerp(lerp(
                     lerp(hash(n +  0.f), hash(n +  1.f), f.x),
                     lerp(hash(n + 57.f), hash(n + 58.f), f.x), f.y),
                lerp(
                    lerp(hash(n + 113.f), hash(n + 114.f), f.x),
                    lerp(hash(n + 170.f), hash(n + 171.f), f.x), f.y), f.z);
}

Vec3f rotate(const Vec3f &v) {
    return Vec3f(Vec3f(0.00,  0.80,  0.60)*v, Vec3f(-0.80,  0.36, -0.48)*v, Vec3f(-0.60, -0.48,  0.64)*v);
}

float fractal_brownian_motion(const Vec3f &x) { // this is a bad noise function with lots of artifacts. TODO: find a better one
    Vec3f p = rotate(x);
    float f = 0;
    f += 0.5000*noise(p); p = p*2.32;
    f += 0.2500*noise(p); p = p*3.03;
    f += 0.1250*noise(p); p = p*2.61;
    f += 0.0625*noise(p);
    return f/0.9375;
}

Vec3f palette_fire(const float d) { // simple linear gradent yellow-orange-red-darkgray-gray. d is supposed to vary from 0 to 1
    const Vec3f   blanc(1.0, 1.0, 1.0); // note that the color is "hot", i.e. has components >1
   
    float x = std::max(0.f, std::min(1.f, d));
    if (x<.25f)
        return lerp(blanc, blanc, x*4.f);
    else if (x<.5f)
        return lerp(blanc, blanc, x*4.f-1.f);
    else if (x<.75f)
        return lerp(blanc, blanc, x*4.f-2.f);
    return lerp(blanc, blanc, x*4.f-3.f);
}


Vec3f palette_bouton(const float d) { // simple linear gradent yellow-orange-red-darkgray-gray. d is supposed to vary from 0 to 1
    const Vec3f   noir(0.0, 0.0, 0.0); // note that the color is "hot", i.e. has components >1
   
    float x = std::max(0.f, std::min(1.f, d));
    if (x<.25f)
        return lerp(noir, noir, x*4.f);
    else if (x<.5f)
        return lerp(noir, noir, x*4.f-1.f);
    else if (x<.75f)
        return lerp(noir, noir, x*4.f-2.f);
    return lerp(noir, noir, x*4.f-3.f);
}

float signed_distance(const Vec3f &p, const Vec3f &sphere_center, float sphere_radius) { // this function defines the implicit surface we render
    float displacement = -fractal_brownian_motion(p*3.4)*noise_amplitude;
    return (p - sphere_center).norm() - (sphere_radius + displacement);
}

bool sphere_trace(const Vec3f &orig, const Vec3f &dir, Vec3f &pos, const Vec3f &sphere_center, float sphere_radius) {
    if ((orig - sphere_center).norm() - std::pow((orig - sphere_center) * dir, 2) > std::pow(sphere_radius, 2))
        return false;

    pos = orig;

    for (size_t i = 0; i < 64; i++) {
        float d = signed_distance(pos, sphere_center, sphere_radius);
        if (d < 0)
            return true;
        pos = pos + dir * std::max(d * 0.1f, 0.01f);
    }

    return false;
}

float signed_distance_cylindre(const Vec3f &p, const Vec3f &a, const Vec3f &b, float radius) { // this function defines the implicit surface we render
    Vec3f ab = b-a;
    float t = dot(ab, p-a) / dot(ab, ab);
    /*if(t>1) t=1;
    if(t<0) t=0;*/
    Vec3f c = a + ab*t;    
    float x = (p-c).norm() - radius;
    float y = (abs(t-0.5f)-0.5f)*ab.norm();
    float e = Vec3f(std::max(x, 0.f), std::max(y, 0.f), 0.f).norm();
    float i = std::min(std::max(x, y), 0.0f);
    return e + i;
}

bool cylindre_trace(const Vec3f &orig, const Vec3f &dir, Vec3f &pos, const Vec3f &a, const Vec3f &b, float radius) {         // Notice the early discard; in fact I know that the noise() function produces non-negative values,
    pos = orig;
    for (size_t i = 0; i < 64; i++) {
        float d = signed_distance_cylindre(pos, a, b, radius);
        if (d < 0) return true;
        pos = pos + dir * std::max(d * 0.1f, 0.01f);
    }
    return false;
}


Vec3f distance_field_normal(const Vec3f &pos, const Vec3f &sphere_center, float sphere_radius) { // simple finite differences, very sensitive to the choice of the eps constant
    const float eps = 0.1;
    float d = signed_distance(pos, sphere_center, sphere_radius);
    float nx = signed_distance(pos + Vec3f(eps, 0, 0), sphere_center, sphere_radius) - d;
    float ny = signed_distance(pos + Vec3f(0, eps, 0), sphere_center, sphere_radius) - d;
    float nz = signed_distance(pos + Vec3f(0, 0, eps), sphere_center, sphere_radius) - d;
    return Vec3f(nx, ny, nz).normalize();
}

Vec3f distance_field_normal_cylindre(const Vec3f &pos, const Vec3f &a, const Vec3f &b, float radius) { // simple finite differences, very sensitive to the choice of the eps constant
    const float eps = 0.1;
    float d = signed_distance_cylindre(pos, a, b, radius);
    float nx = signed_distance_cylindre(pos + Vec3f(eps, 0, 0), a, b, radius) - d;
    float ny = signed_distance_cylindre(pos + Vec3f(0, eps, 0), a, b, radius) - d;
    float nz = signed_distance_cylindre(pos + Vec3f(0, 0, eps), a, b, radius) - d;
    return Vec3f(nx, ny, nz).normalize();
}

int main() {
	const int   width    = 640;     // image width
	const int   height   = 480;     // image height
	const float fov      = M_PI/3.; // field of view angle
	std::vector<Vec3f> framebuffer(width*height);

	#pragma omp parallel for
	for (size_t j = 0; j<height; j++) { // actual rendering loop
		for (size_t i = 0; i<width; i++) {
		    float dir_x =  (i + 0.5) -  width/2.;
		    float dir_y = -(j + 0.5) + height/2.;    // this flips the image at the same time
		    float dir_z = -height/(2.*tan(fov/2.));
		    Vec3f hit;

            // Bouton 
            if (sphere_trace(Vec3f(0, 0, 3), Vec3f(dir_x, dir_y, dir_z).normalize(), hit,  Vec3f(0, 0, -0.5), 0.08)) { 
                float noise_level = (bouton-hit.norm())/noise_amplitude;
                Vec3f light_dir = (Vec3f(10, 10, 10) - hit).normalize();
                float light_intensity  = std::max(0.4f, light_dir*distance_field_normal(hit,  Vec3f(0, 0, -0.5), 0.08));
                framebuffer[i+j*width] = palette_bouton((-.2 + noise_level)*2)*light_intensity;
            }else
		    
			// Boule 1
			if (sphere_trace(Vec3f(0, 0, 3), Vec3f(dir_x, dir_y, dir_z).normalize(), hit, Vec3f(0, 0.7, 0), sphere1_radius)) { 
			    float noise_level = (sphere1_radius-hit.norm())/noise_amplitude;
			    Vec3f light_dir = (Vec3f(10, 10, 10) - hit).normalize();
			    float light_intensity  = std::max(0.4f, light_dir*distance_field_normal(hit,  Vec3f(0, 0.7, 0), sphere1_radius));
			    framebuffer[i+j*width] = palette_fire((-.2 + noise_level)*2)*light_intensity;
			    
			// Boule 2  
			} else if (sphere_trace(Vec3f(0, 0, 3), Vec3f(dir_x, dir_y, dir_z).normalize(), hit,  Vec3f(0, 0, 0), sphere2_radius)) { 
			    float noise_level = (sphere2_radius-hit.norm())/noise_amplitude;
			    Vec3f light_dir = (Vec3f(10, 10, 10) - hit).normalize();
			    float light_intensity  = std::max(0.4f, light_dir*distance_field_normal(hit, Vec3f(0, 0, 0), sphere2_radius));
			    framebuffer[i+j*width] = palette_fire((-.2 + noise_level)*2)*light_intensity;
			    
			// Boule 3
			} else if (sphere_trace(Vec3f(0, 0, 3), Vec3f(dir_x, dir_y, dir_z).normalize(), hit, Vec3f(0, -0.7, 0), sphere3_radius)) { 
			    float noise_level = (sphere3_radius-hit.norm())/noise_amplitude;
			    Vec3f light_dir = (Vec3f(10, 10, 10) - hit).normalize();
			    float light_intensity  = std::max(0.4f, light_dir*distance_field_normal(hit, Vec3f(0, -0.7, 0), sphere3_radius));
			    framebuffer[i+j*width] = palette_fire((-.2 + noise_level)*2)*light_intensity;
 
            } else if(cylindre_trace(Vec3f(0, 0, 3), Vec3f(dir_x, dir_y, dir_z).normalize(), hit, Vec3f(0, 0, 0), Vec3f(-1.2, 0.45f, 0), 0.03)) {
                Vec3f light_dir = (Vec3f(10, 10, 10) - hit).normalize();
			    float light_intensity  = std::max(0.4f, light_dir*distance_field_normal_cylindre(hit, Vec3f(0, 0, 0), Vec3f(-1.2, 0.45f, 0), 0.03));
                framebuffer[i+j*width] = Vec3f(0.5, 0.4, 0)*light_intensity;

            } else if(cylindre_trace(Vec3f(0, 0, 3), Vec3f(dir_x, dir_y, dir_z).normalize(), hit, Vec3f(-1.2, 0.45f, 0), Vec3f(-1.35, 0.6f, 0), 0.03)) {
                Vec3f light_dir = (Vec3f(10, 10, 10) - hit).normalize();
			    float light_intensity  = std::max(0.4f, light_dir*distance_field_normal_cylindre(hit, Vec3f(-1.2, 0.45f, 0), Vec3f(-1.35, 0.6f, 0), 0.03));
                framebuffer[i+j*width] = Vec3f(0.5, 0.4, 0)*light_intensity;
            }else {
			    framebuffer[i+j*width] = Vec3f(0.2, 0.7, 0.8); // Couleur de fond
			}
		}
	}

	std::ofstream ofs("./out.ppm", std::ios::binary); // save the framebuffer to file
	ofs << "P6\n" << width << " " << height << "\n255\n";
	for (size_t i = 0; i < height*width; ++i) {
		for (size_t j = 0; j<3; j++) {
		    ofs << (char)(std::max(0, std::min(255, static_cast<int>(255*framebuffer[i][j]))));
		}
	}
	ofs.close();

    return 0;
}

