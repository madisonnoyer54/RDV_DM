#include "../geometry.h"
#include "Model.h"

#ifndef BOUTON_H
#define BOUTON_H

class Bouton : public Model{

    public:
        Vec3f center;
        float radius;
        float noise_amplitude = 0.15;

        Bouton(Vec3f _center, float _radius){
            center = _center;
            radius = _radius;
        }

        float signed_distance(const Vec3f &p){
            float displacement = -fractal_brownian_motion(p*3.4)*noise_amplitude;
            return (p - center).norm() - (radius + displacement);
        }

        Vec3f distance_field_normal(const Vec3f &hit) {
            const float eps = 0.1;
            float d = signed_distance(hit);
            float nx = signed_distance(hit + Vec3f(eps, 0, 0)) - d;
            float ny = signed_distance(hit + Vec3f(0, eps, 0)) - d;
            float nz = signed_distance(hit + Vec3f(0, 0, eps)) - d;
            return Vec3f(nx, ny, nz).normalize();
        }

        Vec3f getColor(const Vec3f &hit, const Vec3f &lightPos) {
            float noise_level = (radius-hit.norm())/noise_amplitude;
            Vec3f light_dir = (lightPos - hit).normalize();
            float light_intensity  = std::max(0.4f, light_dir*distance_field_normal(hit));
            return palette_fire((-.2 + noise_level)*2)*light_intensity;
        }

    private: 
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
        
};

#endif