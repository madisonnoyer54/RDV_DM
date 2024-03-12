#include "../geometry.h"
#include "Model.h"

#ifndef BOUTON_H
#define BOUTON_H

class Bouton : public Model{

    public:
        Vec3f center;
        float radius=0.02;
        float noise_amplitude = 0;

        Bouton(Vec3f _center){
            center = _center;
    
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
            Vec3f light_dir = (lightPos - hit).normalize();
            float light_intensity  = std::max(0.4f, light_dir*distance_field_normal(hit));
            return Vec3f(0,0,0)*light_intensity;
        }

   
        
};

#endif