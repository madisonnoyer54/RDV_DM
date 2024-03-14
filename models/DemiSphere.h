#include "../geometry.h"
#include "Model.h"

#ifndef DEMISPHERE_H
#define DEMISPHERE_H

class DemiSphere : public Model {
public:
    Vec3f center;
    float radius;
    float noise_amplitude;
    Vec3f couleur;

    DemiSphere(Vec3f _center, float _radius, float _noise_amplitude, Vec3f _couleur) {
        center = _center;
        radius = _radius;
        couleur = _couleur;
        noise_amplitude = _noise_amplitude;
    }

    float signed_distance(const Vec3f &p) {
        float displacement = -fractal_brownian_motion(p * 3.4) * noise_amplitude;
        float distanceToCenter = (p - center).norm();
        
        // Vérifier si le point est en dessous du centre de la demi-sphère
        if (p.y < center.y) {
            return distanceToCenter - (radius + displacement); // Retourner une distance négative
        } else {
            return distanceToCenter - (radius + displacement); // Retourner une distance positive
        }
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
        float noise_level = (radius - hit.norm()) / noise_amplitude;
        Vec3f light_dir = (lightPos - hit).normalize();
        float light_intensity = std::max(0.4f, light_dir * distance_field_normal(hit));
        return palette_fire((-.2 + noise_level) * 2) * light_intensity;
    }

private:
    Vec3f palette_fire(const float d) {
        const Vec3f blanc(1.0, 1.0, 1.0);

        float x = std::max(0.f, std::min(1.f, d));
        if (x < .25f)
            return lerp(couleur, couleur, x * 4.f);
        else if (x < .5f)
            return lerp(couleur, couleur, x * 4.f - 1.f);
        else if (x < .75f)
            return lerp(couleur, couleur, x * 4.f - 2.f);
        return lerp(couleur, couleur, x * 4.f - 3.f);
    }
};
#endif