#include "Model.h"

#ifndef CYLINDER_H
#define CYLINDER_H

class Cylinder : public Model
{
public:
    Vec3f a, b, color;
    float radius;
    
    Cylinder(Vec3f _a, Vec3f _b, float _radius, Vec3f _color){
        a = _a;
        b = _b;
        radius = _radius;
        color = _color;
    }

    float signed_distance(const Vec3f &p){
        Vec3f ab = b-a;
        float t = dot(ab, p-a) / dot(ab, ab);
        Vec3f c = a + ab*t;    
        float x = (p-c).norm() - radius;
        float y = (abs(t-0.5f)-0.5f)*ab.norm();
        float e = Vec3f(std::max(x, 0.f), std::max(y, 0.f), 0.f).norm();
        float i = std::min(std::max(x, y), 0.0f);
        return e + i;
    }

    Vec3f distance_field_normal(const Vec3f &hit) { // simple finite differences, very sensitive to the choice of the eps constant
        const float eps = 0.001;
        float d = signed_distance(hit);
        float nx = signed_distance(hit + Vec3f(eps, 0, 0)) - d;
        float ny = signed_distance(hit + Vec3f(0, eps, 0)) - d;
        float nz = signed_distance(hit + Vec3f(0, 0, eps)) - d;
        return Vec3f(nx, ny, nz).normalize();
    }

    Vec3f getColor(const Vec3f &hit, const Vec3f &lightPos) {
        Vec3f light_dir = (lightPos - hit).normalize();
        float light_intensity  = std::max(0.4f, light_dir * distance_field_normal(hit));
        return color * light_intensity;
    }
};

#endif