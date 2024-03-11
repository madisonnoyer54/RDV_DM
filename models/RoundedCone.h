#include "Model.h"

#ifndef ROUNDEDCONE_H
#define ROUNDEDCONE_H

class RoundedCone : public Model
{
    public: 
        Vec3f a, b;
        float r1, r2;
        Vec3f color;

        RoundedCone(Vec3f _a, Vec3f _b, float _r1, float _r2, Vec3f _color) {
            a = _b;
            b = _a;
            r1 = _r1;
            r2 = _r2;
            color = _color;
        }

        float signed_distance(const Vec3f &p){
            Vec3f ba = b - a;
            float l2 = dot(ba, ba);
            float rr = r1-r2;
            float a2 = l2 - rr*rr;
            float il2 = 1.0/l2;

            Vec3f pa = p - a;
            float y = dot(pa, ba);
            float z = y - l2;
            float x2 = dot(pa*l2 - ba*y, pa*l2 - ba*y);
            float y2 = y*y*l2;
            float z2 = z*z*l2;

            float k = sign(rr)*rr*rr*x2;
            if( sign(z)*a2*z2 > k ) return  sqrt(x2 + z2)        *il2 - r2;
            if( sign(y)*a2*y2 < k ) return  sqrt(x2 + y2)        *il2 - r1;
                                    return (sqrt(x2*a2*il2)+y*rr)*il2 - r1;
        }

        Vec3f distance_field_normal(const Vec3f &hit) {
            const float eps = 0.001;
            float d = signed_distance(hit);
            float nx = signed_distance(hit + Vec3f(eps, 0, 0)) - d;
            float ny = signed_distance(hit + Vec3f(0, eps, 0)) - d;
            float nz = signed_distance(hit + Vec3f(0, 0, eps)) - d;
            return Vec3f(nx, ny, nz).normalize();
        }

        Vec3f getColor(const Vec3f &hit, const Vec3f &lightPos) {
            Vec3f light_dir = (lightPos - hit).normalize();
            float light_intensity  = std::max(0.4f, light_dir * this->distance_field_normal(hit));
            return color * light_intensity;
        }

};

#endif