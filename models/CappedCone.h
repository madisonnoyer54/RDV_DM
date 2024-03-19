#include "Model.h"

#ifndef CAPPEDCONE_H
#define CAPPEDCONE_H

class CappedCone : public Model
{
public:
    Vec3f a, b, color;
    float ra, rb;
    
    CappedCone(Vec3f _a, Vec3f _b, float _ra, float _rb, Vec3f _color){
        a = _a;
        b = _b;
        ra = _ra;
        rb = _rb;
        color = _color;
    }

    float sdCappedCone(const Vec3f &p) {
        float rba  = rb - ra;
        float baba = dot(b - a, b - a);
        float papa = dot(p - a, p - a);
        float paba = dot(p - a, b - a) / baba;

        float x = sqrt(papa - paba * paba * baba);

        float cax = std::max(0.0f, x - ((paba < 0.5f) ? ra : rb));
        float cay = std::abs(paba - 0.5f) - 0.5f;

        float k = rba * rba + baba;
        float f = clamp((rba * (x - ra) + paba * baba) / k, 0.0f, 1.0f);

        float cbx = x - ra - f * rba;
        float cby = paba - f;
        
        float s = (cbx < 0.0f && cay < 0.0f) ? -1.0f : 1.0f;
        
        return s * sqrt(std::min(cax * cax + cay * cay * baba, cbx * cbx + cby * cby * baba));
    }

    float signed_distance(const Vec3f &p) override {
        return sdCappedCone(p);
    }

    Vec3f distance_field_normal(const Vec3f &hit) override {
        const float eps = 0.001f;
        float d = sdCappedCone(hit);
        float nx = sdCappedCone(hit + Vec3f(eps, 0, 0)) - d;
        float ny = sdCappedCone(hit + Vec3f(0, eps, 0)) - d;
        float nz = sdCappedCone(hit + Vec3f(0, 0, eps)) - d;
        return Vec3f(nx, ny, nz).normalize();
    }

    Vec3f getColor(const Vec3f &hit, const Vec3f &lightPos) override {
        Vec3f light_dir = (lightPos - hit).normalize();
        float light_intensity  = std::max(0.4f, light_dir * distance_field_normal(hit));
        return color * light_intensity;
    }

    template<typename T>
    constexpr const T& clamp(const T& v, const T& lo, const T& hi) {
        return (v < lo) ? lo : (hi < v) ? hi : v;
    }

};

#endif
