#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__
#include <cmath>
#include <cassert>

template <size_t DIM, typename T> struct vec {
    vec() { for (size_t i=DIM; i--; data_[i] = T()); }
          T& operator[](const size_t i)       { assert(i<DIM); return data_[i]; }
    const T& operator[](const size_t i) const { assert(i<DIM); return data_[i]; }
private:
    T data_[DIM];
};

template <typename T> struct vec<3,T> {
    vec() : x(T()), y(T()), z(T()) {}
    vec(T X, T Y, T Z) : x(X), y(Y), z(Z) {}
          T& operator[](const size_t i)       { assert(i<3); return i<=0 ? x : (1==i ? y : z); }
    const T& operator[](const size_t i) const { assert(i<3); return i<=0 ? x : (1==i ? y : z); }
    float norm() const { return std::sqrt(x*x+y*y+z*z); }
    vec<3,T> & normalize(T l=1) { *this = (*this)*(l/norm()); return *this; }
    T x,y,z;
};

template<size_t DIM,typename T> T operator*(const vec<DIM,T>& lhs, const vec<DIM,T>& rhs) {
    T ret = T();
    for (size_t i=DIM; i--; ret+=lhs[i]*rhs[i]);
    return ret;
}

template<size_t DIM,typename T>vec<DIM,T> operator+(vec<DIM,T> lhs, const vec<DIM,T>& rhs) {
    for (size_t i=DIM; i--; lhs[i]+=rhs[i]);
    return lhs;
}

template<size_t DIM,typename T>vec<DIM,T> operator-(vec<DIM,T> lhs, const vec<DIM,T>& rhs) {
    for (size_t i=DIM; i--; lhs[i]-=rhs[i]);
    return lhs;
}

template<size_t DIM,typename T,typename U> vec<DIM,T> operator*(const vec<DIM,T> &lhs, const U& rhs) {
    vec<DIM,T> ret;
    for (size_t i=DIM; i--; ret[i]=lhs[i]*rhs);
    return ret;
}

template<size_t DIM,typename T> vec<DIM,T> operator-(const vec<DIM,T> &lhs) {
    return lhs*T(-1);
}

typedef vec<3, float> Vec3f;

template <typename T> inline T lerp(const T &v0, const T &v1, float t) {
    return v0 + (v1-v0)*std::max(0.f, std::min(1.f, t));
}

float hash(const float n) {
    float x = sin(n)*43758.5453f;
    return x-floor(x);
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

float dot(const Vec3f& a, const Vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

float sign(const float x) {
    if(x>0) return 1;
    if(x<0) return -1;
    return 0;
}

#endif //__GEOMETRY_H__

