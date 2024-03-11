#include <algorithm>
#include "../geometry.h"

#ifndef MODEL_H
#define MODEL_H

class Model {
    public:
        virtual float signed_distance(const Vec3f &p) = 0;
        virtual Vec3f distance_field_normal(const Vec3f &hit) = 0;
        virtual Vec3f getColor(const Vec3f &hit, const Vec3f &lightPos) = 0;
};

#endif