#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "geometry.h"
#include "models/Sphere.h"
#include "models/RoundedCone.h"
#include "models/Cylinder.h"

int envmap_width, envmap_height;
std::vector<Vec3f> envmap;

const float sphere1_radius   = 0.5; // Rayon pour la petite boule
const float sphere2_radius  = 0.6; // Rayon pour la boule moyenne
const float sphere3_radius   = 0.7; // Rayon pour la grande boule



void render(Model** models) {
    
	const int   width    = 960;     // image width
	const int   height   = 720;     // image height
	const float fov      = M_PI/3.; // field of view angle
	std::vector<Vec3f> framebuffer(width*height);

	#pragma omp parallel for
	for (size_t j = 0; j<height; j++) { // actual rendering loop
		for (size_t i = 0; i<width; i++) {
		    float dir_x =  (i + 0.5) -  width/2.;
		    float dir_y = -(j + 0.5) + height/2.;    // this flips the image at the same time
		    float dir_z = -height/(2.*tan(fov/2.));

            // avoir
		    Vec3f hit;

            Vec3f pos = Vec3f(0, 0, 3);
            float minDist = std::numeric_limits<float>::max();

            //Parcours du rayon
            for (size_t k = 0; k < 100; k++) {

                //On teste la distance avec chaque modèle de la scène
                //Et on colorie si il y a une collision
                for(int l = 0; l < 6; l++) {
                    minDist = std::min(models[l]->signed_distance(pos), minDist);
                    if(minDist < 0){
                        framebuffer[i+j*width] = models[l]->getColor(pos, Vec3f(10, 10, 10));
                        break;
                    }
                }
                if(minDist < 0) break;

                //Avancement de la position
                pos = pos + Vec3f(dir_x, dir_y, dir_z).normalize() * std::max(minDist * 0.1f, 0.01f);
            }
            if(minDist >= 0) framebuffer[i+j*width] = Vec3f(0.2, 0.7, 0.8); // Couleur de fond
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
}

int main() {
    
    int n = -1;

    unsigned char *pixmap = stbi_load("../envmap.jpg", &envmap_width, &envmap_height, &n, 0);
    if (!pixmap || 3!=n) {
        std::cerr << "Error: can not load the environment map" << std::endl;
        return -1;
    }
    envmap = std::vector<Vec3f>(envmap_width*envmap_height);
    for (int j = envmap_height-1; j>=0 ; j--) {
        for (int i = 0; i<envmap_width; i++) {
            envmap[i+j*envmap_width] = Vec3f(pixmap[(i+j*envmap_width)*3+0], pixmap[(i+j*envmap_width)*3+1], pixmap[(i+j*envmap_width)*3+2])*(1/255.);
        }
    }
    stbi_image_free(pixmap);

    //Initialisation des objets
    Sphere boule1 = Sphere(Vec3f(0, 0.7, 0), sphere1_radius);
    Sphere boule2 = Sphere(Vec3f(0, 0, 0), sphere2_radius);
    Sphere boule3 = Sphere(Vec3f(0, -0.7, 0), sphere3_radius);

    RoundedCone carrote = RoundedCone(Vec3f(0, 0.7, 0.4), Vec3f(0, 0.7, 0.8), 0.04, 0.1, Vec3f(0.8,0.6,0));

    Cylinder bras1 = Cylinder(Vec3f(0, 0, 0), Vec3f(-1.2, 0.45f, 0), 0.03, Vec3f(0.5, 0.4, 0));
    Cylinder bras2 = Cylinder(Vec3f(-1.2, 0.45f, 0), Vec3f(-1.35, 0.6f, 0), 0.03,Vec3f(0.5, 0.4, 0));

    Model ** models = new Model*[6]{
        &boule1,
        &boule2,
        &boule3,
        &carrote,
        &bras1,
        &bras2
    };

    render(models);

    return 0;
}

