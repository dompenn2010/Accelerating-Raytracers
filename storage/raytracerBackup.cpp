// [header]
// A very basic raytracer example.
// [/header]
// [compile]
// c++ -o raytracer -O3 -Wall raytracer.cpp
// [/compile]
// [ignore]
// Copyright (C) 2012  www.scratchapixel.com
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// [/ignore]

#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <cstdlib>

//#include <cstdio>
//#include <cmath>
//#include <iostream>
//#include <cassert>
//#include <chrono>

#include <fstream>
#include <vector>

#if defined __linux__ || defined __APPLE__
// "Compiled for Linux
#else
// Windows doesn't define these values by default, Linux does
#define M_PI 3.141592653589793
#define INFINITY 1e8
#endif

float findMin(float a, float b)
{
    return ((a) < (b) ? (a) : (b));
}

template<typename T>

struct Vec3
{

    T x, y, z;
    Vec3() : x(T(0)), y(T(0)), z(T(0)) {}
    Vec3(T xx) : x(xx), y(xx), z(xx) {}
    Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}
  /*
    Vec3& normalize()
    {
        T nor2 = length2();
        if (nor2 > 0) {
            T invNor = 1 / sqrt(nor2);

            x *= invNor, y *= invNor, z *= invNor;

        }
                    //printf("FIRST: x: %.2f\ny: %.2f\nz: %.2f\n", x, y, z);

        return *this;
    }
*/

//    Vec3<T> operator * (const float &f) const { return Vec3<T>(x * f, y * f, z * f); }
/*
    Vec3<T> operator * (const Vec3<T> &v) const { return Vec3<T>(x * v.x, y * v.y, z * v.z); }
    T dot(const Vec3<T> &v) const { return x * v.x + y * v.y + z * v.z; }
    Vec3<T> operator - (const Vec3<T> &v) const { return Vec3<T>(x - v.x, y - v.y, z - v.z); }
    Vec3<T> operator + (const Vec3<T> &v) const { return Vec3<T>(x + v.x, y + v.y, z + v.z); }
    Vec3<T>& operator += (const Vec3<T> &v) { x += v.x, y += v.y, z += v.z; return *this; }
*/
//    Vec3<T>& operator *= (const Vec3<T> &v) { x *= v.x, y *= v.y, z *= v.z; return *this; }
//    Vec3<T> operator - () const { return Vec3<T>(-x, -y, -z); }
//    T length2() const { return x * x + y * y + z * z; }
//    T length() const { return sqrt(length2()); }

// Do I need this at all?
//    friend std::ostream & operator << (std::ostream &os, const Vec3<T> &v)
//    {
//        os << "[" << v.x << " " << v.y << " " << v.z << "]";
//        return os;
//    }
};

typedef Vec3<float> Vec3f;

/////////////////////////////////////////////////////////////
/*                  TESTING CODE BELOW HERE                */ 
/////////////////////////////////////////////////////////////

struct Vector3f
{
    float x, y, z;
};

void mul_Vec3_Const(const Vec3f * v1, Vec3f * v2, float f)
{
  v2->x = v1->x *f;
  v2->y = v1->y *f;
  v2->z = v1->z *f;

}

void mul_Vec3_Vec3(const Vec3f * v1, const Vec3f * v2, Vec3f * v3)
{
  v3->x = v1->x * v2->x;
  v3->y = v1->y * v2->y;
  v3->z = v1->z * v2->z;

}

void dot_Vec3(const Vec3f * v1, const Vec3f * v2, float * f)
{
  * f = (v1->x * v2->x) + (v1->y * v2->y) + (v1->z * v2->z);

}


void sub_Vec3_Vec3(const Vec3f * v1, const Vec3f * v2, Vec3f * v3)
{
  v3->x = v1->x - v2->x;
  v3->y = v1->y - v2->y;
  v3->z = v1->z - v2->z;

}


void add_Vec3_Vec3(const Vec3f * v1, const Vec3f * v2, Vec3f * v3)
{
  v3->x = v1->x + v2->x;
  v3->y = v1->y + v2->y;
  v3->z = v1->z + v2->z;

}

void addeq_Vec3_Vec3(Vec3f * v1, Vec3f * v2)
{
  v1->x += v2->x;
  v1->y += v2->y;
  v1->z += v2->z;

}

void muleq_Vec3_Vec3(Vec3f * v1, Vec3f * v2)
{
  v1->x *= v2->x;
  v1->y *= v2->y;
  v1->z *= v2->z;

}

void negate_Vec3(Vec3f * v1)
{
  v1->x = v1->x * -1;
  v1->y = v1->y * -1;
  v1->z = v1->z * -1;
}


void length2_Vec3(Vec3f * v1, float * f)
{
  * f = ((v1->x * v1->x) + (v1->y * v1->y) + (v1->z * v1->z));
}


void length_Vec3(Vec3f * v1, float * f)
{
  * f = sqrt((v1->x * v1->x) + (v1->y * v1->y) + (v1->z * v1->z));
}

//     T length2() const { return x * x + y * y + z * z; }

// Need to pass in a vector to be mutated that is = to original
void normalize(Vec3f * v1)
    {
        //T nor2 = length2();
        float nor2;
        length2_Vec3(v1, &nor2);

        if (nor2 > 0) {
 //           T invNor = 1 / sqrt(nor2);
            float invNor = 1 / sqrt(nor2);

            v1->x *= invNor;
            v1->y *= invNor;
            v1->z *= invNor;
        }
             //printf("SECOND x: %.2f\ny: %.2f\nz: %.2f\n",v2->x, v2->y, v2->z);

    }

/////////////////////////////////////////////////////////////
/*                  TESTING CODE ABOVE HERE                */ 
/////////////////////////////////////////////////////////////


struct Sphere
{

    Vec3f center;                           /// position of the sphere
    float radius, radius2;                  /// sphere radius and radius^2
    Vec3f surfaceColor, emissionColor;      /// surface color and emission (light)
    float transparency, reflection;         /// surface transparency and reflectivity
    Sphere(
        const Vec3f &c,
        const float &r,
        const Vec3f &sc,
        const float &refl = 0,
        const float &transp = 0,
        const Vec3f &ec = 0) :
        center(c), radius(r), radius2(r * r), surfaceColor(sc), emissionColor(ec),
        transparency(transp), reflection(refl)
    { /* empty */ }
    //[comment]
    // Compute a ray-sphere intersection using the geometric solution
    //[/comment]
bool intersect(const Vec3f &rayorig, const Vec3f &raydir, float &t0, float &t1) const
    {
        //Vec3f l = center - rayorig;
        
        Vec3f l;
        sub_Vec3_Vec3(&center, &rayorig, &l);

        //float tca = l.dot(raydir);

        float tca;
        dot_Vec3(&l, &raydir, &tca);


        if (tca < 0) return false;
        //float d2 = l.dot(l) - tca * tca;

        float d2;
        dot_Vec3(&l, &l, &d2);
        d2 = d2 - tca * tca;


        if (d2 > radius2) return false;
        float thc = sqrt(radius2 - d2);
        t0 = tca - thc;
        t1 = tca + thc;

        return true;
    }
};

void intersect1(const Vec3f * rayorig, const Vec3f * raydir, float * t0, float * t1, bool * hit, const Sphere * sphere) //const
    {

        //Vec3f l = sphere -> center - rayorig;
        Vec3f l;

        sub_Vec3_Vec3(&sphere ->center, rayorig, &l);


        //float tca = l.dot(raydir);

        float tca;
        dot_Vec3(&l, raydir, &tca);
        
        //printf("FIRST: %.2f\n", tca);
        //printf("SECOND: %.2f\n\n", tca1);



        if (tca < 0)
        {
          *hit = false;
          return;
        }

        //float d2 = l.dot(l) - tca * tca;
        float d2;
        dot_Vec3(&l, &l, &d2);
        d2 = d2 - tca * tca;

        float dd = 1.0;
        //printf("FIRST: %.2f\n", d2);


        //printf("FIRST: %.2f\n", d2);
        //printf("SECOND: %.2f\n\n", d21);

        if (d2 > (sphere ->radius2)) {
          *hit = false;
          return;
        }

        float thc = sqrt((sphere ->radius2) - d2);
        *t0 = tca - thc;
        *t1 = tca + thc;

        *hit = true;
        return;
        
    }

//[comment]
// This variable controls the maximum recursion depth
//[/comment]
#define MAX_RAY_DEPTH 1

float mix(const float &a, const float &b, const float &mix)
{
    return b * mix + a * (1 - mix);
}

//[comment]
// This is the main trace function. It takes a ray as argument (defined by its origin
// and direction). We test if this ray intersects any of the geometry in the scene.
// If the ray intersects an object, we compute the intersection point, the normal
// at the intersection point, and shade this point using this information.
// Shading depends on the surface property (is it transparent, reflective, diffuse).
// The function returns a color for the ray. If the ray intersects an object that
// is the color of the object at the intersection point, otherwise it returns
// the background color.
//[/comment]


//void trace(
Vec3f trace(  
    const Vec3f &rayorig,
    const Vec3f &raydir,
    //const std::vector<Sphere> &spheres,
    const Sphere* spheres,
    int sphereCount,
    const int depth
    //, Vec3f * result
    )
{
    //if (raydir.length() != 1) std::cerr << "Error " << raydir << std::endl;
    float tnear = INFINITY;
    const Sphere* sphere = NULL;
    // find intersection of this ray with the sphere in the scene
    //for (unsigned i = 0; i < spheres.size(); ++i) {
      for (unsigned i = 0; i < sphereCount; ++i) {
        float t0 = INFINITY, t1 = INFINITY;
       
       
        bool sphereIntersect;       
        //printf("1\n");

        intersect1(&rayorig, &raydir, &t0, &t1, &sphereIntersect, &spheres[i]);
        // causes segment fault
        //printf("2\n");
        //printf("1: %d\n", sphereIntersect);
        //fflush(stdout);        //  Flush the stream.
        //printf("2: %d\n\n", spheres[i].intersect(rayorig, raydir, t0, t1));

        //if (spheres[i].intersect(rayorig, raydir, t0, t1)) {
          if (sphereIntersect){
            if (t0 < 0) t0 = t1;
            if (t0 < tnear) {
                tnear = t0;
                sphere = &spheres[i];
            }
        }
    }
    // if there's no intersection return black or background color
    if (!sphere){
        return Vec3f(2);

        //*result = Vec3f(2);

        //return;



    }

    Vec3f surfaceColor = 0; // color of the ray/surfaceof the object intersected by the ray

    //Vec3f phit = rayorig + raydir * tnear; // point of intersection

    Vec3f phit;
    Vec3f raydirInf;
    mul_Vec3_Const(&raydir, &raydirInf, tnear);
    add_Vec3_Vec3(&rayorig, &raydirInf, &phit);

    // printf("x: %.2f\ny: %.2f\nz: %.2f\n",phit.x, phit.y, phit.z);


    //Vec3f nhit = phit - sphere->center; // normal at the intersection point
    //printf("FIRST x: %.2f\ny: %.2f\nz: %.2f\n",nhit.x, nhit.y, nhit.z);
    
    //nhit.x = 0;
    //nhit.y = 0;
    //nhit.z = 0;

    Vec3f nhit;
    sub_Vec3_Vec3(&phit, &sphere->center, &nhit);
    

    //nhit.normalize(); 
   // printf("FIRST\nx: %.2f\ny: %.2f\nz: %.2f\n",nhit.x, nhit.y, nhit.z);

    normalize(&nhit);
    //printf("SECOND\nx: %.2f\ny: %.2f\nz: %.2f\n",nhit2.x, nhit2.y, nhit2.z);

    
    // normalize normal direction
    // If the normal and the view direction are not opposite to each other
    // reverse the normal direction. That also means we are inside the sphere so set
    // the inside bool to true. Finally reverse the sign of IdotN which we want
    // positive.
    float bias = 1e-4; // add some bias to the point from which we will be tracing
    bool inside = false;
   // if (raydir.dot(nhit) > 0){} nhit = -nhit, inside = true;
    float normView;
    dot_Vec3(&raydir, &nhit, &normView);
    
    if (normView > 0)
    {
        //nhit = -nhit;
        negate_Vec3(&nhit);
        inside = true;
    }


    if ((sphere->transparency > 0 || sphere->reflection > 0) && depth < MAX_RAY_DEPTH) {
        //float facingratio = -raydir.dot(nhit);
        float facingratio;
        
        dot_Vec3(&raydir, &nhit, &facingratio);
        facingratio = -facingratio;


        // change the mix value to tweak the effect
        float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1);
        // compute reflection direction (not need to normalize because all vectors
        // are already normalized)

        // v - v * const * const
        //Vec3f refldir = raydir - nhit * 2 * raydir.dot(nhit);
        Vec3f refldir;

        float dotRaydirNhit;
        dot_Vec3(&raydir, &nhit, &dotRaydirNhit);
        
        Vec3f mulNhit;
        mul_Vec3_Const(&nhit, &mulNhit, (dotRaydirNhit * 2));
        sub_Vec3_Vec3(&raydir, &mulNhit, &refldir);

        //printf("FIRST x: %.2f\ny: %.2f\nz: %.2f\n",refldir.x, refldir.y, refldir.z);
        //printf("SECOND x: %.2f\ny: %.2f\nz: %.2f\n",refldir1.x, refldir1.y, refldir1.z);

        //refldir.normalize();

        normalize(&refldir);

            Vec3f rayOrigMul;
{
            Vec3f rayOrigScaled;
            mul_Vec3_Const(&nhit, &rayOrigScaled, bias);  
            add_Vec3_Vec3(&phit, &rayOrigScaled, &rayOrigMul);
            
}
            Vec3f reflection = trace(rayOrigMul, refldir, spheres, sphereCount, depth + 1); 

/*
            Vec3f reflection = trace(phit + nhit * bias, refldir, spheres, depth + 1);
            Vec3f reflection2 = trace(phit + (nhit * bias), refldir, spheres, depth + 1);


            printf("FIRST x: %.2f\ny: %.2f\nz: %.2f\n",reflection1.x, reflection1.y, reflection1.z);
            printf("SECOND x: %.2f\ny: %.2f\nz: %.2f\n",reflection.x, reflection.y, reflection.z);
            printf("THIRD x: %.2f\ny: %.2f\nz: %.2f\n\n",reflection2.x, reflection2.y, reflection2.z);
*/

        
/*
        Vec3f reflection;
        
        {
            Vec3f rayOrigMul;
            add_Vec3_Vec3(&phit, &nhit, &rayOrigMul);
            
            Vec3f rayOrigScaled;
            mul_Vec3_Const(&rayOrigMul, &rayOrigScaled, bias);  

            trace(rayOrigScaled, refldir, spheres, depth + 1, &reflection);  

        }
  */      
        

         
         
        
        Vec3f refraction = 0;
        // if the sphere is also transparent compute refraction ray (transmission)
        if (sphere->transparency) {
            
            float ior = 1.1, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface?
            
            //float cosi = -nhit.dot(raydir);
            float cosi;
        
            dot_Vec3(&nhit, &raydir, &cosi);
            cosi = -cosi;

            float k = 1 - eta * eta * (1 - cosi * cosi);
            
            //                 v    * f   + v   * (f   *   f   -   f)
//            Vec3f refrdir = raydir * eta + nhit * (eta *  cosi - sqrt(k));
            Vec3f refrdir;
            Vec3f rayEta;
            mul_Vec3_Const(&raydir, &rayEta, eta);

            Vec3f nhitEta;
            mul_Vec3_Const(&nhit, &nhitEta, (eta *  cosi - sqrt(k)));
            
            add_Vec3_Vec3(&rayEta, &nhitEta, &refrdir);




            //refrdir.normalize();
            normalize(&refrdir);
            //refraction = trace(phit - nhit * bias, refrdir, spheres, depth + 1);
            
            Vec3f v1;
{
            Vec3f v2;
            mul_Vec3_Const(&nhit, &v2, bias);  
            sub_Vec3_Vec3(&phit, &v2, &v1);

            
}
            refraction = trace(v1, refrdir, spheres, sphereCount, depth + 1); 
            //Vec3f refraction2 = trace(phit - (nhit * bias), refrdir, spheres, depth + 1);

            /*
            Vec3f refraction;
            
            {
                Vec3f rayOrigMul;
                add_Vec3_Vec3(&phit, &nhit, &rayOrigMul);
                
                Vec3f rayOrigScaled;
                mul_Vec3_Const(&rayOrigMul, &rayOrigScaled, bias);  

                trace(rayOrigScaled, refldir, spheres, depth + 1, &refraction);  
            }
            */

        }
        // the result is a mix of reflection and refraction (if the sphere is transparent)
  //      surfaceColor = (
  //          //(Vec3 * float) + (Vec3 * (1 - float) * float) * Vec3
  //          reflection * fresneleffect +
  //          refraction * (1 - fresneleffect) * sphere->transparency) * sphere->surfaceColor;
            
           
        // Lets try this with added scopes
        {
            Vec3f temp1;
            Vec3f temp2;
            Vec3f temp3;
            Vec3f temp4;
            Vec3f temp5;
            float tempF = (1 - fresneleffect) * sphere->transparency;

            mul_Vec3_Const(&reflection, &temp1, fresneleffect);
            mul_Vec3_Const(&refraction, &temp2, tempF);

            add_Vec3_Vec3(&temp1, &temp2, &temp4);
            mul_Vec3_Vec3(&temp4, &sphere->surfaceColor, &surfaceColor);
        }

    }
    else {
        // it's a diffuse object, no need to raytrace any further
        for (unsigned i = 0; i < sphereCount ;++i) {
            if (spheres[i].emissionColor.x > 0) {
                // this is a light
                Vec3f transmission = 1;

                //Vec3f lightDirection = spheres[i].center - phit;
                //printf("FIRST x: %.2f\ny: %.2f\nz: %.2f\n",lightDirection.x, lightDirection.y, lightDirection.z);

                Vec3f lightDirection;

                sub_Vec3_Vec3(&spheres[i].center, &phit, &lightDirection);

                //printf("SECOND x: %.2f\ny: %.2f\nz: %.2f\n",lightDirection1.x, lightDirection1.y, lightDirection1.z);



                //lightDirection.normalize();
                normalize(&lightDirection);

                for (unsigned j = 0; j < sphereCount; ++j) {
                    if (i != j) {
                        float t0, t1;
                        //if (spheres[j].intersect(phit + nhit * bias, lightDirection, t0, t1)) {

                        Vec3f rayOrig;

                    {
                        Vec3f rayOrigScaled;
                        mul_Vec3_Const(&nhit, &rayOrigScaled, bias);  
                        add_Vec3_Vec3(&phit, &rayOrigScaled, &rayOrig);
                                    
                    }
                                          

                        if (spheres[j].intersect(rayOrig, lightDirection, t0, t1)) {
                            transmission = 0;
                            break;
                        }
                    }
                }


// Re-write me
                //Vec3f surfaceColor2;
                //Vec3f surfaceColor = surfaceColor;
/*
                surfaceColor2 = ((sphere->surfaceColor 
                                * transmission) 
                                * std::max(float(0), nhit.dot(lightDirection))) 
                                * spheres[i].emissionColor;
                surfaceColor2 = surfaceColor2 + surfaceColor;

                surfaceColor += sphere->surfaceColor * transmission 
                                *std::max(float(0), nhit.dot(lightDirection)) 
                                * spheres[i].emissionColor;
*/                  
                {
                    Vec3f s1;
                    mul_Vec3_Vec3(&sphere->surfaceColor, &transmission, &s1);

                    Vec3f s2;
                    float nhitDot;
                    dot_Vec3(&nhit, &lightDirection, &nhitDot);


                    //mul_Vec3_Const(&s1, &s2, std::max(float(0), nhit.dot(lightDirection)));
                    mul_Vec3_Const(&s1, &s2, std::max(float(0), nhitDot));

                    Vec3f s3;
                    mul_Vec3_Vec3(&s2, &spheres[i].emissionColor, &s3);
                    //add_Vec3_Vec3(&surfaceColor, &s3, &surfaceColor);
                    //surfaceColor = surfaceColor + s3;
                    add_Vec3_Vec3(&surfaceColor, &s3, &surfaceColor);
                }

                //printf("\nFIRST x: %.2f\ny: %.2f\nz: %.2f\n",surfaceColor.x, surfaceColor.y, surfaceColor.z);
                //printf("SECOND x: %.2f\ny: %.2f\nz: %.2f\n",surfaceColor2.x, surfaceColor2.y, surfaceColor2.z);
                //printf("THIRD x: %.2f\ny: %.2f\nz: %.2f\n\n",surfaceColor1.x, surfaceColor1.y, surfaceColor1.z);



                


            }
        }
    }

    Vec3f result;
    add_Vec3_Vec3(&surfaceColor, &sphere->emissionColor, &result);
    return result;
    //return surfaceColor + sphere->emissionColor;
}


//[comment]
// Main rendering function. We compute a camera ray for each pixel of the image
// trace it and return a color. If the ray hits a sphere, we return the color of the
// sphere at the intersection point, else we return the background color.
//[/comment]
//void render(const std::vector<Sphere> &spheres)
    void render(const Sphere* spheres, int sphereCount)
{
    //unsigned width = 640, height = 480;
    unsigned width = 4096, height = 2160;
    //unsigned width = 2048, height = 1024;
    Vec3f *image = new Vec3f[width * height], *pixel = image;
    float invWidth = 1 / float(width), invHeight = 1 / float(height);
    float fov = 30, aspectratio = width / float(height);
    float angle = tan(M_PI * 0.5 * fov / 180.);
    // Trace rays

//#paragma omp declare target map(tofrom:spheres)
//#pragma omp target teams distribute parallel for

//#pragma omp declare target
//#pragma omp parallel for
 #pragma acc kernels
 for (unsigned y = 0; y < height; y++) {
/*
        if (y % 50 == 0)
        {
            std::cerr << "rendering " << y << " of " << height << std::endl;
        }
*/
       for (unsigned int x = 0; x < width; x++) {
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
            Vec3f raydir(xx, yy, -1);
            
            //raydir.normalize();
            normalize(&raydir);

            
            // *pixel = trace(Vec3f(0), raydir, spheres, 0);

            image[y*width+x] = trace(Vec3f(0), raydir, spheres, sphereCount, 0);

        /*
            {
                trace(Vec3f(0), raydir, spheres, 0, &image[y*width+x]);  
            }
*/

        }
    }
//#pragma omp end declare target
   // Save result to a PPM image (keep these flags if you compile under Windows)
    std::ofstream ofs("./untitled.ppm", std::ios::out | std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (unsigned i = 0; i < width * height; ++i) {
        ofs << (unsigned char)(findMin(float(1), image[i].x) * 255) <<
            (unsigned char)(findMin(float(1), image[i].y) * 255) <<
            (unsigned char)(findMin(float(1), image[i].z) * 255);

/*
// Test to make sure my findMin() works
            std::cout << (std::min(float(1), image[i].x) * 255) - (findMin(float(1), image[i].x) * 255) << std::endl;
            std::cout << (std::min(float(1), image[i].y) * 255) - (findMin(float(1), image[i].y) * 255) << std::endl;
            std::cout << (std::min(float(1), image[i].z) * 255) - (findMin(float(1), image[i].z) * 255) << std::endl;
*/


    }
    ofs.close();
    delete [] image;

}



//[comment]
// In the main function, we will create the scene which is composed of 5 spheres
// and 1 light (which is also a sphere). Then, once the scene description is complete
// we render that scene, by calling the render() function.
//[/comment]
int main(int argc, char **argv)
{
//    srand48(time(NULL) ^ getpid());
//    static const int NUM_SPHERES=8;
/*
    struct Sphere spheressssss[NUM_SPHERES];

    for (int s=0; s<NUM_SPHERES; s++)
    {
        spheressssss[s].center = Vec3f(drand48()*30-15, drand48()*30-15, -50+drand48()*10-5);
        spheressssss[s].radius = 5;
        spheressssss[s].surfaceColor = Vec3f(drand48(), drand48(), drand48());
        spheressssss[s].reflection = .2;
        spheressssss[s].transparency = .5;
    }

*/
/*
/////////// Randomized ///////////////
    
    std::vector<Sphere> spheres;

    for (int s=0; s<NUM_SPHERES; s++)
    {
        spheres.push_back(Sphere(Vec3f(drand48()*30-15, drand48()*30-15, -50+drand48()*10-5), 5, Vec3f(drand48(), drand48(), drand48()), .2, 0.5));
    }
    
/////////// Randomized ///////////////
*/
/////////// NON RANDOM ///////////////
//    srand48(13);
    //std::vector<Sphere> spheres;
    Sphere * spheres;
    int sphereCount = 1000;

    spheres = (Sphere*)malloc(sphereCount * sizeof(Sphere));
    for (int s=0; s<sphereCount; s++)
        {
            spheres[s]=(Sphere(Vec3f(drand48()*30-15, drand48()*30-15, -50+drand48()*10-5), 5, Vec3f(drand48(), drand48(), drand48()), .2, 0.5));
        }

    // position, radius, surface color, reflectivity, transparency, emission color
    /*spheres.push_back(Sphere(Vec3f( 0.0, -10004, -20), 10000, Vec3f(0.20, 0.20, 0.20), 0, 0.0));
    spheres.push_back(Sphere(Vec3f( 0.0,      0, -20),     4, Vec3f(1.00, 0.32, 0.36), 1, 0.5));
    spheres.push_back(Sphere(Vec3f( 5.0,     -1, -15),     2, Vec3f(0.90, 0.76, 0.46), 1, 0.0));
    spheres.push_back(Sphere(Vec3f( 5.0,      0, -25),     3, Vec3f(0.65, 0.77, 0.97), 1, 0.0));
    spheres.push_back(Sphere(Vec3f(-5.5,      0, -15),     3, Vec3f(0.90, 0.90, 0.90), 1, 0.0));
    // light
    spheres.push_back(Sphere(Vec3f( 0.0,     20, -30),     3, Vec3f(0.00, 0.00, 0.00), 0, 0.0, Vec3f(3)));
    */
/*
    spheres[0]=(Sphere(Vec3f( 0.0, -10004, -20), 10000, Vec3f(0.20, 0.20, 0.20), 0, 0.0));
    spheres[1]=(Sphere(Vec3f( 0.0,      0, -20),     4, Vec3f(1.00, 0.32, 0.36), 1, 0.5));
    spheres[2]=(Sphere(Vec3f( 5.0,     -1, -15),     2, Vec3f(0.90, 0.76, 0.46), 1, 0.0));
    spheres[3]=(Sphere(Vec3f( 5.0,      0, -25),     3, Vec3f(0.65, 0.77, 0.97), 1, 0.0));
    spheres[4]=(Sphere(Vec3f(-5.5,      0, -15),     3, Vec3f(0.90, 0.90, 0.90), 1, 0.0));
    // light
    spheres[5]=(Sphere(Vec3f( 0.0,     20, -30),     3, Vec3f(0.00, 0.00, 0.00), 0, 0.0, Vec3f(3)));
    */
/////////// NON RANDOM ///////////////


    //render(spheres);

/*
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    render(spheres);

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1);

    //std::cout <<"The execution time of " << NUM_SPHERES << " spheres was: " << time_span.count() << " seconds." << std::endl;

    std::printf("The execution time of %d spheres was: %.03f seconds.\n" ,NUM_SPHERES, time_span.count());
  */

    clock_t t;
    t = clock();
    render(spheres, sphereCount);
    t = clock() - t;

    double runTime = ((double)t)/ CLOCKS_PER_SEC;

//    std::printf("The execution time of %d spheres was: %.03f seconds.\n" ,NUM_SPHERES, runTime);



    return 0;
}
