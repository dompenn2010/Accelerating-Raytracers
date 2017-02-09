// [/ignore]

#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>


//#include <cstdlib>
//#include <cstdio>
//#include <cmath>
//#include <iostream>
//#include <cassert>
//#include <chrono>

//#include <fstream>
#include <vector>

struct Vector3f
{
    float x, y, z;
};

void mul_Vec3_Const(Vector3f * v1, Vector3f * v2, float * f)
{
  v2->x = v1->x *f;
  v2->y = v1->y *f;
  v2->z = v1->z *f;

}

void mul_Vec3_Vec3(Vector3f * v1, Vector3f * v2, Vector3f * v3)
{
  v3->x = v1->x * v2->x;
  v3->y = v1->y * v2->y;
  v3->z = v1->z * v2->z;

}

void dot_Vec3(Vector3f * v1, Vector3f * v2, float * f)
{
  f = (v1->x * v2->x) + (v1->y * v2->y) + (v1->z * v2->z);

}


void sub_Vec3_Vec3(Vector3f * v1, Vector3f * v2, Vector3f * v3)
{
  v3->x = v1->x - v2->x;
  v3->y = v1->y - v2->y;
  v3->z = v1->z - v2->z;

}


void add_Vec3_Vec3(Vector3f * v1, Vector3f * v2, Vector3f * v3)
{
  v3->x = v1->x + v2->x;
  v3->y = v1->y + v2->y;
  v3->z = v1->z + v2->z;

}

void addeq_Vec3_Vec3(Vector3f * v1, Vector3f * v2)
{
  v1->x += v2->x;
  v1->y += v2->y;
  v1->z += v2->z;

}

void muleq_Vec3_Vec3(Vector3f * v1, Vector3f * v2)
{
  v1->x *= v2->x;
  v1->y *= v2->y;
  v1->z *= v2->z;

}

void negate_Vec3(Vector3f * v1, Vector3f * v2)
{
  v2->x = v1->x * -1;
  v2->y = v1->y * -1;
  v2->z = v1->z * -1;
}


void length2_Vec3(Vector3f * v1, float * f)
{
  f = ((v1->x * v1->x) + (v1->y * v1->y) + (v1->z * v1->z));
}


void length_Vec3(Vector3f * v1, float * f)
{
  f = sqrt((v1->x * v1->x) + (v1->y * v1->y) + (v1->z * v1->z));
}


int main(int argc, char **argv)
{


  Vector3f vect1 = {.x = 1, .y = 2, .z = 3};
  Vector3f vect2 = {};

  printf("Vect1:\n%d: %.1f\n", 0,vect1.x);
  printf("%d: %.1f\n", 1,vect1.y);
  printf("%d: %.1f\n\n", 2,vect1.z);

  vect2 = mul_Vec3_Const(&vect1, &vect2, 2);

  printf("Vect2:\n%d: %.1f\n", 0,vect2.x);
  printf("%d: %.1f\n", 1,vect2.y);
  printf("%d: %.1f\n", 2,vect2.z);


  return 0;
}
