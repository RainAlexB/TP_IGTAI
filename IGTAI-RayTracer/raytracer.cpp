#include "image.h"
#include "kdtree.h"
#include "ray.h"
#include "raytracer.h"
#include "scene_types.h"
#include <stdio.h>

#define PI 3.141592653f
#define MAX_DEPTH 10

/* acne_eps is a small constant used to prevent acne when computing
intersection
or bouncing (add this amount to the position before casting a new ray ! */
const float acne_eps = 1e-4;

bool intersectPlane(Ray *ray, Intersection *intersection, Object *obj) {

  point3 rayOrigin = ray->orig;
  vec3 rayDir = ray->dir;
  vec3 planeNorm = obj->geom.plane.normal;
  float planeDist = obj->geom.plane.dist;
  float t = 0.f;

  if ((dot(rayDir, planeNorm)) != 0)
    t = -((dot(rayOrigin, planeNorm) + planeDist) / dot(rayDir, planeNorm));
  else
    return false;

  if (ray->tmin <= t && ray->tmax >= t){
    ray->tmax = t;
    intersection->normal = planeNorm;
    intersection->position = rayAt(*ray, t);
    intersection->mat = &(obj->mat);
    return true;
  }
  
  return false;

}

bool quadSolver(float * t, float A, float B, float C){
  bool res = false;
  float delta = (B * B) - (4 * A * C);
  if (delta >= 0){
    float sol1 = (-B - sqrtf(delta))/(2.0f * A);
    if (sol1 > 0)
      *t = sol1;
    else
      *t = (-B + sqrtf(delta))/(2.0f * A);
    res = true;
  }
  return res;
}

bool intersectSphere(Ray *ray, Intersection *intersection, Object *obj) {

  point3 rayOrigin = ray->orig;
  vec3 rayDir = ray->dir;
  vec3 sphCenter = obj->geom.sphere.center;
  float sphRadius = obj->geom.sphere.radius;
  float t = 0.0f;
  vec3 diffOC = rayOrigin - sphCenter;

  float B = 2.0f * (dot(rayDir, diffOC));
  float C = dot(diffOC, diffOC) - (sphRadius * sphRadius);
  bool res = quadSolver(&t, 1.0f, B, C);

  if (res){
    if (ray->tmin < t && ray->tmax > t){
      point3 pointInter = rayAt(*ray, t);
      ray->tmax = t;
      intersection->normal = normalize(pointInter - sphCenter);
      intersection->position = pointInter;
      intersection->mat = &(obj->mat);
    } 
    else res = false;
  }

  return res;
}

bool intersectTriangle(Ray *ray, Intersection *intersection, Object *obj){
  point3 rayOrigin = ray->orig;
  vec3 rayDir = ray->dir;
  point3 v0 = obj->geom.triangle.a;
  point3 v1 = obj->geom.triangle.b;
  point3 v2 = obj->geom.triangle.c;

  vec3 edge1, edge2, T, P, Q;
  float t, u, v, det, inv_det;

  // je cherche les vecteurs pour les deux bords de triangle qui partagent
  // le point v0
  edge1 = v1 - v0;
  edge2 = v2 - v0;

  vec3 triNorm = normalize(cross(edge1, edge2));
  float d_n = dot(ray->dir, triNorm);

  if (d_n > -acne_eps && d_n < acne_eps) return false;
  if (d_n >= 0) triNorm *= -1.f;

  // calcul du determinant ainsi que le paramatre u
  // si le det est pres de zero le rayon est parallele au triangle
  P = cross(rayDir, edge2);
  det = dot(edge1, P);
  if (det > -acne_eps && det < acne_eps) return false;

  // calcul de la distance entre v0 et l'origine du rayon
  T = rayOrigin - v0;

  inv_det = 1.0 / det;
  u = dot(T, P) * inv_det;
  if (u < 0.0 || u > 1.0) return false;

  // calcul du parametre v et son appartenance a l'intervalle
  Q = cross(T, edge1);
  v = dot(rayDir, Q) * inv_det;
  if (v < 0.0 || u + v > 1.0) return false;

  t = dot(edge2, Q) * inv_det;
  if (ray->tmin < t && ray->tmax > t){
    ray->tmax = t;
    intersection->normal = triNorm;
    intersection->position = rayAt(*ray, t);
    intersection->mat = &(obj->mat);
    return true;
  }

  return false;
}

bool intersectPyramid(Ray *ray, Intersection *intersection, Object *obj){
  int nbIntersections = 0;
  int triNb = 6;
  for (int i = 0; i < triNb; i++)
    if (intersectTriangle(ray, intersection, obj->geom.pyramid.tri[i])) nbIntersections++;

  return nbIntersections > 0;
}

bool intersectMesh(Ray *ray, Intersection *intersection, Object *obj){
  int nbIntersections = 0;
  int triNb = obj->geom.mesh.nbTriangles;
  for (int i = 0; i < triNb; i++)
    if (intersectTriangle(ray, intersection, obj->geom.mesh.tri[i])) nbIntersections++;

  return nbIntersections > 0;
}

bool intersectScene(const Scene *scene, Ray *ray, Intersection *intersection) {
  bool hasIntersection = false;
  size_t objectCount = scene->objects.size();
  Objects objects = scene->objects;
  size_t i;
  int nbIntersections = 0;

  for (i = 0; i < objectCount; i++){
    if (objects[i]->geom.type == SPHERE) 
      hasIntersection = intersectSphere(ray, intersection, objects[i]);
    else if (objects[i]->geom.type == PLANE)
      hasIntersection = intersectPlane(ray, intersection, objects[i]);
    else if (objects[i]->geom.type == TRIANGLE)
      hasIntersection = intersectTriangle(ray, intersection, objects[i]);
    else if (objects[i]->geom.type == PYRAMID)
      hasIntersection = intersectPyramid(ray, intersection, objects[i]);
    else if (objects[i]->geom.type == MESH)
      hasIntersection = intersectMesh(ray, intersection, objects[i]);
    if (hasIntersection) nbIntersections++;
  }

  return nbIntersections > 0;
}

/* ---------------------------------------------------------------------------
 */
/*
 *	The following functions are coded from Cook-Torrance bsdf model
 *description and are suitable only
 *  for rough dielectrics material (RDM. Code has been validated with Mitsuba
 *renderer)
 */

// Shadowing and masking function. Linked with the NDF. Here, Smith function,
// suitable for Beckmann NDF
float RDM_chiplus(float c) { return (c > 0.f) ? 1.f : 0.f; }

/** Normal Distribution Function : Beckmann
 * NdotH : Norm . Half
 */
float RDM_Beckmann(float NdotH, float alpha) {
  //x2 <=> x^2
  float beck = 0.0f;
  float cos2 = NdotH * NdotH;
  float cos4 = cos2 * cos2;
  float tan2 = (1 - cos2)/cos2;
  float alpha2 = alpha * alpha;
  float div = exp(-tan2/alpha2);
  float denom = PI * alpha2 * cos4;
  beck = div / denom;

  return beck;

}

// Fresnel term computation. Implantation of the exact computation. we can use
// the Schlick approximation
// LdotH : Light . Half
float RDM_Fresnel(float LdotH, float extIOR, float intIOR) {

  float cos_i = LdotH;
  float sin2_t = ((extIOR/intIOR) * (extIOR/intIOR)) * (1 - (cos_i * cos_i));
  if (sin2_t > 1)
    return 1;
  float cos_t = sqrtf(1 - sin2_t);

  float rs_div = ((extIOR * cos_i) - (intIOR * cos_t)) * ((extIOR * cos_i) - (intIOR * cos_t));
  float rs_denom = ((extIOR * cos_i) + (intIOR * cos_t)) * ((extIOR * cos_i) + (intIOR * cos_t));
  float rp_div = ((extIOR * cos_t) - (intIOR * cos_i)) * ((extIOR * cos_t) - (intIOR * cos_i));
  float rp_denom = ((extIOR * cos_t) + (intIOR * cos_i)) * ((extIOR * cos_t) + (intIOR * cos_i));
  
  float rs = rs_div / rs_denom;
  float rp = rp_div / rp_denom;

  float F = (rs + rp) / 2;

  return F;

}

// DdotH : Dir . Half
// HdotN : Half . Norm
float RDM_G1(float DdotH, float DdotN, float alpha) {

  float cos = DdotN;
  float cos2 = cos * cos;
  float tan = sqrtf(1 - cos2)/cos;

  float b = 1 / (alpha * tan);
  float k = DdotH / DdotN;

  float G1;
  if (k > 0 && b < 1.6){
    G1 = (3.535 * b + 2.181 * b * b)/(1 + 2.276 * b + 2.577 * b * b);
  } else G1 = 1;

  G1 *= RDM_chiplus(k);
  return G1;

}

// LdotH : Light . Half
// LdotN : Light . Norm
// VdotH : View . Half
// VdotN : View . Norm
float RDM_Smith(float LdotH, float LdotN, float VdotH, float VdotN,
                float alpha) {

  float smith = RDM_G1(LdotH, LdotN, alpha) * RDM_G1(VdotH, VdotN, alpha);
  return smith;

}

// Specular term of the Cook-torrance bsdf
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdotN : View . Norm
color3 RDM_bsdf_s(float LdotH, float NdotH, float VdotH, float LdotN,
                  float VdotN, Material *m) {

  color3 bsdf_s;
  float alpha = m->roughness;
  float D = RDM_Beckmann(NdotH, alpha);
  float F = RDM_Fresnel(LdotH, 1, m->IOR);
  float G = RDM_Smith(LdotH, LdotN, VdotH, VdotN, alpha);

  float tmp = 4 * LdotN * VdotN;
  bsdf_s = (m->specularColor * D * F * G) / tmp;

  return bsdf_s;

}
// diffuse term of the cook torrance bsdf
color3 RDM_bsdf_d(Material *m) {

  color3 bsdf_d = m->diffuseColor / PI;

  return bsdf_d;

}

// The full evaluation of bsdf(wi, wo) * cos (thetai)
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdtoN : View . Norm
// compute bsdf * cos(Oi)
color3 RDM_bsdf(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN,
                Material *m) {

  color3 bsdf_d = RDM_bsdf_d(m);
  color3 bsdf_s = RDM_bsdf_s(LdotH, NdotH, VdotH, LdotN, VdotN, m);
  color3 bsdf = bsdf_d + bsdf_s;

  return bsdf;

}



color3 shade(vec3 n, vec3 v, vec3 l, color3 lc, Material *mat) {
  color3 ret = color3(0.f);
  // float cos = dot(l, n);

  vec3 h = normalize((v + l)/length(v + l));
  float LdotH = dot(l, h);
  float NdotH = dot(n, h);
  float VdotH = dot(v, h);
  float LdotN = dot(l, n);
  float VdotN = dot(v, n);

  if (LdotN > 0)
    // ret = color3((mat->diffuseColor / PI) * cos * lc);
    ret = lc * RDM_bsdf(LdotH, NdotH, VdotH, LdotN, VdotN, mat) * LdotN;
  
  return ret;
}

//! if tree is not null, use intersectKdTree to compute the intersection instead
//! of intersect scene

color3 trace_ray(Scene *scene, Ray *ray, KdTree *tree) {
  color3 color = color3(0,0,0);
  Intersection intersection;

  if (ray->depth > MAX_DEPTH)
    return color;

  if (intersectScene(scene, ray, &intersection)){
    // color = color3(0.5f * intersection.normal + 0.5f); /*tp13*/
    vec3 vue = -(ray->dir);
    vec3 surfNorm = intersection.normal;
    Material *mat = intersection.mat;
    point3 pointInter = intersection.position;
    size_t i;
    size_t lightCount = scene->lights.size();
    Lights lights = scene->lights;
    for (i = 0; i < lightCount; i++){
      color3 lc = lights[i]->color;
      point3 lPos = lights[i]->position;
      vec3 lDir = normalize(lPos - pointInter);
      Ray shadowRay;
      rayInit(&shadowRay, pointInter, lDir, acne_eps, length(lPos - pointInter));
      Intersection shadowInter;
      if (!intersectScene(scene, &shadowRay, &shadowInter))
        color += shade(surfNorm, vue, lDir, lc, mat);
    }

    vec3 reflDir = normalize(reflect(ray->dir, intersection.normal));
    float LdotH = dot(reflDir, normalize(vue + reflDir));
    Ray newRay; 
    rayInit(&newRay, intersection.position, reflDir, acne_eps, 100000, ray->depth + 1);
    color3 cr = trace_ray(scene, &newRay, tree);
    if (cr.x > 1.0f) cr.x = 1.0f;
    if (cr.y > 1.0f) cr.y = 1.0f;
    if (cr.z > 1.0f) cr.z = 1.0f;
    float F = RDM_Fresnel(LdotH, 1, intersection.mat->IOR);
    color += (F * cr * intersection.mat->specularColor);

  } else {
    color = scene->skyColor;
  }

  return color;
}

void renderImage(Image *img, Scene *scene) {

  //! This function is already operational, you might modify it for antialiasing
  //! and kdtree initializaion
  float aspect = 1.f / scene->cam.aspect;

  KdTree *tree = NULL;


//! \todo initialize KdTree

  float delta_y = 1.f / (img->height * 0.5f);   //! one pixel size
  vec3 dy = delta_y * aspect * scene->cam.ydir; //! one pixel step
  vec3 ray_delta_y = (0.5f - img->height * 0.5f) / (img->height * 0.5f) *
                     aspect * scene->cam.ydir;

  float delta_x = 1.f / (img->width * 0.5f);
  vec3 dx = delta_x * scene->cam.xdir;
  vec3 ray_delta_x =
      (0.5f - img->width * 0.5f) / (img->width * 0.5f) * scene->cam.xdir;


  for (size_t j = 0; j < img->height; j++) {
    if (j != 0)
      printf("\033[A\r");
    float progress = (float)j / img->height * 100.f;
    printf("progress\t[");
    int cpt = 0;
    for (cpt = 0; cpt < progress; cpt += 5)
      printf(".");
    for (; cpt < 100; cpt += 5)
      printf(" ");
    printf("]\n");
#pragma omp parallel for
    for (size_t i = 0; i < img->width; i++) {
      color3 *ptr = getPixelPtr(img, i,j);
      
      vec3 ray_dir1 = scene->cam.center + ray_delta_x + ray_delta_y + float(i)*dx + float(j)*dy;
      vec3 ray_dir2 = scene->cam.center + ray_delta_x + ray_delta_y + float(i + 0.25)*dx + float(j + 0.25)*dy;
      vec3 ray_dir3 = scene->cam.center + ray_delta_x + ray_delta_y + float(i - 0.25)*dx + float(j + 0.25)*dy;
      vec3 ray_dir4 = scene->cam.center + ray_delta_x + ray_delta_y + float(i + 0.25)*dx + float(j - 0.25)*dy;
      vec3 ray_dir5 = scene->cam.center + ray_delta_x + ray_delta_y + float(i - 0.25)*dx + float(j - 0.25)*dy;
      
      Ray rx1, rx2, rx3, rx4, rx5;

      rayInit(&rx1, scene->cam.position, normalize(ray_dir1));
      rayInit(&rx2, scene->cam.position, normalize(ray_dir2));
      rayInit(&rx3, scene->cam.position, normalize(ray_dir3));
      rayInit(&rx4, scene->cam.position, normalize(ray_dir4));
      rayInit(&rx5, scene->cam.position, normalize(ray_dir5));

      *ptr += trace_ray(scene, &rx1, tree);
      *ptr += trace_ray(scene, &rx2, tree);
      *ptr += trace_ray(scene, &rx3, tree);
      *ptr += trace_ray(scene, &rx4, tree);
      *ptr += trace_ray(scene, &rx5, tree);
      *ptr /= 5;
    }
  }
}
