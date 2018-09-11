#include "scene.h"
#include "scene_types.h"
#include <string.h>
#include <algorithm>

#define RATIO 1000

/* ROTATION MATRICES */
// x matrix
mat3 rx(float angle){
	mat3 matRx = mat3(1, 	0, 	0,
					0, cos(angle), -(sin(angle)),
					0, sin(angle), cos(angle));
	return matRx;
}

// y matrix
mat3 ry(float angle){
 	mat3 matRy =  mat3(cos(angle), 	0, 	sin(angle),
						0, 	1, 	0,
					-(sin(angle)), 	0, 	cos(angle));
 	return matRy;
}

// z matrix
mat3 rz(float angle){
	mat3 matRz =  mat3(cos(angle), 	-(sin(angle)), 	0,
						sin(angle), 	cos(angle), 	0,
					 	0, 	0, 	1);
	return matRz;
}

/* OVERALL rotation */
// multiplies the x y z matrices with angle for
mat3 totalRot(vec3 angles){
	return rx(angles.x) * ry(angles.y) * rz(angles.z);
}


Object *initSphere(point3 center, float radius, Material mat) {
	Object *ret;
	ret = (Object *)malloc(sizeof(Object));
	ret->geom.type = SPHERE;
	ret->geom.sphere.center = center;
	ret->geom.sphere.radius = radius;
	memcpy(&(ret->mat), &mat, sizeof(Material));
	return ret;
}

Object *initPlane(vec3 normal, float d, Material mat) {
	Object *ret;
	ret = (Object *)malloc(sizeof(Object));
	ret->geom.type = PLANE;
	ret->geom.plane.normal = normalize(normal);
	ret->geom.plane.dist = d;
	memcpy(&(ret->mat), &mat, sizeof(Material));
	return ret;
}

Object *initTriangle(point3 a, point3 b, point3 c, Material mat, vec3 rotation) {
	Object *ret;
	ret = (Object *)malloc(sizeof(Object));
	ret->geom.type = TRIANGLE;

	a = totalRot(rotation) * a;
	b = totalRot(rotation) * b;
	c = totalRot(rotation) * c;

	ret->geom.triangle.a = a;
	ret->geom.triangle.b = b;
	ret->geom.triangle.c = c;

	memcpy(&(ret->mat), &mat, sizeof(Material));
	return ret;
}

Object *initPyramid(float height, float length, point3 center, Material mat, vec3 rotation){
	Object *ret;
	ret = (Object *)malloc(sizeof(Object));
	ret->geom.type = PYRAMID;
	float tmpLen = length / 2;
	point3 tip = point3(center.x, center.y + height, center.z);
	point3 point1 = point3(center.x + tmpLen, center.y, center.z + tmpLen);
	point3 point2 = point3(center.x - tmpLen, center.y, center.z + tmpLen);
	point3 point3_ = point3(center.x - tmpLen, center.y, center.z - tmpLen);
	point3 point4 = point3(center.x + tmpLen, center.y, center.z - tmpLen);
	ret->geom.pyramid.height = height;
	ret->geom.pyramid.length = length;
	ret->geom.pyramid.center = center;
	// les 4 faces
	ret->geom.pyramid.tri[0] = initTriangle(tip, point1, point2, mat, rotation);
	ret->geom.pyramid.tri[1] = initTriangle(tip, point2, point3_, mat, rotation);
	ret->geom.pyramid.tri[2] = initTriangle(tip, point3_, point4, mat, rotation);
	ret->geom.pyramid.tri[3] = initTriangle(tip, point4, point1, mat, rotation);
	//la base
	ret->geom.pyramid.tri[4] = initTriangle(point1, point3_, point4, mat, rotation);
	ret->geom.pyramid.tri[5] = initTriangle(point1, point2, point3_, mat, rotation);
	return ret;
}

Object *initMesh(const char * fileName, Material mat){
	Object *ret;
	ret = (Object *)malloc(sizeof(Object));
	ret->geom.type = MESH;
	ret->geom.mesh.nbVertices = 0;
	ret->geom.mesh.nbTriangles = 0;
	FILE *file;
	char line[128];
	float x, y, z; //coordinates of a point
	int i, j, k; // indexes of points for a face
	
	file = fopen(fileName, "r");
	
	if (file == NULL){
		fprintf(stderr, "%s not found in working directory\n", fileName);
	} else {
		while (fgets(line, 128, file) != NULL){
			if(line[0] == 'v' && line[1] != 't'){
				sscanf(&line[2], "%f %f %f", &x, &y, &z);
				ret->geom.mesh.pt[ret->geom.mesh.nbVertices] = point3(x/RATIO, y/RATIO, z/RATIO);
				ret->geom.mesh.nbVertices++;
			} else if (line[0] == 'f'){
				sscanf(&line[2], "%d/%*d %d/%*d %d/%*d", &i, &j, &k);                
				ret->geom.mesh.tri[ret->geom.mesh.nbTriangles] = initTriangle(ret->geom.mesh.pt[i-1], 
																			ret->geom.mesh.pt[j-1], 
																			ret->geom.mesh.pt[k-1], mat, vec3(0));
				ret->geom.mesh.nbTriangles++;
			}
		}
	}

	fclose(file);
	return ret;
}

void freeObject(Object *obj) {
	free(obj);
}

Light *initLight(point3 position, color3 color) {
	Light *light = (Light*)malloc(sizeof(Light));
	light->position = position;
	light->color = color;
	return light;
}

void freeLight(Light *light) {
	free(light);
}

Scene * initScene() {
	return new Scene;
}

void freeScene(Scene *scene) {
	std::for_each(scene->objects.begin(), scene->objects.end(), freeObject);
	std::for_each(scene->lights.begin(), scene->lights.end(), freeLight);
	delete scene;
}

void setCamera(Scene *scene, point3 position, point3 at, vec3 up, float fov, float aspect) {
	scene->cam.fov = fov;
	scene->cam.aspect = aspect;
	scene->cam.position = position;
	scene->cam.zdir = normalize(at-position);
	scene->cam.xdir = normalize(cross(up, scene->cam.zdir));
	scene->cam.ydir = normalize(cross(scene->cam.zdir, scene->cam.xdir));
	scene->cam.center = 1.f / tanf ((scene->cam.fov * glm::pi<float>() / 180.f) * 0.5f) * scene->cam.zdir;
}

void addObject(Scene *scene, Object *obj) {
	scene->objects.push_back(obj);
}

void addLight(Scene *scene, Light *light) {
	scene->lights.push_back(light);
}

void setSkyColor(Scene *scene, color3 c) {
	scene->skyColor = c;
}
