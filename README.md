## Ray tracing project

This is a project given by the CS dept at Université Paul Sabatier to the final year students. The objective of this project is to implement a 3D rendering by ray-tracing. Image generation is done by calculating the intersection of rays with the surfaces of the objects in the 3d scene. Is object has a geometry and a material. At each intersection point, a lighting simulation is done such that the result is a color going in the direction of the ray.

### Content

This project has includes

* Sphere-Ray intersection
* Plane-Ray intersection
* Triangle-Ray intersection
* Shading
* Basic materials
* Anti aliasing


### Usage

All compilation and image generation is in the `IGTAI-RayTracer` directory

To compile all files
```
make
```

To compile and execute unit tests
```
make unit-test
./unit-test
```

To compile the generator
```
make mrt
```

To generate an image
```
./mrt file_name_without_extenstion [scene_number]
```
If the scene number is not entered, scene 0 will be generated:

<img src="https://raw.githubusercontent.com/RainAlexandra/TP_IGTAI/master/myImgs/tp32-scene0.png" width="400" title="Scene 0">
