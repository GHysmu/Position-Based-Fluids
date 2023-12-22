#include <iostream>

#include "GLSL.h"
#include "Scene.h"
#include "Particle.h"
#include "Cloth.h"
#include "Shape.h"
#include "Program.h"
#include "Fluid.hpp"
#include "TriangleMeshDistance.h"
using namespace std;
using namespace Eigen;

Scene::Scene() :
	t(0.0),
	h(1e-2),
	grav(0.0, 0.0, 0.0)
{
}

Scene::~Scene()
{
}

void Scene::load(const string &RESOURCE_DIR)
{
	// Units: meters, kilograms, seconds
    //8 e -3
	h = 10e-3;
	grav << 0.0, -9.8, 0.0;
	double mass = 1;
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<int, 3>> triangles;
    
	sphereShape = make_shared<Shape>();
	sphereShape->loadMesh(RESOURCE_DIR + "sphere2.obj");
    innerShape = make_shared<Shape>();
    boundaryShape = make_shared<Shape>();

    if(cuboidBoundary) {
        innerShape->loadMesh(RESOURCE_DIR + "bunny.obj");
        boundaryShape->loadMesh(RESOURCE_DIR + "cube.obj");
    }
    else {
        innerShape->loadMesh(RESOURCE_DIR + "prince.obj");
        boundaryShape->loadMesh(RESOURCE_DIR + "sphere2.obj");
    }
       
    auto boundary = make_shared<Particle>(boundaryShape);
    spheres.push_back(boundary);
    boundaryShape->getVerTriInfo(vertices, triangles);
    boundary->constructSDF(vertices,triangles);

    auto inner = make_shared<Particle>(innerShape);
    spheres.push_back(inner);
    innerShape->getVerTriInfo(vertices, triangles);
    inner->constructSDF(vertices,triangles);
    if(cuboidBoundary) {
        boundary->r = 2;
        boundary->x = Vector3d(0, 2, 0);
        inner->r = 1.8;
        inner->x = Vector3d(2.5, -0.8, 0);

        fluid = make_shared<Fluid>(12,mass,sphereShape);
        fluid->set(50, 500, 0.005, 0.002, 4, true, 0.00005, 0.01, 0.00002);
    }
    else {
        h = 1e-3;
        boundary->r = 3.5;
        boundary->x = Vector3d(-0.5, 1, 0);
        inner->r = 0.3;
        inner->x = Vector3d(-0.5, 0, 0);
        fluid = make_shared<Fluid>(12,mass,sphereShape);
        fluid->set(10, 400, 0.01, 0.004, 4, true, 0.00005, 0.01, 0.00002);
    }
    
}

void Scene::init()
{
    fluid->init();
    innerShape->init();
    boundaryShape->init();
}

void Scene::tare()
{

    fluid->tare();
}

void Scene::reset()
{
    fluid->reset();
}

void Scene::step()
{
	t += h;
    if(cuboidBoundary) spheres[1]->x(0) -= 0.05 * sin(2 * t);
    fluid->step(h, grav, spheres,cuboidBoundary);
	
}

void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
    glUniform1f(prog->getUniform("alpha"),0.8f);
	glUniform3fv(prog->getUniform("kdFront"), 1, Vector3f(1.0, 1.0, 1.0).data());
    fluid->draw(MV, prog,cuboidBoundary);
    Vector3d boundaryColor = Vector3d();
    Vector3d innerColor = Vector3d();
    float boundaryAlpha;
    
    
    if(cuboidBoundary) {
        boundaryColor = Vector3d();
        innerColor = Vector3d(1.0f, 228.0/255.0f, 185.0f/225.0f);
        boundaryAlpha = 0.1f;
    }
    else {
        boundaryColor = Vector3d(1.0f,1.0f,1.0f);
        innerColor = Vector3d(0.882, 0.368, 0.0f);
        boundaryAlpha = 0.5f;
    }
    
    if(spheres.size() > 1) {
        glUniform1f(prog->getUniform("alpha"),1.0f);
        glUniform3f(prog->getUniform("kdFront"), innerColor(0), innerColor(1), innerColor(2));
        spheres[1]->draw(MV, prog);
    }
    // boundary
    glUniform1f(prog->getUniform("alpha"),boundaryAlpha);
    glUniform3f(prog->getUniform("kdFront"), boundaryColor(0), boundaryColor(1), boundaryColor(2));

    spheres[0]->r *= 1.2;
    spheres[0]->draw(MV, prog);
    spheres[0]->r /= 1.2;



    
    
}
