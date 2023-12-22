#pragma once
#ifndef Scene_H
#define Scene_H

#include <vector>
#include <memory>
#include <string>
#include "TriangleMeshDistance.h"

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Cloth;
class Particle;
class MatrixStack;
class Program;
class Shape;
class Fluid;

class Scene
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	Scene();
	virtual ~Scene();
	
	void load(const std::string &RESOURCE_DIR);
	void init();
	void tare();
	void reset();
	void step();
	
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog) const;
    void setCuboidBoundary(bool flag) {this->cuboidBoundary = flag; }
	double getTime() const { return t; }
	
private:
	double t;
	double h;
	Eigen::Vector3d grav;
    bool cuboidBoundary = true;
    std::shared_ptr<Shape> sphereShape;
    std::shared_ptr<Shape> innerShape;
    std::shared_ptr<Shape> boundaryShape;
    std::shared_ptr<Fluid> fluid;
	std::vector< std::shared_ptr<Particle> > spheres;

};

#endif
