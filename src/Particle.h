#pragma once
#ifndef Particle_H
#define Particle_H

#include <vector>
#include <memory>
#include "TriangleMeshDistance.h"
#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Shape;
class Program;
class MatrixStack;

class Particle
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	Particle();
	Particle(const std::shared_ptr<Shape> shape);
	virtual ~Particle();
	void tare();
	void reset();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
    void constructSDF(std::vector<std::array<double, 3>> vertices, std::vector<std::array<int, 3>> triangles);
    double getSDF(Eigen::Vector3d pX,Eigen::Vector3d &direction);
	
	double r; // radius
	double m; // mass
	double d; // damping
	Eigen::Vector3d x0; // initial position
	Eigen::Vector3d v0; // initial velocity
	Eigen::Vector3d x;  // position
	Eigen::Vector3d p;  // previous position
	Eigen::Vector3d v;  // velocity
    int i; // index 
	bool fixed;
	
private:
	const std::shared_ptr<Shape> sphere;
    std::shared_ptr<tmd::TriangleMeshDistance> mesh_distance;
};

#endif
