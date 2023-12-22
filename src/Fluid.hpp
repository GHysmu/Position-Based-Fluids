//
//  Fluid.hpp
//  Proj
//
//  Created by 余爽 on 2023/11/18.
//

#ifndef Fluid_hpp
#define Fluid_hpp

#include <stdio.h>
#include <vector>
#include <memory>
#include "TriangleMeshDistance.h"

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

class Particle;
class MatrixStack;
class Program;
class Shape;
class NeighborSpace;
class Fluid
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    Fluid(int n,
          double mass,
          std::shared_ptr<Shape>);
    virtual ~Fluid();
    
    void tare();
    void reset();
    void step(double h, const Eigen::Vector3d &grav, const std::vector< std::shared_ptr<Particle> > spheres,bool rec_boundary);
    void init();
    void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p, bool cuboidBoundary) const;
    void set( double rho0,double CFM_EPSILON, double k, double dq, double nPower, bool TENSILE_INSTABILITY, double VORTICITY_EPS, double DELTA, double C);
    
private:
    void computeIncompressibility(std::vector<double> &lamdas,std::vector<Eigen::Vector3d >& deltaPs);
    void computeCollision(const std::vector< std::shared_ptr<Particle> > spheres,bool cuboidBoundary);
    void computeXSPHviscosity(std::vector<Eigen::Vector3d> & deltaVs);
    void computeVorticity(std::vector<Eigen::Vector3d> & deltaFs);
    double tensileInstability(Eigen::Vector3d pij);
    int n;
    const double damping = 1e-5;
    const double pR = 0.1; // radius of particle;
    double rho0 = 50; //  REST_DENSITY = 68378; 6878 // kg m^-2
    double CFM_EPSILON = 100;// Constraint Force Mixing //50000 //600
    double k = 0.005; // Artificial pressure strength 0.0002
    double dq = 0.002; // Artificial pressure radius  0.01
    double nPower = 4; // Artificial pressure power 4
    bool TENSILE_INSTABILITY = true;
    double VORTICITY_EPS = 0.00005;//0.00005
    double DELTA = 0.01; // 0.01
    double C = 0.00002;// The larger C is, the larger the effect of viscosity 0.0002
    const double cellLen = 0.1;
    std::vector< std::shared_ptr<Particle> > particles;
    std::shared_ptr<Shape> particalShape;
    std::shared_ptr<NeighborSpace> neighborSpace;
};

#endif /* Fluid_hpp */
