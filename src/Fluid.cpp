//
//  Fluid.cpp
//  Proj
//
//  Created by 余爽 on 2023/11/18.
//

#include "Fluid.hpp"
#include <iostream>
#include <fstream>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "Shape.h"
#include "Particle.h"
#include "MatrixStack.h"
#include "Program.h"
#include "GLSL.h"
#include "NeighborSpace.hpp"
#include "Kernal.h"

using namespace std;
using namespace Eigen;

Fluid::Fluid(int n,double mass,std::shared_ptr<Shape> particalShape) {

    this->n = n;
    this->particalShape = particalShape;
    int index = 0;
    for(int i = 0; i < n;i++) {
        for(int j =0; j <  n;j++){
            for(int k = 0; k <n;k++) {
                
                auto p = make_shared<Particle>(particalShape);
                p->r = pR;
                p->d = damping;
//                p->x0 = Vector3d(-0.5 + 1.0/(n-1) * i, 0.5 + 1.0/(n-1) * j, -0.5 + 1.0/(n -1)* k);
//                p->v0 = Vector3d(0,0,0);
                //p->x = Vector3d(-0.5 + 1.0 / (n - 1) * i, 0.5 + 1.0 / (n - 1) * j, -0.5 + 1.0 / (n - 1) * k);
                p->x = Vector3d(-1.5 + 2.0 / (n - 1) * i, 1 + 2.0 / (n - 1) * j, -0.5 + 2.0 / (n - 1) * k);
                p->m = mass;
                p->i = index;
                index++;
                p->fixed = false;
                particles.push_back(p);
            }
        }

    }
    // create neighbor hashmap
    neighborSpace = make_shared<NeighborSpace>(this->particles, cellLen);

}

Fluid::~Fluid() {
}

void Fluid::tare() {
    for(auto p : particles) {
        p->tare();
    }
}

void Fluid::reset() {
    for(auto p : particles) {
        p->reset();
    }
    
    neighborSpace->reset(particles);
}

void Fluid::step(double h, const Vector3d &grav, const vector< shared_ptr<Particle> > spheres,bool rec_boundary) {
    // cout<<"Step"<<endl;
    // update based on force
    int num = particles.size();
    vector<Vector3d> preX(num);
    
    // calculate vorticity confinement based on previous step
    vector<Vector3d> deltaFs;
    computeVorticity(deltaFs);

    for(int i = 0; i < num; i++) {
        shared_ptr<Particle> p = particles[i];
        preX[i] = Vector3d(p->x);
        Vector3d force = p->m * grav - p->d * p->v;
        // sif(i % 1 == 0) force += Vector3d(0,-2000,0);
        if(deltaFs.size() == num) {
            force += deltaFs[i];
        }
        if(p->fixed) continue;
        p->v += (h / p->m) * force;
        p->x += h * p->v;
    }
    deltaFs.clear();
    
    // Incompressibility
    neighborSpace->findAllNeighbors(particles,Kernel::RADIUS);

    for(int k = 0; k < 4; k++) {
        std::vector<double> lamdas(num);
        std::vector<Eigen::Vector3d> deltaPs(num);
        computeIncompressibility(lamdas, deltaPs);
        
        for(int i = 0; i< num;i++) {
            shared_ptr<Particle> p = particles[i];
            p->x += deltaPs[i];
        }
        
        // collision detection and response
        computeCollision(spheres,rec_boundary);

    }

    // update v
    for(int i = 0; i < num; i++) {
        shared_ptr<Particle> p = particles[i];
        p->v = (1/h) * (p->x - preX[i]);
        
    }

    
    // XSPH viscosity
    vector<Vector3d> deltaVXSPH;
    computeXSPHviscosity(deltaVXSPH);

    for(int i = 0; i < num; i++) {
        shared_ptr<Particle> p = particles[i];
        p->v += deltaVXSPH[i];
    }
    computeVorticity(deltaFs);
    

//    // update neighbor
//    for(int i = 0; i < num; i++) {
//        shared_ptr<Particle> p = particles[i];
//        if(p->fixed) continue;
//        neighborSpace->updateParticle(preX[i], p);
//    }
    
    
}

void Fluid::init() {
    particalShape->init();
    
}

double Fluid::tensileInstability(Vector3d pij) {
    double wp = Kernel::poly6(pij.norm());
    double wq = Kernel::poly6(this->dq);
    double scorr = -this->k * std::pow(wp / wq, this->nPower);
    return scorr;
}

void Fluid::computeCollision(const vector< shared_ptr<Particle> > spheres, bool rec_boundary) {
    int num = particles.size();
    for(int i = 0; i < num;i++) {
        shared_ptr<Particle> p = particles[i];
        if(p->fixed) continue;
        
        Vector3d normal;
        // inner object
        if(spheres.size() > 1) {
            double innerD = spheres[1]->getSDF(p->x, normal);
            double d2 = (p->x - spheres[1]->x).norm() - spheres[1]->r -p->r;
            if(innerD < 0) p->x -= 1 * innerD * normal ;

        }
        
        // boundary
        double dis = spheres[0]->getSDF(p->x, normal);
        if(dis > 0) p->x += dis * normal;
    }
    
}

void Fluid::computeIncompressibility(vector<double> &lamdas,vector<Vector3d >& deltaPs) {
// #pragma omp parallel
    int num = particles.size();
    
    vector<vector<Vector3d> > dWdps(num);
    // dWdps.clear();
//#pragma omp parallel for

    double gradSqrNormMax = 0;
    for(int i = 0; i< num; i++) {
        deltaPs[i] = Vector3d(0.0, 0.0, 0.0);
        vector<int> neighbors;
        double gradSqrNorm = 0;
        double rhoi = 0;
        Vector3d dWdpi = Vector3d(0,0,0);
        shared_ptr<Particle> p = particles[i];
        neighbors = neighborSpace->neighborsList[i];
        vector<Vector3d> dWdp(neighbors.size());
        // vector<Vector3d> dWdp;
        // for each neighbor
        for (int j = 0; j < neighbors.size();j++) {
            shared_ptr<Particle>  neighbor = particles[neighbors[j]];
            Vector3d pij = p->x - neighbor->x;
            double dis = pij.norm();
            // gradSqrNorm += Kernel::poly6dr(dis)/(rho0 * rho0);
            rhoi += p->m * Kernel::poly6(dis);
            Vector3d dWdpj = Kernel::poly6dp(pij);
            dWdpi += dWdpj;
            dWdp[j] = dWdpj; 
            gradSqrNorm += pow(dWdpj.norm() / rho0,2);
            }
        rhoi += p->m * Kernel::poly6(0);
        double Ci = rhoi / rho0 - 1;
        //if (Ci < 0)
            //Ci = 0;
        dWdpi /= rho0;
        gradSqrNorm += dWdpi.norm() * dWdpi.norm();
        lamdas[i] = -Ci / (gradSqrNorm + CFM_EPSILON);
        dWdps[i] = dWdp;
        if (gradSqrNorm > gradSqrNormMax)
            gradSqrNormMax = gradSqrNorm;
    }

#pragma omp parallel for
    for(int i = 0; i < num;i++) {
        vector<int> neighbors = neighborSpace->neighborsList[i];
        vector<Vector3d> dWdp = dWdps[i];
        // Vector3d deltaP = Vector3d(0,0,0);
        for(int j = 0; j < neighbors.size();j++) {
            Vector3d dpijRho0 = dWdp[j] / rho0;
            double scorr = 0;
            if (TENSILE_INSTABILITY) {
                        scorr = tensileInstability(particles[i]->x - particles[neighbors[j]]->x);
                    }
            deltaPs[i] += (lamdas[i] + lamdas[neighbors[j]] + scorr) * dpijRho0;
            //deltaPs[neighbors[j]] += (lamdas[i] + scorr) * dpijRho0;
        }
    }
     
}

void Fluid::computeXSPHviscosity(vector<Vector3d> & deltaVs) {
    deltaVs.clear();
    //#pragma omp parallel
    for(int i = 0; i< particles.size();i++) {
        shared_ptr<Particle> p = particles[i];
        vector<int> neighbors;
        neighbors = neighborSpace->neighborsList[i];
        Vector3d deltaV = Vector3d(0,0,0);
        for(int j = 0; j < neighbors.size();j++) {
            int neighbor = neighbors[j];
            Vector3d vRel = particles[neighbor]->v - p->v;
            deltaV += vRel * Kernel::poly6((particles[neighbor]->x - p->x).norm());
        }
        deltaV *= C;
        deltaVs.push_back(deltaV);
        
    }
}

void Fluid::computeVorticity(std::vector<Eigen::Vector3d> & deltaFs) {
    deltaFs.clear();
    for(int i = 0; i < particles.size(); i++) {
        shared_ptr<Particle> p = particles[i];
        vector<int> neighbors;
        vector<double> dists;
        neighbors = neighborSpace->neighborsList[i];
        Vector3d omegai = Vector3d(0,0,0);
        Vector3d omegaiPurturbedX = Vector3d(0,0,0);
        Vector3d omegaiPurturbedY = Vector3d(0,0,0);
        Vector3d omegaiPurturbedZ = Vector3d(0,0,0);
        
        for(int j = 0; j<dists.size();j++) {
            int neighbor = neighbors[j];
            Vector3d pij = particles[neighbor]->x - p->x;
            Vector3d dWdpj = Kernel::poly6dp(pij);
            Vector3d vij = particles[neighbor]->v - p->v;
            Vector3d crossij = dWdpj.cross(vij);
            omegai += crossij;
            
            Vector3d pijPurturbedX = Vector3d(pij);
            Vector3d pijPurturbedY = Vector3d(pij);
            Vector3d pijPurturbedZ = Vector3d(pij);
            
            pijPurturbedX(0) += DELTA;
            pijPurturbedY(1) += DELTA;
            pijPurturbedZ(2) += DELTA;
            
            omegaiPurturbedX += Kernel::poly6dp(pijPurturbedX).cross(vij);
            omegaiPurturbedY += Kernel::poly6dp(pijPurturbedY).cross(vij);
            omegaiPurturbedZ += Kernel::poly6dp(pijPurturbedZ).cross(vij);
            
            
        }
        if (omegai(0) == 0 && omegai(1) == 0 && omegai(2) == 0) {
          return;
        }
        Vector3d eta = Vector3d(
             (omegaiPurturbedX.norm() - omegai.norm()) / DELTA,
             (omegaiPurturbedY.norm() - omegai.norm()) / DELTA,
             (omegaiPurturbedZ.norm() - omegai.norm()) / DELTA);
        eta.normalize();
        Vector3d deltaF = Vector3d(0,0,0);
        deltaF = VORTICITY_EPS * eta.cross(omegai);
        deltaFs.push_back(deltaF);
    }
    
}

void Fluid::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, bool cuboidBoundary) const {
    
    Vector3d fluidColor;
    if(cuboidBoundary) fluidColor = Vector3d(30.0/255.0f, 144.0f/255.0f, 1.0f);
    else fluidColor = Vector3d(1.0,1.0f, 1.0f);
    for(auto p : particles) {
        if(!cuboidBoundary  &&  p->i % 10 != 0) continue;
        else {
            glUniform3f(prog->getUniform("kdFront"), fluidColor(0),fluidColor(1), fluidColor(2));
            p->draw(MV, prog);
        }
    }


}

void Fluid::set( double rho0,double CFM_EPSILON, double k, double dq, double nPower, bool TENSILE_INSTABILITY, double VORTICITY_EPS, double DELTA, double C) {
    this->rho0 = rho0;
    this->CFM_EPSILON = CFM_EPSILON;
    this->k = k;
    this->dq = dq;
    this->nPower = nPower;
    this->TENSILE_INSTABILITY = TENSILE_INSTABILITY;
    this->DELTA = DELTA;
    this->C = C;
    
}
