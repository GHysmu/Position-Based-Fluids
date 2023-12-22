//
//  NeighborSpace.cpp
//  Proj
//
//  Created by 余爽 on 2023/11/19.
//

#include "NeighborSpace.hpp"
#include "Particle.h"

// using namespace Eigen;
using namespace std;


NeighborSpace::NeighborSpace(vector<shared_ptr<Particle> > particles,double cellLen) {
    this->cellLen = cellLen;
    cellTable.clear();
    for(auto p : particles) {
        key_i key = computeKey(p->x);
        addParticle(key, p);
    }
    neighborsList.resize(particles.size());

}
void NeighborSpace::reset(vector<shared_ptr<Particle> > particles) {
    cellTable.clear();
    for(auto p : particles) {
        key_i key = computeKey(p->x);
        addParticle(key, p);
    }

}

NeighborSpace::~NeighborSpace() {
    
}
void NeighborSpace::addParticle(key_i key, shared_ptr<Particle> p) {
    
    if(cellTable.find(key) != cellTable.end()) {
        // contains
         cellTable.find(key)->second.insert(p->i);
    }
    else {
        set<int> curPs;
        curPs.insert(p->i);
        cellTable.insert(std::make_pair(key, curPs));
    }
    
}

void NeighborSpace::updateParticle(Eigen::Vector3d preX,shared_ptr<Particle> newP) {
    key_i preKey = computeKey(preX);
    if(cellTable.find(preKey) != cellTable.end())
        cellTable.find(preKey)->second.erase(newP->i);
    key_i newKey = computeKey(newP->x);
    addParticle(newKey, newP);
    
}

key_i NeighborSpace::computeKey(Eigen::Vector3d x) {
    return key_i((int)x(0)/cellLen,(int)x(1)/cellLen,(int)x(2)/cellLen);
}

void NeighborSpace::findAllNeighbors(vector<shared_ptr<Particle> > particles, double radius) {
#pragma omp parallel for
    for(int i = 0; i < particles.size();i++) {
        vector<int> neighbors;
        vector<double> dists;
        findNeighbors(particles[i], particles, radius, neighbors, dists);
        neighborsList[i] = neighbors;
    }
}

void NeighborSpace::findNeighbors(shared_ptr<Particle> p,vector<shared_ptr<Particle> > particles, double radius,vector<int> &neighbors, vector<double> &dists) {
    neighbors.clear();
    dists.clear();
    for (int i = 0; i < particles.size(); i++) {
        double dist = (p->x - particles[i]->x).norm();
        if (particles[i]->i != p->i && dist < radius) {
            neighbors.push_back(i);
            dists.push_back(dist);
        }
    }
    /*
    key_i key = computeKey(p->x);
    int w = (int)(radius / cellLen) + 1;
    auto [pX, pY, pZ] = key;
    for (int keyX = pX - w; keyX <= pX + w; keyX++) {
            for (int keyY = pY - w; keyY <= pY + w; keyY++) {
                for (int keyZ = pZ - w; keyZ <= pZ+ w; keyZ++) {
                    key_i currKey = key_i(keyX,keyY,keyZ);
                    if (cellTable.size() != 0 &&  cellTable.find(currKey)!= cellTable.end()) {
                        set<int> curNeighbors = cellTable.find(currKey)->second;
                        for (const auto &neighbor : curNeighbors) {
                            double dist = (p->x - particles[neighbor]->x).norm();
                            if (neighbor != p->i && dist < radius) {
                                neighbors.push_back(neighbor);
                                dists.push_back(dist);
                            }
                        }
                    }
                }
            }
        }
    */
    
}
    
