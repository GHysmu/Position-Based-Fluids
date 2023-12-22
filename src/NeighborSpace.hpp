//
//  NeighborSpace.hpp
//  Proj
//
//  Created by 余爽 on 2023/11/19.
//

#ifndef NeighborSpace_hpp
#define NeighborSpace_hpp

#include <stdio.h>
#include <stdio.h>
#include <vector>
#include <memory>
#include <set>


#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <map>


using namespace std;
typedef tuple <int,int,int> key_i;
class Particle;


class NeighborSpace
{
public:
    NeighborSpace(std::vector<std::shared_ptr<Particle> > particles,double cellLen);
    virtual ~NeighborSpace();
    void findNeighbors(shared_ptr<Particle> p,vector<shared_ptr<Particle> > particles, double radius,vector<int> &neighbors, vector<double> &dists);
    void findAllNeighbors(vector<shared_ptr<Particle> > particles, double radius);
    void updateParticle(Eigen::Vector3d preX,std::shared_ptr<Particle> newP);
    void addParticle(key_i key,std::shared_ptr<Particle> P);
    void reset(std::vector<std::shared_ptr<Particle> > particles);
    key_i computeKey(Eigen::Vector3d x);
    std::vector<std::vector<int> > neighborsList;
    
private:
    std::map<key_i , std::set<int>> cellTable;
    
    double cellLen;
    
};

#endif /* NeighborSpace_hpp */
