#pragma once
#ifndef SHAPE_H
#define SHAPE_H

#include <string>
#include <vector>
#include <memory>
#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Program;

/**
 * A shape defined by a list of triangles
 * - posBuf should be of length 3*ntris
 * - norBuf should be of length 3*ntris (if normals are available)
 * - texBuf should be of length 2*ntris (if texture coords are available)
 * posBufID, norBufID, and texBufID are OpenGL buffer identifiers.
 */
class Shape
{
public:
	Shape();
	virtual ~Shape();
	void loadMesh(const std::string &meshName);
	void init();
	void draw(const std::shared_ptr<Program> prog) const;
    void getVerTriInfo(std::vector<std::array<double, 3>> &vertices,
                   std::vector<std::array<int, 3>> &triangles);
    Eigen::Vector3d getNormal(int face_id);
	
private:
	std::vector<float> posBuf;
	std::vector<float> norBuf;
	std::vector<float> texBuf;
    
    std::vector<float> verticesSDF;
    std::vector<int> facesSDF;
    std::vector<float> normalsSDF;
	unsigned posBufID;
	unsigned norBufID;
	unsigned texBufID;
};

#endif
