#include "Shape.h"
#include <iostream>
#include <array>
#include "GLSL.h"
#include "Program.h"

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

using namespace std;

Shape::Shape() :
	posBufID(0),
	norBufID(0),
	texBufID(0)
{
}

Shape::~Shape()
{
}

void Shape::loadMesh(const string &meshName)
{
	// Load geometry
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	string errStr;
	bool rc = tinyobj::LoadObj(&attrib, &shapes, &materials, &errStr, meshName.c_str());
	if(!rc) {
		cerr << errStr << endl;
	} else {
		// Some OBJ files have different indices for vertex positions, normals,
		// and texture coordinates. For example, a cube corner vertex may have
		// three different normals. Here, we are going to duplicate all such
		// vertices.
		// Loop over shapes
		for(size_t s = 0; s < shapes.size(); s++) {
			// Loop over faces (polygons)
			size_t index_offset = 0;
			for(size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
				size_t fv = shapes[s].mesh.num_face_vertices[f];
				// Loop over vertices in the face.
				for(size_t v = 0; v < fv; v++) {
					// access to vertex
					tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];

                    
					posBuf.push_back(attrib.vertices[3*idx.vertex_index+0]);
					posBuf.push_back(attrib.vertices[3*idx.vertex_index+1]);
					posBuf.push_back(attrib.vertices[3*idx.vertex_index+2]);
					if(!attrib.normals.empty()) {
						norBuf.push_back(attrib.normals[3*idx.normal_index+0]);
						norBuf.push_back(attrib.normals[3*idx.normal_index+1]);
						norBuf.push_back(attrib.normals[3*idx.normal_index+2]);
					}
					if(!attrib.texcoords.empty()) {
						texBuf.push_back(attrib.texcoords[2*idx.texcoord_index+0]);
						texBuf.push_back(attrib.texcoords[2*idx.texcoord_index+1]);
					}
                    // SDF
                    facesSDF.push_back(idx.vertex_index);
 
				}
				index_offset += fv;
				// per-face material (IGNORE)
				//shapes[s].mesh.material_ids[f];
			}
		}
        
        // for SDF
        verticesSDF = attrib.vertices;
        normalsSDF = attrib.normals;
	}
}

void Shape::init()
{
	// Send the position array to the GPU
	glGenBuffers(1, &posBufID);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size()*sizeof(float), &posBuf[0], GL_STATIC_DRAW);
	
	// Send the normal array to the GPU
	if(!norBuf.empty()) {
		glGenBuffers(1, &norBufID);
		glBindBuffer(GL_ARRAY_BUFFER, norBufID);
		glBufferData(GL_ARRAY_BUFFER, norBuf.size()*sizeof(float), &norBuf[0], GL_STATIC_DRAW);
	}
	
	// Send the texture array to the GPU
	if(!texBuf.empty()) {
		glGenBuffers(1, &texBufID);
		glBindBuffer(GL_ARRAY_BUFFER, texBufID);
		glBufferData(GL_ARRAY_BUFFER, texBuf.size()*sizeof(float), &texBuf[0], GL_STATIC_DRAW);
	}
	
	// Unbind the arrays
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	
	GLSL::checkError(GET_FILE_LINE);
}

void Shape::draw(const shared_ptr<Program> prog) const
{
	GLSL::checkError(GET_FILE_LINE);
	// Bind position buffer
	int h_pos = prog->getAttribute("aPos");
	glEnableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);
	
	// Bind normal buffer
	int h_nor = prog->getAttribute("aNor");
	if(h_nor != -1 && norBufID != 0) {
		glEnableVertexAttribArray(h_nor);
		glBindBuffer(GL_ARRAY_BUFFER, norBufID);
		glVertexAttribPointer(h_nor, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);
	}
	
	// Bind texcoords buffer
	int h_tex = prog->getAttribute("aTex");
	if(h_tex != -1 && texBufID != 0) {
		glEnableVertexAttribArray(h_tex);
		glBindBuffer(GL_ARRAY_BUFFER, texBufID);
		glVertexAttribPointer(h_tex, 2, GL_FLOAT, GL_FALSE, 0, (const void *)0);
	}
	
	// Draw
	int count = (int)posBuf.size()/3; // number of indices to be rendered
	glDrawArrays(GL_TRIANGLES, 0, count);
	
	// Disable and unbind
	if(h_tex != -1) {
		glDisableVertexAttribArray(h_tex);
	}
	if(h_nor != -1) {
		glDisableVertexAttribArray(h_nor);
	}
	glDisableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	
	GLSL::checkError(GET_FILE_LINE);
}

void Shape::getVerTriInfo(std::vector<std::array<double, 3>> &vertices,
                   std::vector<std::array<int, 3>> &triangles) {
    vertices.clear();
    triangles.clear();
    int num = verticesSDF.size();
    cout<<"vertex size"<<num<<endl;
    int index = 0;
    for(int i = 0; i < num; i+=3) {
        std::array<double, 3> ver = {verticesSDF[i],verticesSDF[i+1],verticesSDF[i+2]};
        vertices.push_back(ver);
    }
    for(int i = 0; i < facesSDF.size();i+=3) {
        std::array<int, 3> tri = {facesSDF[i],facesSDF[i + 1],facesSDF[i+2]};
        triangles.push_back(tri);
    }
    cout<<"vertices size"<<vertices.size()<<endl;
    cout<<"triangles size"<<triangles.size()<<endl;

}
Eigen::Vector3d Shape::getNormal(int face_id) {
//    Eigen::Vector3d normal = Eigen::Vector3d(norBuf[face_id *3],norBuf[face_id*3 + 1],norBuf[face_id*3 + 2]);
//    face_id++;
//    normal += Eigen::Vector3d(norBuf[face_id *3],norBuf[face_id*3 + 1],norBuf[face_id*3 + 2]);
//    face_id++;
//    normal += Eigen::Vector3d(norBuf[face_id *3],norBuf[face_id*3 + 1],norBuf[face_id*3 + 2]);
    int v0 = facesSDF[face_id * 3];
    int v1 = facesSDF[face_id * 3 + 1];
    int v2 = facesSDF[face_id * 3 + 2];
    Eigen::Vector3d normal = Eigen::Vector3d(normalsSDF[v0 * 3],normalsSDF[v0 * 3 + 1],normalsSDF[v0 * 3 + 2]);
    
    normal += Eigen::Vector3d(normalsSDF[v1 * 3],normalsSDF[v1 * 3 + 1],normalsSDF[v1 * 3 + 2]);
    normal += Eigen::Vector3d(normalsSDF[v2 * 3],normalsSDF[v2 * 3 + 1],normalsSDF[v2 * 3 + 2]);
    normal.normalize();
    return normal;
}
