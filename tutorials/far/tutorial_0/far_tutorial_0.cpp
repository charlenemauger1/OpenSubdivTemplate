//------------------------------------------------------------------------------
// Tutorial description:
//
// This tutorial presents in a very succint way the requisite steps to
// instantiate and refine a mesh with Far from simple topological data.
//

#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/far/primvarRefiner.h>
#include <stdio.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <time.h>
using namespace std;

//------------------------------------------------------------------------------
// Vertex container implementation.
// Struct Vertex contain 3 attributs (3 3D coordinates: x,y,z)
struct Vertex {

	// Minimal required interface ----------------------
	Vertex() { }

	Vertex(Vertex const & src) {
		_position[0] = src._position[0];
		_position[1] = src._position[1];
		_position[2] = src._position[2];
	}

	void Clear(void * = 0) {
		_position[0] = _position[1] = _position[2] = 0.0f;
	}

	void AddWithWeight(Vertex const & src, float weight) {
		_position[0] += weight*src._position[0];
		_position[1] += weight*src._position[1];
		_position[2] += weight*src._position[2];
	}

	// Public interface ------------------------------------
	void SetPosition(float x, float y, float z) {
		_position[0] = x;
		_position[1] = y;
		_position[2] = z;
	}
	
	const float * GetPosition() const {
		return _position;
	}

private:
	float _position[3];
};

//------------------------------------------------------------------------------
// Cube geometry from catmark_cube.h

// coordinates of the nodes
static float g_verts[8][3] = { { -0.5f, -0.5f,  0.5f },
{ 0.5f, -0.5f,  0.5f },
{ -0.5f,  0.5f,  0.5f },
{ 0.5f,  0.5f,  0.5f },
{ -0.5f,  0.5f, -0.5f },
{ 0.5f,  0.5f, -0.5f },
{ -0.5f, -0.5f, -0.5f },
{ 0.5f, -0.5f, -0.5f } };

static int g_nverts = 8,  // number of vertices
g_nfaces = 6;  // number of faces

static int g_vertsperface[6] = { 4, 4, 4, 4, 4, 4 }; // number of vertices per face

// nodes index for each face
static int g_vertIndices[24] = { 0, 1, 3, 2,
2, 3, 5, 4,
4, 5, 7, 6,
6, 7, 1, 0,
1, 7, 5, 3,
6, 0, 2, 4 };


using namespace OpenSubdiv;

// ------------------------------------------------------ MAIN -----------------------------------------------------
//------------------------------------------------------------------------------
int main(int, char **) {

	// Populate a topology descriptor with our raw data
	clock_t tStart = clock();
	typedef Far::TopologyDescriptor Descriptor;

	Sdc::SchemeType type = OpenSubdiv::Sdc::SCHEME_CATMARK;
	Sdc::Options options;
	options.SetVtxBoundaryInterpolation(Sdc::Options::VTX_BOUNDARY_EDGE_ONLY);

	Descriptor desc;
	desc.numVertices = g_nverts;
	desc.numFaces = g_nfaces;
	desc.numVertsPerFace = g_vertsperface;
	desc.vertIndicesPerFace = g_vertIndices;

	// FarTopologyRefiner = refinement of an arbitrary mesh and provides access to the refined mesh topology (building bloc for many useful classes in Far). Can be used 
	// for primvar refinement using PrimvarRefiner (directly) or (indirectly) being used to create a stencil table () or a patch table.

	// Instantiate a FarTopologyRefiner from the descriptor (TopologyRefinerFactory converts mesh topology into the required Far:Topologyrefiner. TopologyRefinerFactory<MESH> )
	Far::TopologyRefiner * refiner = Far::TopologyRefinerFactory<Descriptor>::Create(desc,Far::TopologyRefinerFactory<Descriptor>::Options(type, options));

	int maxlevel = 1;   // number of subdivision

	// Uniformly refine the topolgy up to 'maxlevel'

	refiner->RefineUniform(Far::TopologyRefiner::UniformOptions(maxlevel));


	// Allocate a buffer for vertex primvar data. The buffer length is set to
	// be the sum of all children vertices up to the highest level of refinement.
	std::vector<Vertex> vbuffer(refiner->GetNumVerticesTotal());
	Vertex * verts = &vbuffer[0];


	// Initialize coarse mesh positions
	int nCoarseVerts = g_nverts;
	for (int i = 0; i<nCoarseVerts; ++i) {
		verts[i].SetPosition(g_verts[i][0], g_verts[i][1], g_verts[i][2]);
	}


	// Interpolate vertex primvar data = subdivision (ici utilise une interpolation classique. le poids est renseigne dans la fonction)
	Far::PrimvarRefiner primvarRefiner(*refiner);

	Vertex * src = verts;
	for (int level = 1; level <= maxlevel; ++level) {
		Vertex * dst = src + refiner->GetLevel(level - 1).GetNumVertices();
		primvarRefiner.Interpolate(level, src, dst);  // Interpolate(level,source,destination): apply vertex interpolation weights. Could be Loop, 
		// Catmull-Clark or Bilinear, depends on the sdc (Sdc::SCHEME_LOOP, Sdc::SCHEME_CATMARK, Sdc::SCHEME_BILINEAR)
		src = dst;
	}

    // Write output
	{ // Output OBJ of the highest level refined -----------

		ofstream myfile;
		myfile.open("C:\\Users\\cmau619\\Desktop\\Tuto0.obj");

		if (myfile.is_open())
		{

			Far::TopologyLevel const & refLastLevel = refiner->GetLevel(maxlevel);

			int nverts = refLastLevel.GetNumVertices();
			int nfaces = refLastLevel.GetNumFaces();

			// Print vertex positions
			int firstOfLastVerts = refiner->GetNumVerticesTotal() - nverts;

			for (int vert = 0; vert < nverts; ++vert) 
			{
				float const * pos = verts[firstOfLastVerts + vert].GetPosition();
			//	printf("v %f %f %f\n", pos[0], pos[1], pos[2]);
				myfile << "v " << pos[0]<< " " << pos[1] << " " << pos[2] << "\n";
	//			fprintf(pFile, "v %f %f %f\n", pos[0], pos[1], pos[2]);

			}

			// Print faces
			for (int face = 0; face < nfaces; ++face) 
			{

				Far::ConstIndexArray fverts = refLastLevel.GetFaceVertices(face);

				// all refined Catmark faces should be quads
				assert(fverts.size() == 4);

			//	printf("f ");
				myfile << "f ";
				for (int vert = 0; vert<fverts.size(); ++vert) 
				{
//					printf("%d ", fverts[vert] + 1); // OBJ uses 1-based arrays...
					myfile  << fverts[vert] + 1 << " ";
	//				fprintf(pFile, "%d ", fverts[vert] + 1);
				}
			//	printf("\n");
				myfile << "\n";
			}

			myfile.close();
		}
		else cout << "Unable to open file";
	}
	printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
}

//------------------------------------------------------------------------------