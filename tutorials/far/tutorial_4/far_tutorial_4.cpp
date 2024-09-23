#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/far/stencilTable.h>
#include <opensubdiv/far/stencilTableFactory.h>
#include <time.h>
#include <opensubdiv/far/primvarRefiner.h>
#include <opensubdiv/far/topologyRefinerFactory.h>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <opensubdiv/far/topologyLevel.h>
#include <opensubdiv/vtr/level.h>
#include <algorithm> // for std::find
#include <iterator> // for std::begin, std::end
#include <opensubdiv/far/patchTableFactory.h>
#include <opensubdiv/far/patchMap.h>
#include "opensubdiv/far/patchTable.h"
#include "opensubdiv/far/patchBasis.h"
#include <opensubdiv/far/ptexIndices.h>
using namespace std;

//------------------------------------------------------------------------------
// Vertex container implementation.
//
struct Vertex {

	// Minimal required interface ----------------------
	Vertex() { }

	Vertex(Vertex const& src) {
		_position[0] = src._position[0];
		_position[1] = src._position[1];
		_position[2] = src._position[2];
	}

	void Clear(void* = 0) {
		_position[0] = _position[1] = _position[2] = 0.0f;
	}

	void AddWithWeight(Vertex const& src, float weight) {
		_position[0] += weight * src._position[0];
		_position[1] += weight * src._position[1];
		_position[2] += weight * src._position[2];
	}

	// Public interface ------------------------------------
	void SetPosition(float x, float y, float z) {
		_position[0] = x;
		_position[1] = y;
		_position[2] = z;
	}

	float const* GetPosition() const {
		return _position;
	}

	float _position[3];
};

struct LimitFrame {

	void Clear(void* = 0) {
		point[0] = point[1] = point[2] = 0.0f;
		deriv1[0] = deriv1[1] = deriv1[2] = 0.0f;
		deriv2[0] = deriv2[1] = deriv2[2] = 0.0f;
	}

	void AddWithWeight(Vertex const& src,
		float weight, float d1Weight, float d2Weight) {
		point[0] += weight * src._position[0];
		point[1] += weight * src._position[1];
		point[2] += weight * src._position[2];

		deriv1[0] += d1Weight * src._position[0];
		deriv1[1] += d1Weight * src._position[1];
		deriv1[2] += d1Weight * src._position[2];

		deriv2[0] += d2Weight * src._position[0];
		deriv2[1] += d2Weight * src._position[1];
		deriv2[2] += d2Weight * src._position[2];

	}


	float point[3],
		deriv1[3],
		deriv2[3];
};

using namespace OpenSubdiv;

Far::TopologyRefiner* createTopologyRefiner(float* g_verts);


//------------------------------------------------------------------------------
int main(int, char**) {

	clock_t tStart = clock();

	float g_verts[100000];


	Far::TopologyRefiner* refiner = createTopologyRefiner(g_verts);
	ofstream myfile;
	myfile.open("C:\\Users\\charl\\Desktop\\bi-a_model.obj");
	/*Far::Index crease[10] = {0,3,12,14,17,18,21,30,32,35};

	for (int i = 0; i < 10; i++)
	{
	refiner->_levels[0]->_vertTags[(Far::Index)(crease[i])]._rule = Sdc::Crease::RULE_CORNER;
	}*/

	// Toy problem

	//Far::Index crease[16] = { 5,25,7,27,6,26,3,23,24,28,29,30,8,4,9,10 };
	//for (int i = 0; i < 16; i++)
	//{
	//	refiner->_levels[0]->_vertTags[(Far::Index)(crease[i])]._rule = Sdc::Crease::RULE_CREASE;
	//}

	// biventricular model
	//Vtr::internal::Level::VTag vtag = refiner->_levels[0]->getVertexTag((Far::Index)(0));    // to get the Tag for a smooth vertex 90
	//Vtr::internal::Level::ETag etag = refiner->_levels[0]->getEdgeTag(refiner->_levels[0]->findEdge((Far::Index)(90), (Far::Index)(69))); // assign smooth edge 69,90


	/*Far::Index tab_ventricle1[21] = { 98,93,102,95,100,118,117,106,92,91,83,128,24,88,87,85,86,89,116,97,98 };
	Far::Index tab_ventricle2[21] = { 82,127,84,81,65,23,74,22,112,17,27,119,122,120,121,104,105,49,50,44,82 };

	//Far::Index tab_epi[11] = { 329,328,290,288};// 309, 288, 287, 315, 315, 336, 329, 328

	Far::Index tab_epi[11] = { 328,290,321 };// 309, 288, 287, 315, 315, 336, 329, 328};
	Far::ConstIndexArray array1;
	Far::ConstIndexArray array2;*/
	//Far::ConstIndexArray array1;


	//Vtr::internal::Level::VTag vtag = refiner->_levels[0]->getVertexTag((Far::Index)(105));

	// mitral 
	/*Vtr::internal::Level::VTag vtag = refiner->_levels[0]->getVertexTag((Far::Index)(20));
	refiner->_levels[0]->_vertTags[(Far::Index)(1)] = vtag;
	refiner->_levels[0]->_vertTags[(Far::Index)(2)] = vtag;
	refiner->_levels[0]->_vertTags[(Far::Index)(53)] = vtag;
	refiner->_levels[0]->_vertTags[(Far::Index)(3)] = vtag;
	refiner->_levels[0]->_vertTags[(Far::Index)(4)] = vtag;
	refiner->_levels[0]->_vertTags[(Far::Index)(56)] = vtag;
	refiner->_levels[0]->_vertTags[(Far::Index)(32)] = vtag;
	refiner->_levels[0]->_vertTags[(Far::Index)(67)] = vtag;
	refiner->_levels[0]->_vertTags[(Far::Index)(20)] = vtag;
	refiner->_levels[0]->_vertTags[(Far::Index)(5)] = vtag;
	refiner->_levels[0]->_vertTags[(Far::Index)(21)] = vtag;
	refiner->_levels[0]->_vertTags[(Far::Index)(69)] = vtag;

	refiner->_levels[0]->_vertTags[(Far::Index)(199)] = vtag;
	refiner->_levels[0]->_vertTags[(Far::Index)(206)] = vtag;
	refiner->_levels[0]->_vertTags[(Far::Index)(192)] = vtag;


	array1 = refiner->_levels[0]->getVertexEdges((Far::Index)(109));
	Vtr::internal::Level::ETag etag = refiner->_levels[0]->getEdgeTag(array1[0]);

	float v = refiner->_levels[0]->getVertexSharpness((Far::Index)(127));
	float e = refiner->_levels[0]->getEdgeSharpness((Far::Index)(array1[0]));

	for (int i = 0; i < refiner->_levels[0]->getNumVertices(); i++)
	{
		if (refiner->_levels[0]->getNumVertexEdges((Far::Index)(i)) > 3) // extraordinary vertex
		{
			refiner->_levels[0]->_vertTags[(Far::Index)(i)] = refiner->_levels[0]->getVertexTag((Far::Index)(127)); //205
			refiner->_levels[0]->setVertSmooth((Far::Index)(i), v); // set sharpness a zero

			array1 = refiner->_levels[0]->getVertexEdges((Far::Index)(i));
			refiner->_levels[0]->_edgeTags[(Far::Index)(array1[0])] = etag;
			refiner->_levels[0]->_edgeTags[(Far::Index)(array1[1])] = etag;
			refiner->_levels[0]->_edgeTags[(Far::Index)(array1[2])] = etag;
			refiner->_levels[0]->_edgeTags[(Far::Index)(array1[3])] = etag;

			refiner->_levels[0]->setEdgeSmooth(array1[0], e);
			refiner->_levels[0]->setEdgeSmooth(array1[1], e);
			refiner->_levels[0]->setEdgeSmooth(array1[2], e);
			refiner->_levels[0]->setEdgeSmooth(array1[3], e);

		}
	}*/

	int maxlevel = 2; // 2 subdivision levels


	refiner->RefineUniform(Far::TopologyRefiner::UniformOptions(maxlevel));   // Refine the topology RefineUniform
	cout << "refinement " << endl;
	//refiner->RefineAdaptive(Far::TopologyRefiner::AdaptiveOptions(maxlevel));

	// PatchTable
	// ---------------------
	// Generate a set of Far::PatchTable that we will use to evaluate the surface limit
	Far::PatchTableFactory::Options patchOptions;
	patchOptions.endCapType = Far::PatchTableFactory::Options::ENDCAP_BSPLINE_BASIS;

	Far::PatchTable const* patchTable = Far::PatchTableFactory::Create(*refiner, patchOptions);           // constructor
	cout << "fin de la patchTable factory " << endl;

	Far::PatchTable::PatchVertsTable const& table = patchTable->GetPatchControlVerticesTable();            // Indices of the control vertices of the patches

																										   // Compute the total number of points we need to evaluate patchtable. We use local points around extraordinary features.
	int nRefinerVertices = refiner->GetNumVerticesTotal();  // gives the total number of vertices in all levels

	int nLocalPoints = patchTable->GetNumLocalPoints();     // returns the number of points of the change of basis patches

															// Create a buffer to hold the position of the refined verts and local points, then copy the coarse positions at the beginning.
	std::vector<Vertex> verts(nRefinerVertices + nLocalPoints);

	cout << "nRefinerVertices = " << nRefinerVertices << endl;
	cout << "nLocalPoint = " << nLocalPoints << endl;

	//	static int const g_nverts = 338;
	//static int const g_nverts = 40;
	static int const g_nverts = 1152;// 210; // control points

	memcpy(&verts[0], g_verts, g_nverts * 3 * sizeof(float));

	// Adaptive refinement may result in fewer levels than maxIsolation.
	int nRefinedLevels = refiner->GetNumLevels();      // returns the number of refinement levels


													   // Interpolate vertex primvar data : they are the control vertices of the limit patches
	Vertex* src = &verts[0];    // verts includes the refined vertices and the local points

	for (int level = 1; level < nRefinedLevels; ++level)
	{
		Vertex* dst = src + refiner->GetLevel(level - 1).GetNumVertices();
		Far::PrimvarRefiner(*refiner).Interpolate(level, src, dst);     // apply vertex interpolation weights to the primvar buffer for a single level of refinement
		src = dst;
	}


	// Evaluate local points from interpolated vertex primvars.
	cout << "evaluate local point... " << endl;
	//patchTable->ComputeLocalPointValues(&verts[0], &verts[nRefinerVertices]);    // updates local point values based on the refined values src=&vert[0] --> Buffer with primvar data for the control vertices and refined vertices. dst=&verts[nRefinerVertices] --> Destination buffer for the computed local points

	Far::StencilTable const* stenciltab = patchTable->GetLocalPointStencilTable();   // Returns the stencil table to get change of basis patch points
																			 // Creation of a Far::PatchMap to help locating patches in the table

	ofstream mytable1;

	mytable1.open("C:\\Users\\charl\\Desktop\\stencil.obj");

	// Stencil Table without Patches

	//cout << stenciltab->GetNumControlVertices() <<  endl;
	//cout << stenciltab->GetNumStencils() << endl;

	// get local point stencil table
	// -----------------------------------


	const int num_vertices1 = stenciltab->GetNumControlVertices();
	const int num_stencil1 = stenciltab->GetNumStencils();
	cout << " stenciltab->GetNumStencils() = " << stenciltab->GetNumStencils() << endl;
	const float* weights1;
	Far::Stencil Stencil1;
	const Far::Index* indices1;
	std::cout << "stencil table" << endl;
	float** array1 = new float*[num_stencil1];
	for (int i = 0; i < num_stencil1; ++i)
	array1[i] = new float[num_vertices1];


	// initialization
	for (int i = 0; i < num_stencil1; ++i)
	for (int j = 0; j < num_vertices1; ++j)
	array1[i][j] = 0.0;


	// calcul uniquement des poids et des indices	cout << "creation of the table" << endl;
	for (int i = 0; i < stenciltab->GetNumStencils(); i++)
	{
	Stencil1 = stenciltab->GetStencil((Far::Index)(i));
	int size1 = Stencil1.GetSize();
	indices1 = Stencil1.GetVertexIndices();
	weights1 = Stencil1.GetWeights();


	for (int j = 0; j < size1; j++)
	mytable1 << weights1[j] << ' ';

	mytable1 << ' ';
	for (int k = 0; k < size1; k++)
	mytable1 << indices1[k] << ' ';
	mytable1 << endl;
	}


	// The PatchMap provides a quad-tree based lookup structure that, given a singular parametric location, can efficiently return a handle to the sub - patch that contains this location. 
	// (An quadtree-based map connecting coarse faces to their sub-patches).
	Far::PatchMap patchmap(*patchTable);

	// Create a Far::PtexIndices to help find indices of ptex faces.
	Far::PtexIndices ptexIndices(*refiner); // Object used to compute and query ptex face indices. Given a refiner, constructing a PtexIndices object builds the mapping from coarse faces 
											// to ptex ids.Once built, the object can be used to query the mapping.
											// Generate random samples on each ptex face
	int nsamples = 25;

	int	nfaces = ptexIndices.GetNumFaces();

	cout << "std vector " << endl;
	std::vector<LimitFrame> samples(nsamples * nfaces);

	srand(static_cast<int>(2147483647));

	//float pWeights[20], dsWeights[20], dtWeights[20];

	// patch parameters

	ofstream file;
	file.open("C:\\Users\\charl\\Desktop\\testMatlab.obj");
	ofstream patch;
	patch.open("C:\\Users\\charl\\Desktop\\patch_param.obj");

	int u;
	cout << "evaluation... " << endl;
	for (int face = 0, count = 0; face < nfaces; ++face)
	{
		cout << "faces =" << face;
		for (float s = 0.0f; s <= 1.0f; s = s + 0.25f)
		{
			cout << " s =" << s;
			for (float t = 0.0f; t <= 1.0f; t = t + 0.25f, ++count)
			{

				// Locate the patch corresponding to the face ptex idx and (s,t)
				Far::PatchTable::PatchHandle const* handle = patchmap.FindPatch(face, s, t);

				if (face < 96)
				{
					u = 0;
				}
				else
				{
					u = 1;
				}
				//file << face + 1 << " " << s << " " << t << " " << u << " " << handle->patchIndex +1 << endl;
				Far::PatchParam const& param = patchTable->_paramTable[handle->patchIndex];
				assert(handle);

				// Evaluate the patch weights, identify the CVs and compute the limit frame:
				//patchTable->EvaluateBasis(*handle, s, t, pWeights, dsWeights, dtWeights);
				Far::ConstIndexArray cvs = patchTable->GetPatchVertices(*handle);  // Returns the control vertex indices for the patch identified by 'handle'

																				   //patch << handle->patchIndex << " ";
				//LimitFrame & dst = samples[count];
				//dst.Clear();
				file << face + 1 << " " << handle->patchIndex << " ";
				patch << handle->patchIndex << " ";
				for (int cv = 0; cv < 16; cv++)
				{
					file << cvs[cv] + 1 << " ";
					patch << cvs[cv] << " ";
					//dst.AddWithWeight(verts[cvs[cv]], pWeights[cv], dsWeights[cv], dtWeights[cv]);
				}

				// passe en coordonee patch
				float frac = param.GetParamFraction();

				// top left corner
				float pu = (float)param.GetU() * frac;
				float pv = (float)param.GetV() * frac;

				// normalize u,v coordinates
				float s_patch = (s - pu) / frac;
				float t_patch = (t - pv) / frac;

				file << param.GetBoundary() << " " << s_patch << " " << t_patch << " " << u << " " << s << " " << t << " " << u << " " << param.GetParamFraction() << " " << param.GetU() << " " << param.GetV() << endl;
				patch << endl;
			}
		}
	}
	cout << "end parameters" << endl;
	//  End Patch
	// __________________________


	cout << "Stencil Table" << endl;
	Far::StencilTableFactory::Options options;
	options.generateIntermediateLevels = false;// false;
	options.generateOffsets = true;
	Far::StencilTable const* stencilTable = Far::StencilTableFactory::Create(*refiner, options);

	ofstream mytable;
	mytable.open("C:\\Users\\charl\\Desktop\\StencilTableessai.obj");

	// Stencil Table without Patches

	const int num_vertices = stencilTable->GetNumControlVertices();
	const int num_stencil = stencilTable->GetNumStencils();
	const float* weights;
	Far::Stencil Stencil;
	const Far::Index* indices;
	std::cout << "stencil table" << endl;
	float** array = new float* [num_stencil];
	for (int i = 0; i < num_stencil; ++i)
		array[i] = new float[num_vertices];


	// initialization
	for (int i = 0; i < num_stencil; ++i)
		for (int j = 0; j < num_vertices; ++j)
			array[i][j] = 0.0;


	// calcul uniquement des poids et des indices	cout << "creation of the table" << endl;
	/*for (int i = 0; i < stencilTable->GetNumStencils(); i++)
	{
	Stencil = stencilTable->GetStencil((Far::Index)(i));
	int size = Stencil.GetSize();
	indices = Stencil.GetVertexIndices();
	weights = Stencil.GetWeights();


	for (int j = 0; j < size; j++)
	mytable << weights[j] << ' ';

	mytable << ' ' ;
	for (int k = 0; k < size; k++)
	mytable << indices[k] << ' ';
	mytable << endl;
	}*/


	for (int i = 0; i < stencilTable->GetNumStencils(); i++)
	{
		Stencil = stencilTable->GetStencil((Far::Index)(i));
		int size = Stencil.GetSize();
		indices = Stencil.GetVertexIndices();
		weights = Stencil.GetWeights();

		for (int j = 0; j < size; j++)
		{
			int m = (int)(indices[j]);
			array[i][m] = weights[j];
		}

	}

	for (int i = 0; i < stencilTable->GetNumStencils(); i++)
	{
		for (int j = 0; j < num_vertices; j++)
		{
			mytable << array[i][j] << ' ';
		}

		mytable << endl;

	}
	mytable.close();

	// Allocate vertex primvar buffer (1 stencil for each vertex)
	cout << "Allocate vertex primvar buffer" << endl;
	int nstencils = stencilTable->GetNumStencils();
	std::vector<Vertex> vertexBuffer(nstencils);

	Vertex* controlValues = reinterpret_cast<Vertex*>(g_verts);

	{
		stencilTable->UpdateValues(controlValues, &vertexBuffer[0]);             //Updates point values based on the control values
	}

	cout << "ecriture..." << endl;
	if (myfile.is_open())
	{
		Far::TopologyLevel const& refLastLevel = refiner->GetLevel(maxlevel);
		int nfaces = refLastLevel.GetNumFaces();

		// print the vertex position
		for (int i = 0; i < (int)vertexBuffer.size(); ++i)
		{
			float const* pos = vertexBuffer[i].GetPosition();

			myfile << "v " << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
		}

		// Print faces
		for (int face = 0; face < nfaces; ++face)
		{

			Far::ConstIndexArray fverts = refLastLevel.GetFaceVertices(face);


			myfile << "f ";
			for (int vert = 0; vert < fverts.size(); ++vert)
			{
				myfile << fverts[vert] + 1 << " "; // OBJ uses 1-based arrays...
			}
			myfile << "\n";
		}
		myfile.close();
	}
	else std::cout << "Unable to open file";


	std::printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);


	delete refiner;
	delete stencilTable;

}

//------------------------------------------------------------------------------

Far::TopologyRefiner * createTopologyRefiner(float *g_verts)
{
	// Mesh creation
	ifstream infile;
	ofstream infile1;

	//	Read the file
	infile.open("C:\\Users\\charl\\Desktop\\bi-a_model.obj");// file containing numbers in 3 columns
	infile1.open("C:\\Users\\charl\\Desktop\\verif.obj");// file containing numbers in 3 columns

	int pos = 0;
	int pos1 = 0;
	int pos2 = 0;
	int line = 0;
	int g_nverts = 0;
	int g_nfaces = 0;
	int g_vertsperface[10096];
	int g_vertIndices[10096];

	while (!infile.eof())
	{
		string line2;
		string temp;
		while (getline(infile, line2))
		{
			// counter le nombre colonne
			int ncols = 0;
			// calcul le nombre de colonne:

			istringstream iss(line2);
			do
			{
				string sub;
				iss >> sub;
				if (sub.length())
					ncols++;
			} while (iss);
			g_vertsperface[pos2] = ncols;
			pos2++;

			if (ncols == 5)       // so it is a face with 4 indices
			{
				istringstream ss(line2);
				char c;
				ss >> c;
				int num1, num2, num3, num4;
				while (ss >> num1 >> num2 >> num3 >> num4)
				{
					g_vertIndices[pos1] = num1 - 1;                  // les indices de .obj crees par Blender comment par 1 et non par 0. Or la fonction veut qu'ils commencent par 1
					pos1++;
					g_vertIndices[pos1] = num2 - 1;
					pos1++;
					g_vertIndices[pos1] = num3 - 1;
					pos1++;
					g_vertIndices[pos1] = num4 - 1;
					pos1++;
					infile1 << "f " << num1 << " " << num2 << " " << num3 << " " << num4 << endl;
					g_nfaces++;
				}
			}

			if (ncols == 4)       // could be either a vertex or face with 3 indices
			{
				istringstream ss(line2);
				char c;
				ss >> c;

				if (c == 'v')     // This is a vertex point
				{
					float num1, num2, num3;
					while (ss >> num1 >> num2 >> num3)
					{
						g_verts[pos] = num1;
						pos++;
						g_verts[pos] = num2;
						pos++;
						g_verts[pos] = num3;
						pos++;
						infile1 << "v " << num1 << " " << num2 << " " << num3 << "" << endl;
						g_nverts++;
					}
				}

				if (c == 'f')     // This is a face
				{
					int num1, num2, num3;
					while (ss >> num1 >> num2 >> num3)
					{
						g_vertIndices[pos1] = num1 - 1;
						pos1++;
						g_vertIndices[pos1] = num2 - 1;
						pos1++;
						g_vertIndices[pos1] = num3 - 1;
						pos1++;
						infile1 << "f " << num1 << " " << num2 << " " << num3 << endl;
						g_nfaces++;
					}
				}


			}
			line++;
		}
	}

	// Populate a topology descriptor with our raw data.
	typedef Far::TopologyDescriptor Descriptor;
	Sdc::SchemeType type = OpenSubdiv::Sdc::SCHEME_CATMARK;
	Sdc::Options options;
	options.SetVtxBoundaryInterpolation(Sdc::Options::VTX_BOUNDARY_EDGE_ONLY);
	Descriptor desc;

	desc.numVertices = g_nverts;
	desc.numFaces = g_nfaces;
	desc.numVertsPerFace = g_vertsperface;
	desc.vertIndicesPerFace = g_vertIndices;

	// Instantiate a FarTopologyRefiner from the descriptor.
	return Far::TopologyRefinerFactory<Descriptor>::Create(desc,
		Far::TopologyRefinerFactory<Descriptor>::Options(type, options));

}

