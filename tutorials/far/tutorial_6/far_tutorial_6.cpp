//------------------------------------------------------------------------------
// Tutorial description:
//
// This tutorial shows how to interpolate surface limits at arbitrary
// parametric locations using feature adaptive Far::PatchTables.
//
// The evaluation of the limit surface at arbitrary locations requires the
// adaptive isolation of topological features. This process converts the
// input polygonal control cage into a collection of bi-cubic patches.
//
// We can then evaluate the patches at random parametric locations and
// obtain analytical positions and tangents on the limit surface.
//
// The results are dumped into a MEL script that draws 'streak' particle
// systems that show the tangent and bi-tangent at the random samples locations.
//

#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/far/primvarRefiner.h>
#include <opensubdiv/far/patchTableFactory.h>
#include <opensubdiv/far/patchMap.h>
#include <opensubdiv/far/ptexIndices.h>

#include <cassert>
#include <cstdio>
#include <cstring>
#include <cfloat>

using namespace OpenSubdiv;

// pyramid geometry from catmark_pyramid_crease0.h
//static int const g_nverts = 5;
//static float const g_verts[24] = { 0.0f,  0.0f, 2.0f,
//0.0f, -2.0f, 0.0f,
//2.0f,  0.0f, 0.0f,
//0.0f,  2.0f, 0.0f,
//-2.0f,  0.0f, 0.0f, };
//
//
//static int const g_vertsperface[5] = { 3, 3, 3, 3, 4 };
//
//static int const g_nfaces = 5;
//static int const g_faceverts[16] = { 0, 1, 2,
//0, 2, 3,
//0, 3, 4,
//0, 4, 1,
//4, 3, 2, 1 };

// Creates a Far::TopologyRefiner from the pyramid shape above
static Far::TopologyRefiner * createTopologyRefiner();

//------------------------------------------------------------------------------
// Vertex container implementation.
//
struct Vertex {

	// Minimal required interface ----------------------
	Vertex() { }

	void Clear(void * = 0) {
		point[0] = point[1] = point[2] = 0.0f;
	}

	void AddWithWeight(Vertex const & src, float weight) {
		point[0] += weight * src.point[0];
		point[1] += weight * src.point[1];
		point[2] += weight * src.point[2];
	}

	float point[3];
};

//------------------------------------------------------------------------------
// Limit frame container implementation -- this interface is not strictly
// required but follows a similar pattern to Vertex.
//
struct LimitFrame {

	void Clear(void * = 0) {
		point[0] = point[1] = point[2] = 0.0f;
		deriv1[0] = deriv1[1] = deriv1[2] = 0.0f;
		deriv2[0] = deriv2[1] = deriv2[2] = 0.0f;
	}

	void AddWithWeight(Vertex const & src,
		float weight, float d1Weight, float d2Weight) {

		point[0] += weight * src.point[0];
		point[1] += weight * src.point[1];
		point[2] += weight * src.point[2];

		deriv1[0] += d1Weight * src.point[0];
		deriv1[1] += d1Weight * src.point[1];
		deriv1[2] += d1Weight * src.point[2];

		deriv2[0] += d2Weight * src.point[0];
		deriv2[1] += d2Weight * src.point[1];
		deriv2[2] += d2Weight * src.point[2];
	}

	float point[3],
		deriv1[3],
		deriv2[3];
};

Far::TopologyRefiner* createTopologyRefiner_2(float* g_verts, const char* input_file);

//------------------------------------------------------------------------------
int main(int, char **) {

	// Generate a FarTopologyRefiner (see far_tutorial_0 for details).
	//Far::TopologyRefiner * refiner = createTopologyRefiner();

	const char* input_file = "C:\\Users\\charl\\Desktop\\bi-a_model.obj";

	float g_verts[100000];
	Far::TopologyRefiner* refiner = createTopologyRefiner_2(g_verts, input_file);

	// Adaptively refine the topology with an isolation level capped at 3
	// because the sharpest crease in the shape is 3.0f (in g_creaseweights[])
	int maxIsolation = 2;

	refiner->RefineAdaptive(Far::TopologyRefiner::AdaptiveOptions(maxIsolation));

	// Generate a set of Far::PatchTable that we will use to evaluate the
	// surface limit
	Far::PatchTableFactory::Options patchOptions;
	patchOptions.endCapType = Far::PatchTableFactory::Options::ENDCAP_BSPLINE_BASIS;

	Far::PatchTable const * patchTable = Far::PatchTableFactory::Create(*refiner, patchOptions);

	// Compute the total number of points we need to evaluate patchtable.
	// we use local points around extraordinary features.
	int nRefinerVertices = refiner->GetNumVerticesTotal();
	int nLocalPoints = patchTable->GetNumLocalPoints();

	// Create a buffer to hold the position of the refined verts and
	// local points, then copy the coarse positions at the beginning.
	std::vector<Vertex> verts(nRefinerVertices + nLocalPoints);

	static int const g_nverts = refiner->GetLevel(0).GetNumVertices(); // control points
	memcpy(&verts[0], g_verts, g_nverts * 3 * sizeof(float));

	// Adaptive refinement may result in fewer levels than maxIsolation.
	int nRefinedLevels = refiner->GetNumLevels();

	// Interpolate vertex primvar data : they are the control vertices
	// of the limit patches (see far_tutorial_0 for details)
	Vertex * src = &verts[0];
	for (int level = 1; level < nRefinedLevels; ++level) {
		Vertex * dst = src + refiner->GetLevel(level - 1).GetNumVertices();
		Far::PrimvarRefiner(*refiner).Interpolate(level, src, dst);
		src = dst;
	}

	cout << "Evaluate local points from interpolated vertex primvars" << endl;
	// Evaluate local points from interpolated vertex primvars.
	patchTable->ComputeLocalPointValues(&verts[0], &verts[nRefinerVertices]);

	cout << "Create a Far::PatchMap to help locating patches in the table" << endl;
	// Create a Far::PatchMap to help locating patches in the table
	Far::PatchMap patchmap(*patchTable);

	cout << "Create a Far::PtexIndices to help find indices of ptex faces" << endl;
	// Create a Far::PtexIndices to help find indices of ptex faces.
	Far::PtexIndices ptexIndices(*refiner);

	// Generate random samples on each ptex face
	int nsamples = 25;// 200,
	int nfaces = ptexIndices.GetNumFaces();

	std::vector<LimitFrame> samples(nsamples * nfaces);

	srand(static_cast<int>(2147483647));

	float pWeights[20], dsWeights[20], dtWeights[20];

	cout << "Beginning of the for loop" << endl;
	for (int face = 0, count = 0; face<nfaces; ++face) {

		//for (int sample = 0; sample<nsamples; ++sample, ++count) {
		for (float s = 0; s <= 1.0f; s = s + 0.25f) {
			for (float t = 0.0f; t <= 1.0f; t = t + 0.25f, ++count) {
			//float s = (float)rand() / (float)RAND_MAX,
			//	t = (float)rand() / (float)RAND_MAX;

			// Locate the patch corresponding to the face ptex idx and (s,t)
			Far::PatchTable::PatchHandle const * handle = patchmap.FindPatch(face, s, t);
			assert(handle);

			// Evaluate the patch weights, identify the CVs and compute the limit frame:
			patchTable->EvaluateBasis(*handle, s, t, pWeights, dsWeights, dtWeights);

			Far::ConstIndexArray cvs = patchTable->GetPatchVertices(*handle);

			LimitFrame & dst = samples[count];
			dst.Clear();
			for (int cv = 0; cv < cvs.size(); ++cv) {
				dst.AddWithWeight(verts[cvs[cv]], pWeights[cv], dsWeights[cv], dtWeights[cv]);
			}
		}

		}
	}

	cout << "writing..." << endl;
	{ // Visualization with Maya : print a MEL script that generates particles
	  // at the location of the limit vertices
		ofstream myfile;
		myfile.open("C:\\Users\\charl\\Desktop\\tangent.obj");

		int nsamples = (int)samples.size();
		std::cout << nsamples << endl;
		// Output particle positions for the tangent
		for (int sample = 0; sample<nsamples; ++sample) {
			float const * pos = samples[sample].point;
			myfile << "v " << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
		}
	}

}

//------------------------------------------------------------------------------
//static Far::TopologyRefiner *
//createTopologyRefiner() {
//
//
//	typedef Far::TopologyDescriptor Descriptor;
//	Sdc::SchemeType type = OpenSubdiv::Sdc::SCHEME_CATMARK;
//	Sdc::Options options;
//	options.SetVtxBoundaryInterpolation(Sdc::Options::VTX_BOUNDARY_EDGE_ONLY);
//	Descriptor desc;
//
//	desc.numVertices = g_nverts;
//	desc.numFaces = g_nfaces;
//	desc.numVertsPerFace = g_vertsperface;
//	desc.vertIndicesPerFace = g_faceverts;
//	//desc.numCreases = g_ncreases;
//	//desc.creaseVertexIndexPairs = g_creaseverts;
//	//desc.creaseWeights = g_creaseweights;
//
//	// Instantiate a FarTopologyRefiner from the descriptor.
//	Far::TopologyRefiner * refiner =
//		Far::TopologyRefinerFactory<Descriptor>::Create(desc,
//			Far::TopologyRefinerFactory<Descriptor>::Options(type, options));
//
//	return refiner;
//}

Far::TopologyRefiner * createTopologyRefiner_2(float *g_verts, const char* input_file)
{
	// Mesh creation
	ifstream infile;

	//	Read the file
	infile.open(input_file);// file containing numbers in 3 columns

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
					g_nfaces++;
				}

				g_vertsperface[pos2] = ncols - 1;
				pos2++;


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
						g_nfaces++;
					}
					g_vertsperface[pos2] = ncols - 1;
					pos2++;
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