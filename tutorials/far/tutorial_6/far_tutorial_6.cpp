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

#include <opensubdiv/far/stencilTable.h>
#include <opensubdiv/far/stencilTableFactory.h>

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

// Creates a Far::TopologyRefiner from the pyramid shape above
static Far::TopologyRefiner * createTopologyRefiner();

//------------------------------------------------------------------------------
// Vertex container implementation.
//
struct Vertex {

	// Minimal required interface ----------------------
	Vertex() { }

	void Clear(void * = 0) {
		_position[0] = _position[1] = _position[2] = 0.0f;
	}

	void AddWithWeight(Vertex const & src, float weight) {
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

//------------------------------------------------------------------------------
// Limit frame container implementation -- this interface is not strictly
// required but follows a similar pattern to Vertex.
//
struct LimitFrame {

	void Clear(void * = 0) {
		_position[0] = _position[1] = _position[2] = 0.0f;
		deriv1[0] = deriv1[1] = deriv1[2] = 0.0f;
		deriv2[0] = deriv2[1] = deriv2[2] = 0.0f;
	}

	void AddWithWeight(Vertex const & src,
		float weight, float d1Weight, float d2Weight) {

		_position[0] += weight * src._position[0];
		_position[1] += weight * src._position[1];
		_position[2] += weight * src._position[2];

		deriv1[0] += d1Weight * src._position[0];
		deriv1[1] += d1Weight * src._position[1];
		deriv1[2] += d1Weight * src._position[2];

		deriv2[0] += d2Weight * src._position[0];
		deriv2[1] += d2Weight * src._position[1];
		deriv2[2] += d2Weight * src._position[2];
	}

	float _position[3],
		deriv1[3],
		deriv2[3];
};

void export_subdivived_mesh(Far::TopologyRefiner* refiner, const char* output_obj, float* g_verts);

Far::TopologyRefiner* createTopologyRefiner(float* g_verts, const char* input_file);
void export_subdivived_matrices(Far::TopologyRefiner* refiner, float* g_verts, const char* output_sub_matrix, bool generate_intermediate_levels, bool generate_offset);

//------------------------------------------------------------------------------
int main(int, char **) {

	// Generate a FarTopologyRefiner (see far_tutorial_0 for details).
	//Far::TopologyRefiner * refiner = createTopologyRefiner();

	const char* input_file = "C:\\Users\\charl\\Desktop\\bi-a_model.obj";
	const char* output_mesh = "C:\\Users\\charl\\Desktop\\bi-a_model-sub.obj";
	const char* subdivision_matrix = "C:\\Users\\charl\\Desktop\\bi-a-template\\subdivision_matrix.txt";
	float g_verts[100000];
	Far::TopologyRefiner* refiner = createTopologyRefiner(g_verts, input_file);

	// generate subdivision matrix
	Far::TopologyRefiner::UniformOptions options(2);
	refiner->RefineUniform(options);   // Refine the topology RefineUniform
	export_subdivived_mesh(refiner, output_mesh, g_verts);
	bool generate_offset = true;
	bool generate_intermediate_levels = false; // we want the final matrix
	//bool factorizeIntermediateLevels = true;
	export_subdivived_matrices(refiner, g_verts, subdivision_matrix, generate_intermediate_levels, generate_offset);

	// ADAPTIVE REFINEMENT
	refiner = createTopologyRefiner(g_verts, input_file);
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

	std::cout << "nRefinerVertices = " << nRefinerVertices << endl;
	std::cout << "nLocalPoint = " << nLocalPoints << endl;
	std::cout << "n control points = " << g_nverts << endl;
	std::cout << "nRefinedLevels = " << nRefinedLevels << endl;
	std::cout << "total for basis calculation " << nRefinerVertices + nLocalPoints << endl;

	// Interpolate vertex primvar data : they are the control vertices
	// of the limit patches (see far_tutorial_0 for details)
	Vertex * src = &verts[0];
	for (int level = 1; level < nRefinedLevels; ++level) {
		Vertex * dst = src + refiner->GetLevel(level - 1).GetNumVertices();
		Far::PrimvarRefiner(*refiner).Interpolate(level, src, dst);
		src = dst;
	}

	std::cout << "Evaluate local points from interpolated vertex primvars" << endl;
	// Evaluate local points from interpolated vertex primvars.
	patchTable->ComputeLocalPointValues(&verts[0], &verts[nRefinerVertices]);

	// get local vpoints basis
	Far::StencilTable const* stenciltab = patchTable->GetLocalPointStencilTable();   // Returns the stencil table to get change of basis patch points
	//const int num_vertices1 = stenciltab->GetNumControlVertices();
	//const int num_stencil1 = stenciltab->GetNumStencils();
	//std::cout << "Export local points sub matrix " << num_vertices1 << " " << num_stencil1 << endl;

	std::cout << "Create a Far::PatchMap to help locating patches in the table" << endl;
	// Create a Far::PatchMap to help locating patches in the table
	Far::PatchMap patchmap(*patchTable);

	std::cout << "Create a Far::PtexIndices to help find indices of ptex faces" << endl;
	// Create a Far::PtexIndices to help find indices of ptex faces.
	Far::PtexIndices ptexIndices(*refiner);

	// Generate random samples on each ptex face
	int nsamples = 25;// 200,
	int nfaces = ptexIndices.GetNumFaces();

	std::vector<LimitFrame> samples(nsamples * nfaces);
	srand(static_cast<int>(2147483647));

	float pWeights[20], dsWeights[20], dtWeights[20];

	ofstream file, patch, boundary, bspline, etvertexelemnum, etvertexXi, fraction, patch_coordinates, patch_index;
	file.open("C:\\Users\\charl\\Desktop\\bi-a-template\\testMatlab.txt");
	patch.open("C:\\Users\\charl\\Desktop\\bi-a-template\\patch_param.txt");
	boundary.open("C:\\Users\\charl\\Desktop\\bi-a-template\\boundary.txt");
	bspline.open("C:\\Users\\charl\\Desktop\\bi-a-template\\control_points_patches.txt");
	etvertexelemnum.open("C:\\Users\\charl\\Desktop\\bi-a-template\\etVertexElementNum.txt");
	etvertexXi.open("C:\\Users\\charl\\Desktop\\bi-a-template\\etVertexXi.txt");
	fraction.open("C:\\Users\\charl\\Desktop\\bi-a-template\\fraction.txt");
	patch_coordinates.open("C:\\Users\\charl\\Desktop\\bi-a-template\\patch_coordinates.txt");
	patch_index.open("C:\\Users\\charl\\Desktop\\bi-a-template\\patch_index.txt");

	int u = 0; //WHEN ADDIND EPICARDIUM

	for (int face = 0, count = 0; face<nfaces; ++face) {

		//for (int sample = 0; sample<nsamples; ++sample, ++count) {
		for (float s = 0; s <= 1.0f; s = s + 0.25f) {
			for (float t = 0.0f; t <= 1.0f; t = t + 0.25f, ++count) {

				// Locate the patch corresponding to the face ptex idx and (s,t)
				Far::PatchTable::PatchHandle const * handle = patchmap.FindPatch(face, s, t);
			
				Far::PatchParam const& param = patchTable->_paramTable[handle->patchIndex];
				assert(handle);

				// Evaluate the patch weights, identify the CVs and compute the limit frame:
				patchTable->EvaluateBasis(*handle, s, t, pWeights, dsWeights, dtWeights);

				Far::ConstIndexArray cvs = patchTable->GetPatchVertices(*handle);

				file << face + 1 << " " << handle->patchIndex << " ";
				patch << handle->patchIndex << " ";

				LimitFrame & dst = samples[count];
				dst.Clear();
				for (int cv = 0; cv < cvs.size(); ++cv) {
					file << cvs[cv] + 1 << " ";
					patch << cvs[cv] << " ";
					bspline << cvs[cv] + 1 << " ";
					dst.AddWithWeight(verts[cvs[cv]], pWeights[cv], dsWeights[cv], dtWeights[cv]);
				}
				bspline << endl;

				// passe en coordonee patch
				float frac = param.GetParamFraction();

				// top left corner
				float pu = (float)param.GetU() * frac;
				float pv = (float)param.GetV() * frac;

				// normalize u,v coordinates
				float s_patch = (s - pu) / frac;
				float t_patch = (t - pv) / frac;

				patch_index << handle->patchIndex + 1 << endl;
				etvertexelemnum << face + 1 << " " << endl;
				boundary << param.GetBoundary() << endl;
				patch_coordinates << s_patch << " " << t_patch << " " << u << endl;
				etvertexXi << s << " " << t << " " << u << endl;
				fraction << param.GetParamFraction() << endl;

				file << param.GetBoundary() << " " << s_patch << " " << t_patch << " " << u << " " << s << " " << t << " " << u << " " << param.GetParamFraction() << " " << param.GetU() << " " << param.GetV() << endl;
				patch << endl;

			}
		}
	}

	std::cout << "Writing the obj" << endl;
	{ // Visualization with Maya : print a MEL script that generates particles
	  // at the location of the limit vertices
		ofstream myfile;
		myfile.open("C:\\Users\\charl\\Desktop\\evaluated_points.obj");

		int nsamples = (int)samples.size();
		std::cout << nsamples << endl;
		// Output particle positions for the tangent
		for (int sample = 0; sample<nsamples; ++sample) {
			float const * pos = samples[sample]._position;
			myfile << "v " << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
		}
	}

}

//------------------------------------------------------------------------------
void export_subdivived_mesh(Far::TopologyRefiner* refiner, const char* output_obj, float* g_verts)
{
	ofstream output_file;
	output_file.open(output_obj);

	Far::StencilTableFactory::Options options;
	options.generateIntermediateLevels = false;
	options.generateOffsets = true;

	Far::StencilTable const* stencil_table = Far::StencilTableFactory::Create(*refiner, options);

	// export mesh
	// Allocate vertex primvar buffer (1 stencil for each vertex)
	int nstencils = stencil_table->GetNumStencils();
	std::vector<Vertex> vertexBuffer(nstencils);

	Vertex* controlValues = reinterpret_cast<Vertex*>(g_verts);
	{
		stencil_table->UpdateValues(controlValues, &vertexBuffer[0]);             //Updates point values based on the control values
	}

	if (output_file.is_open())
	{
		Far::TopologyLevel const& last_refinement_level = refiner->GetLevel(refiner->GetMaxLevel());
		int nfaces = last_refinement_level.GetNumFaces();

		// print the vertex position
		for (int i = 0; i < (int)vertexBuffer.size(); ++i)
		{
			float const* pos = vertexBuffer[i].GetPosition();

			output_file << "v " << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
		}

		// faces
		Far::ConstIndexArray fverts;
		for (int face = 0; face < nfaces; ++face)
		{
			fverts = last_refinement_level.GetFaceVertices(face);
			output_file << "f ";
			for (int vert = 0; vert < fverts.size(); ++vert)
			{
				output_file << fverts[vert] + 1 << " "; // OBJ uses 1-based arrays...
			}
			output_file << "\n";
		}
		output_file.close();
	}
	else std::cout << "Unable to open file";
}

Far::TopologyRefiner * createTopologyRefiner(float *g_verts, const char* input_file)
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

void export_subdivived_matrices(Far::TopologyRefiner* refiner, float* g_verts, const char* output_sub_matrix, bool generate_intermediate_levels, bool generate_offset)
{
	ofstream matrix_file;
	matrix_file.open(output_sub_matrix);

	Far::StencilTableFactory::Options options;
	options.generateIntermediateLevels = generate_intermediate_levels;
	options.generateOffsets = generate_offset;
	Far::StencilTable const* stencil_table = Far::StencilTableFactory::Create(*refiner, options);

	if (matrix_file.is_open())
	{
		// export matrix
		const int num_vertices = stencil_table->GetNumControlVertices();
		const int num_stencil = stencil_table->GetNumStencils();
		const float* weights;
		Far::Stencil stencil;
		const Far::Index* indices;
		std::cout << "Export subdivision matrix " << num_vertices << " " << num_stencil << endl;
		float** array = new float*[num_stencil];
		for (int i = 0; i < num_stencil; ++i)
			array[i] = new float[num_vertices];

		// initialization
		for (int i = 0; i < num_stencil; ++i)
			for (int j = 0; j < num_vertices; ++j)
				array[i][j] = 0.0;

		for (int i = 0; i < num_stencil; i++)
		{
			stencil = stencil_table->GetStencil((Far::Index)(i));
			int size = stencil.GetSize();
			indices = stencil.GetVertexIndices();
			weights = stencil.GetWeights();

			for (int j = 0; j < size; j++)
			{
				int m = (int)(indices[j]);
				array[i][m] = weights[j];
			}
		}

		for (int i = 0; i < stencil_table->GetNumStencils(); i++)
		{
			for (int j = 0; j < num_vertices; j++)
			{
				matrix_file << array[i][j] << ' ';
			}

			matrix_file << endl;

		}
		matrix_file.close();
	}
	else std::cout << "Unable to open file";
}