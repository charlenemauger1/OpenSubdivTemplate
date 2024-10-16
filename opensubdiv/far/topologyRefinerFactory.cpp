//
//   Copyright 2014 DreamWorks Animation LLC.
//
//   Licensed under the Apache License, Version 2.0 (the "Apache License")
//   with the following modification; you may not use this file except in
//   compliance with the Apache License and the following modification to it:
//   Section 6. Trademarks. is deleted and replaced with:
//
//   6. Trademarks. This License does not grant permission to use the trade
//      names, trademarks, service marks, or product names of the Licensor
//      and its affiliates, except as required to comply with Section 4(c) of
//      the License and to reproduce the content of the NOTICE file.
//
//   You may obtain a copy of the Apache License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the Apache License with the above modification is
//   distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
//   KIND, either express or implied. See the Apache License for the specific
//   language governing permissions and limitations under the Apache License.
//
#include "../far/topologyRefinerFactory.h"
#include "../far/topologyRefiner.h"
#include "../sdc/types.h"
#include "../vtr/level.h"

#include <cstdio>
#ifdef _MSC_VER
    #define snprintf _snprintf
#endif


namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

namespace Far {

//
//  Methods for the Factory base class -- general enough to warrant including in
//  the base class rather than the subclass template (and so replicated for each
//  usage)
//
//
bool
TopologyRefinerFactoryBase::prepareComponentTopologySizing(TopologyRefiner& refiner) {

    Vtr::internal::Level& baseLevel = refiner.getLevel(0);

    //
    //  At minimum we require face-vertices (the total count of which can be determined
    //  from the offsets accumulated during sizing pass) and we need to resize members
    //  related to them to be populated during assignment:
    //
    int vCount = baseLevel.getNumVertices();
    int fCount = baseLevel.getNumFaces();

    if (vCount == 0) {
        Error(FAR_RUNTIME_ERROR, "Failure in TopologyRefinerFactory<>::Create() -- "
                "mesh contains no vertices.");
        return false;
    }
    if (fCount == 0) {
        Error(FAR_RUNTIME_ERROR, "Failure in TopologyRefinerFactory<>::Create() -- "
                "meshes without faces not yet supported.");
        return false;
    }

    //  Make sure no face was defined that would lead to a valence overflow -- the max
    //  valence has been initialized with the maximum number of face-vertices:
    if (baseLevel.getMaxValence() > Vtr::VALENCE_LIMIT) {
        char msg[1024];
        snprintf(msg, 1024, "Failure in TopologyRefinerFactory<>::Create() -- "
                "face with %d vertices > %d max.",
                baseLevel.getMaxValence(), Vtr::VALENCE_LIMIT);
        Error(FAR_RUNTIME_ERROR, msg);
        return false;
    }

    int fVertCount = baseLevel.getNumFaceVertices(fCount - 1) +
                     baseLevel.getOffsetOfFaceVertices(fCount - 1);

    if (fVertCount == 0) {
        Error(FAR_RUNTIME_ERROR, "Failure in TopologyRefinerFactory<>::Create() -- "
                "mesh contains no face-vertices.");
        return false;
    }
    if ((refiner.GetSchemeType() == Sdc::SCHEME_LOOP) && (fVertCount != (3 * fCount))) {
        Error(FAR_RUNTIME_ERROR, "Failure in TopologyRefinerFactory<>::Create() -- "
                "non-triangular faces not supported by Loop scheme.");
        return false;
    }
    baseLevel.resizeFaceVertices(fVertCount);

    //
    //  If edges were sized, all other topological relations must be sized with it, in
    //  which case we allocate those members to be populated.  Otherwise, sizing of the
    //  other topology members is deferred until the face-vertices are assigned and the
    //  resulting relationships determined:
    //
    int eCount = baseLevel.getNumEdges();

    if (eCount > 0) {
        baseLevel.resizeFaceEdges(baseLevel.getNumFaceVerticesTotal());
        baseLevel.resizeEdgeVertices();
        baseLevel.resizeEdgeFaces(  baseLevel.getNumEdgeFaces(eCount-1)   + baseLevel.getOffsetOfEdgeFaces(eCount-1));
        baseLevel.resizeVertexFaces(baseLevel.getNumVertexFaces(vCount-1) + baseLevel.getOffsetOfVertexFaces(vCount-1));
        baseLevel.resizeVertexEdges(baseLevel.getNumVertexEdges(vCount-1) + baseLevel.getOffsetOfVertexEdges(vCount-1));

        assert(baseLevel.getNumFaceEdgesTotal() > 0);
        assert(baseLevel.getNumEdgeVerticesTotal() > 0);
        assert(baseLevel.getNumEdgeFacesTotal() > 0);
        assert(baseLevel.getNumVertexFacesTotal() > 0);
        assert(baseLevel.getNumVertexEdgesTotal() > 0);
    }
    return true;
}

bool
TopologyRefinerFactoryBase::prepareComponentTopologyAssignment(TopologyRefiner& refiner, bool fullValidation,
                                                               TopologyCallback callback, void const * callbackData) {
    Vtr::internal::Level& baseLevel = refiner.getLevel(0);
    bool completeMissingTopology = (baseLevel.getNumEdges() == 0);

    if (completeMissingTopology) 
	{
        if (not baseLevel.completeTopologyFromFaceVertices()) 
		{
            char msg[1024];
            snprintf(msg, 1024, "Failure in TopologyRefinerFactory<>::Create() -- "
                    "vertex with valence %d > %d max.",
                    baseLevel.getMaxValence(), Vtr::VALENCE_LIMIT);
            Error(FAR_RUNTIME_ERROR, msg);
            return false;
        }
    } 
	else 
	{
        if (baseLevel.getMaxValence() == 0) 
		{
            Error(FAR_RUNTIME_ERROR, "Failure in TopologyRefinerFactory<>::Create() -- "
                "maximum valence not assigned.");
            return false;
        }
    }

	
    if (fullValidation) {
        if (not baseLevel.validateTopology(callback, callbackData)) {
            if (completeMissingTopology) {
                Error(FAR_RUNTIME_ERROR, "Failure in TopologyRefinerFactory<>::Create() -- "
                    "invalid topology detected from partial specification.");
            } else {
                Error(FAR_RUNTIME_ERROR, "Failure in TopologyRefinerFactory<>::Create() -- "
                    "invalid topology detected as fully specified.");
            }
            return false;
        }
    }

    //  Now that we have a valid base level, initialize the Refiner's component inventory:
    refiner.initializeInventory();
    return true;
}

bool
TopologyRefinerFactoryBase::prepareComponentTagsAndSharpness(TopologyRefiner& refiner) {

    //
    //  This method combines the initialization of internal component tags with the sharpening
    //  of edges and vertices according to the given boundary interpolation rule in the Options.
    //  Since both involve traversing the edge and vertex lists and noting the presence of
    //  boundaries -- best to do both at once...
    //
    Vtr::internal::Level&  baseLevel = refiner.getLevel(0);

    Sdc::Options options = refiner.GetSchemeOptions();
    Sdc::Crease  creasing(options);

    bool makeBoundaryFacesHoles = (options.GetVtxBoundaryInterpolation() == Sdc::Options::VTX_BOUNDARY_NONE
                                && Sdc::SchemeTypeTraits::GetLocalNeighborhoodSize(refiner.GetSchemeType()) > 0);
    bool sharpenCornerVerts     = (options.GetVtxBoundaryInterpolation() == Sdc::Options::VTX_BOUNDARY_EDGE_AND_CORNER);
    bool sharpenNonManFeatures  = true; //(options.GetNonManifoldInterpolation() == Sdc::Options::NON_MANIFOLD_SHARP);

    //
    //  Process the Edge tags first, as Vertex tags (notably the Rule) are dependent on
    //  properties of their incident edges.
    //
    for (Vtr::Index eIndex = 0; eIndex < baseLevel.getNumEdges(); ++eIndex) {
        Vtr::internal::Level::ETag& eTag       = baseLevel.getEdgeTag(eIndex);
        float&                      eSharpness = baseLevel.getEdgeSharpness(eIndex);

        eTag._boundary = (baseLevel.getNumEdgeFaces(eIndex) < 2);
        if (eTag._boundary || (eTag._nonManifold && sharpenNonManFeatures)) {
            eSharpness = Sdc::Crease::SHARPNESS_INFINITE;
        }
        eTag._infSharp  = Sdc::Crease::IsInfinite(eSharpness);
        eTag._semiSharp = Sdc::Crease::IsSharp(eSharpness) && !eTag._infSharp;
    }

    //
    //  Process the Vertex tags now -- for some tags (semi-sharp and its rule) we need
    //  to inspect all incident edges:
    //
    int schemeRegularInteriorValence = Sdc::SchemeTypeTraits::GetRegularVertexValence(refiner.GetSchemeType());
    int schemeRegularBoundaryValence = schemeRegularInteriorValence / 2;

    for (Vtr::Index vIndex = 0; vIndex < baseLevel.getNumVertices(); ++vIndex) {
        Vtr::internal::Level::VTag& vTag       = baseLevel.getVertexTag(vIndex);
        float&                      vSharpness = baseLevel.getVertexSharpness(vIndex);

        Vtr::ConstIndexArray vEdges = baseLevel.getVertexEdges(vIndex);
        Vtr::ConstIndexArray vFaces = baseLevel.getVertexFaces(vIndex);

        //
        //  Take inventory of properties of incident edges that affect this vertex:
        //
        int boundaryEdgeCount    = 0;
        int infSharpEdgeCount    = 0;
        int semiSharpEdgeCount   = 0;
        int nonManifoldEdgeCount = 0;
        for (int i = 0; i < vEdges.size(); ++i) {
            Vtr::internal::Level::ETag const& eTag = baseLevel.getEdgeTag(vEdges[i]);

            boundaryEdgeCount    += eTag._boundary;
            infSharpEdgeCount    += eTag._infSharp;
            semiSharpEdgeCount   += eTag._semiSharp;
            nonManifoldEdgeCount += eTag._nonManifold;
        }
        int sharpEdgeCount = infSharpEdgeCount + semiSharpEdgeCount;

        //
        //  Sharpen the vertex before using it in conjunction with incident edge
        //  properties to determine the semi-sharp tag and rule:
        //
        bool isTopologicalCorner = (vFaces.size() == 1) && (vEdges.size() == 2);
        bool isSharpenedCorner =  isTopologicalCorner && sharpenCornerVerts;
        if (isSharpenedCorner) {
            vSharpness = Sdc::Crease::SHARPNESS_INFINITE;
        } else if (vTag._nonManifold && sharpenNonManFeatures) {
            //
            //  We avoid sharpening non-manifold vertices when they occur on interior
            //  non-manifold creases, i.e. a pair of opposing non-manifold edges with
            //  more than two incident faces.  In these cases there are more incident
            //  faces than edges (1 more for each additional "fin") and no boundaries.
            //
            if (not ((nonManifoldEdgeCount == 2) && (boundaryEdgeCount == 0) && (vFaces.size() > vEdges.size()))) {
                vSharpness = Sdc::Crease::SHARPNESS_INFINITE;
            }
        }

        vTag._infSharp       = Sdc::Crease::IsInfinite(vSharpness);
        vTag._semiSharp      = Sdc::Crease::IsSemiSharp(vSharpness);
        vTag._semiSharpEdges = (semiSharpEdgeCount > 0);

        vTag._rule = (Vtr::internal::Level::VTag::VTagSize)creasing.DetermineVertexVertexRule(vSharpness, sharpEdgeCount);

        //
        //  Assign topological tags -- note that the "xordinary" tag is not strictly
        //  correct (or relevant) if non-manifold:
        //
        vTag._boundary = (boundaryEdgeCount > 0);
        vTag._corner = isSharpenedCorner;
        if (vTag._corner) {
            vTag._xordinary = false;
        } else if (vTag._boundary) {
            vTag._xordinary = (vFaces.size() != schemeRegularBoundaryValence);
        } else {
            vTag._xordinary = (vFaces.size() != schemeRegularInteriorValence);
        }
        vTag._incomplete = 0;

        //
        //  Having just decided if a vertex is on a boundary, and with its incident faces
        //  available, mark incident faces as holes.
        //
        if (makeBoundaryFacesHoles && vTag._boundary) {
            for (int i = 0; i < vFaces.size(); ++i) {
                baseLevel.getFaceTag(vFaces[i])._hole = true;

                //  Don't forget this -- but it will eventually move to the Level
                refiner._hasHoles = true;
            }
        }
    }
    return true;
}

bool
TopologyRefinerFactoryBase::prepareFaceVaryingChannels(TopologyRefiner& refiner) {

    Vtr::internal::Level& baseLevel = refiner.getLevel(0);

    int regVertexValence   = Sdc::SchemeTypeTraits::GetRegularVertexValence(refiner.GetSchemeType());
    int regBoundaryValence = regVertexValence / 2;

    for (int channel=0; channel<refiner.GetNumFVarChannels(); ++channel) {
        if (baseLevel.getNumFVarValues(channel) == 0) {
            char msg[1024];
            snprintf(msg, 1024, "Failure in TopologyRefinerFactory<>::Create() -- "
                    "face-varying channel %d has no values.", channel);
            Error(FAR_RUNTIME_ERROR, msg);
            return false;
        }
        baseLevel.completeFVarChannelTopology(channel, regBoundaryValence);
    }
    return true;
}

} // end namespace Far

} // end namespace OPENSUBDIV_VERSION
} // end namespace OpenSubdiv
