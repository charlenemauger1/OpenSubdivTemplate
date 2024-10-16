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
#ifndef OPENSUBDIV3_FAR_TOPOLOGY_REFINER_H
#define OPENSUBDIV3_FAR_TOPOLOGY_REFINER_H

#include "../version.h"

#include "../sdc/types.h"
#include "../sdc/options.h"
#include "../far/types.h"
#include "../far/topologyLevel.h"

#include <vector>


namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

namespace Vtr { namespace internal { class SparseSelector; } }

namespace Far {

template <class MESH> class TopologyRefinerFactory;

///
///  \brief Stores topology data for a specified set of refinement options.
///
class TopologyRefiner {

public:

    /// \brief Constructor
    TopologyRefiner(Sdc::SchemeType type, Sdc::Options options = Sdc::Options());

    /// \brief Destructor
    ~TopologyRefiner();

    /// \brief Returns the subdivision scheme
    Sdc::SchemeType GetSchemeType() const    { return _subdivType; }

    /// \brief Returns the subdivision options
    Sdc::Options GetSchemeOptions() const { return _subdivOptions; }

    /// \brief Returns true if uniform subdivision has been applied
    bool IsUniform() const   { return _isUniform; }

    /// \brief Returns the number of refinement levels
    int  GetNumLevels() const { return (int)_levels.size(); }

    /// \brief Returns the highest level of refinement
    int  GetMaxLevel() const { return _maxLevel; }

    /// \brief Returns the maximum vertex valence in all levels
    int  GetMaxValence() const { return _maxValence; }

    /// \ brief Returns true if faces have been tagged as holes
    bool HasHoles() const { return _hasHoles; }

    /// \brief Returns the total number of vertices in all levels
    int GetNumVerticesTotal() const { return _totalVertices; }

    /// \brief Returns the total number of edges in all levels
    int GetNumEdgesTotal() const { return _totalEdges; }

    /// \brief Returns the total number of edges in all levels
    int GetNumFacesTotal() const { return _totalFaces; }

    /// \brief Returns the total number of face vertices in all levels
    int GetNumFaceVerticesTotal() const { return _totalFaceVertices; }

    /// \brief Returns a handle to access data specific to a particular level
    TopologyLevel const & GetLevel(int level) const { return _farLevels[level]; }

    //@{
    ///  @name High-level refinement and related methods
    ///

    //
    // Uniform refinement
    //

    /// \brief Uniform refinement options
    struct UniformOptions {

        UniformOptions(int level) :
            refinementLevel(level),
            orderVerticesFromFacesFirst(false),
            fullTopologyInLastLevel(false) { }

        unsigned int refinementLevel:4,             ///< Number of refinement iterations
                     orderVerticesFromFacesFirst:1, ///< Order child vertices from faces first
                                                    ///< instead of child vertices of vertices
                     fullTopologyInLastLevel:1;     ///< Skip topological relationships in the last
                                                    ///< level of refinement that are not needed for
                                                    ///< interpolation (keep false if using limit).
    };

    /// \brief Refine the topology uniformly
    ///
    /// @param options   Options controlling uniform refinement
    ///
    void RefineUniform(UniformOptions options);

    /// \brief Returns the options specified on refinement
    UniformOptions GetUniformOptions() const { return _uniformOptions; }

    //
    // Adaptive refinement
    //

    /// \brief Adaptive refinement options
    struct AdaptiveOptions {

        AdaptiveOptions(int level) :
            isolationLevel(level),
            useSingleCreasePatch(false),
            orderVerticesFromFacesFirst(false) { }

        unsigned int isolationLevel:4,              ///< Number of iterations applied to isolate
                                                    ///< extraordinary vertices and creases
                     useSingleCreasePatch:1,        ///< Use 'single-crease' patch and stop
                                                    ///< isolation where applicable
                     orderVerticesFromFacesFirst:1; ///< Order child vertices from faces first
                                                    ///< instead of child vertices of vertices
    };

    /// \brief Feature Adaptive topology refinement (restricted to scheme Catmark)
    ///
    /// @param options   Options controlling adaptive refinement
    ///
    void RefineAdaptive(AdaptiveOptions options);

    /// \brief Returns the options specified on refinement
    AdaptiveOptions GetAdaptiveOptions() const { return _adaptiveOptions; }

    /// \brief Unrefine the topology (keep control cage)
    void Unrefine();


    //@{
    /// @name Number and properties of face-varying channels:
    ///

    /// \brief Returns the number of face-varying channels in the tables
    int GetNumFVarChannels() const;

    /// \brief Returns the face-varying interpolation rule-set for a given channel
    Sdc::Options::FVarLinearInterpolation GetFVarLinearInterpolation(int channel = 0) const;

    /// \brief Returns the total number of face-varying values in all levels
    int GetNumFVarValuesTotal(int channel = 0) const;

    //@}

protected:

    //
    //  Lower level protected methods intended strictly for internal use:
    //
    template <class MESH>
    friend class TopologyRefinerFactory;
    friend class TopologyRefinerFactoryBase;
    friend class PatchTableFactory;
    friend class EndCapGregoryBasisPatchFactory;
    friend class EndCapLegacyGregoryPatchFactory;
    friend class PtexIndices;
    friend class PrimvarRefiner;

    Vtr::internal::Level & getLevel(int l) { return *_levels[l]; }
    Vtr::internal::Level const & getLevel(int l) const { return *_levels[l]; }

    Vtr::internal::Refinement & getRefinement(int l) { return *_refinements[l]; }
    Vtr::internal::Refinement const & getRefinement(int l) const { return *_refinements[l]; }

private:
    //  Not default constructible or copyable:
    TopologyRefiner() : _uniformOptions(0), _adaptiveOptions(0) { }
    TopologyRefiner(TopologyRefiner const &) : _uniformOptions(0), _adaptiveOptions(0) { }
    TopologyRefiner & operator=(TopologyRefiner const &) { return *this; }

    void selectFeatureAdaptiveComponents(Vtr::internal::SparseSelector& selector);

    void initializeInventory();
    void updateInventory(Vtr::internal::Level const & newLevel);

    void appendLevel(Vtr::internal::Level & newLevel);
    void appendRefinement(Vtr::internal::Refinement & newRefinement);
    void assembleFarLevels();

public:

    Sdc::SchemeType _subdivType;
    Sdc::Options    _subdivOptions;

    unsigned int _isUniform : 1,
                 _hasHoles : 1,
                 _maxLevel : 4;

    //  Options assigned on refinement:
    UniformOptions  _uniformOptions;
    AdaptiveOptions _adaptiveOptions;

    //  Cumulative properties of all levels:
    int _totalVertices;
    int _totalEdges;
    int _totalFaces;
    int _totalFaceVertices;
    int _maxValence;

    //  There is some redundancy here -- to be reduced later
    std::vector<Vtr::internal::Level *>      _levels;
    std::vector<Vtr::internal::Refinement *> _refinements;

    std::vector<TopologyLevel> _farLevels;
};


inline int
TopologyRefiner::GetNumFVarChannels() const {

    return _levels[0]->getNumFVarChannels();
}
inline Sdc::Options::FVarLinearInterpolation
TopologyRefiner::GetFVarLinearInterpolation(int channel) const {

    return _levels[0]->getFVarOptions(channel).GetFVarLinearInterpolation();
}

} // end namespace Far

} // end namespace OPENSUBDIV_VERSION
using namespace OPENSUBDIV_VERSION;
} // end namespace OpenSubdiv

#endif /* OPENSUBDIV3_FAR_TOPOLOGY_REFINER_H */
