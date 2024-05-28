// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#pragma once
// Project include(s)

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/ProtoDetector.hpp"
//#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/TypeList.hpp"
#include "Acts/Utilities/GridAxisGenerators.hpp"
#include "Acts/Utilities/IAxis.hpp"


#include "Acts/Plugins/Json/AlgebraJsonConverter.hpp"
#include "Acts/Plugins/Json/SurfaceJsonConverter.hpp"
#include "Acts/Plugins/Json/DetrayJsonHelper.hpp"
#include "Acts/Plugins/Json/DetectorVolumeJsonConverter.hpp"

#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/DetectorVolumeUpdaters.hpp"

#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"

//debug
#include "detray/utils/type_list.hpp"


// detray geometry
#include "detray/builders/detector_builder.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/io/frontend/detector_reader_config.hpp"
#include "detray/io/frontend/implementation/json_readers.hpp"
#include "detray/io/frontend/utils/detector_components_reader.hpp"

#include "detray/io/common/geometry_reader.hpp"
#include "detray/io/common/geometry_writer.hpp"
#include "detray/io/common/surface_grid_writer.hpp"
#include "detray/io/common/surface_grid_reader.hpp"
#include "detray/io/common/detail/grid_reader.hpp"



#include "detray/io/json/json_reader.hpp"
#include "detray/io/json/json_writer.hpp"
#include "detray/io/json/json_io.hpp"
#include "detray/io/frontend/writer_interface.hpp"

// surface grid
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Detector/detail/IndexedSurfacesGenerator.hpp"

// checks - prints
#include "detray/utils/consistency_checker.hpp"
#include "detray/io/frontend/detector_reader.hpp"
#include "detray/navigation/volume_graph.hpp"

// System include(s)
#include <fstream>
#include <optional>
#include <initializer_list>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <ios>


namespace Acts {
class IMaterialDecorator;
}  // namespace Acts
namespace ActsExamples {
class IMaterialWriter;
class IWriter;
}  // namespace ActsExamples

using namespace Acts;
using namespace ActsExamples;
using namespace detray;

using detector_t = detector<default_metadata>;

//
//ACTS GEOMETRY TO DETRAY CONVERSION
//use acts detector data to fill in all the payloads and convert these to a detray detector

namespace {

/// Find the position of the volume to point to
///
/// @param volume the volume to find
/// @param the collection of volumes
///
/// @note return -1 if not found, to be interpreted by the caller
int findVolume(
    const Acts::Experimental::DetectorVolume* volume,
    const std::vector<const Acts::Experimental::DetectorVolume*>& volumes) {
  auto candidate = std::find(volumes.begin(), volumes.end(), volume);
  if (candidate != volumes.end()) {
    return std::distance(volumes.begin(), candidate);
  }
  return -1;
}
}  // namespace

namespace detray{

    /// These functions perform the conversion of the ACTS detector to the detray detector. 
    /// The payloads are created and populated by traversing the initial ACTS detector object.
        

    // debug print json differences
    void print_json_diff(const nlohmann::json& j1, const nlohmann::json& j2, const std::string& prefix = "") {
        // If j1 and j2 are not of the same type, they are different
        if (j1.type() != j2.type()) {
            std::cout << prefix << "Type mismatch: " << j1.type_name() << " vs " << j2.type_name() << std::endl;
            return;
        }

        // If j1 and j2 are objects, compare each key-value pair
        if (j1.is_object() && j2.is_object()) {
            for (auto it = j1.begin(); it != j1.end(); ++it) {
                std::string key = it.key();
                if (j2.find(key) != j2.end()) {
                    print_json_diff(it.value(), j2.at(key), prefix + key + ".");
                } else {
                    std::cout << prefix << key << " is missing in second JSON" << std::endl;
                }
            }
            for (auto it = j2.begin(); it != j2.end(); ++it) {
                std::string key = it.key();
                if (j1.find(key) == j1.end()) {
                    std::cout << prefix << key << " is missing in first JSON" << std::endl;
                }
            }
        }
        // If j1 and j2 are arrays, compare each element
        else if (j1.is_array() && j2.is_array()) {
            size_t size = std::max(j1.size(), j2.size());
            for (size_t i = 0; i < size; ++i) {
                if (i < j1.size() && i < j2.size()) {
                    print_json_diff(j1[i], j2[i], prefix + "[" + std::to_string(i) + "].");
                } else if (i < j1.size()) {
                    std::cout << prefix << "[" << i << "] is missing in second JSON" << std::endl;
                } else {
                    std::cout << prefix << "[" << i << "] is missing in first JSON" << std::endl;
                }
            }
        }
        // If j1 and j2 are primitives (string, number, boolean, or null), compare their values
        else if (j1 != j2) {
            std::cout << prefix << "Value mismatch: " << j1 << " vs " << j2 << std::endl;
        }

    }





    //debug function printing out grids in json file 
    void detray_grids_print( io::detector_grids_payload<std::size_t, io::accel_id>& grids_pd){

        nlohmann::json j;
        const std::string filename = "debug_grids.json";
        //debugging loop for grids
        for (size_t i = 0; i < grids_pd.grids.size(); ++i) {
            for ( auto& grid : grids_pd.grids[i]) {
                nlohmann::json jGrid;
                jGrid["owner_link"] = grid.owner_link.link;
                jGrid["grid_link"] = {{"type", grid.grid_link.type}, {"index", grid.grid_link.index}};

                for (const auto& axis : grid.axes) {
                    nlohmann::json jAxis;
                    jAxis["binning"] = axis.binning;
                    jAxis["bounds"] = axis.bounds;
                    jAxis["label"] = axis.label;
                    jAxis["bins"] = axis.bins;
                    jAxis["edges"] = axis.edges;
                    std::cout << "jAxis: " << jAxis.dump(4) << std::endl;
                    jGrid["axes"].push_back(jAxis);
                }

                for (const auto& bin : grid.bins) {
                    nlohmann::json jBin;
                    jBin["loc_index"] = bin.loc_index;
                    jBin["content"] = bin.content;
                    std::cout << "jBins: " << jBin.dump(4) << std::endl;
                    jGrid["bins"].push_back(jBin);
                }

                j["grids"].push_back(jGrid);
                
            }
        }

        std::ofstream file(filename);
        file << j.dump(4); // 4 spaces for indentation

        for (auto& element : j.items()) {
            if (element.key() == "grid_data") {
                for (auto& array : element.value()) {
                    if (array.size() > 1) {
                        std::cout << "Array with more than one element found: " << array << std::endl;
                    }
                }
            }
        }


        //diff 2 jsons
        /*std::ifstream goldenFile("odd-detray_surface_grids_detray.json");
        nlohmann::json gold_grids;
        goldenFile >> gold_grids;
        print_json_diff(j["grids"], gold_grids["data"]["grids"]["grid_data"]);
        std::cout<<"check debug_grids.json and odd-detray_surface_grids_detray.json"<<std::endl;*/

        return;
    }
    /// detray geometry writer function, debug purposes 
    void detray_detector_print(const detector_t& det){

        std::ofstream outputFile("data_try.json");
        nlohmann::ordered_json out_json;
        typename detector_t::name_map names{};
        out_json["data"] = detray::io::geometry_writer::convert(det, names);
        outputFile << out_json << std::endl;
        return;
    }
    
    /// @return the transform_payload(translation, rotation) of each surface/volume
    /// @param Transform3 acts object, Transform3JsonConverter::Options
    /// @brief convert the acts transform to detray transform payload
    static io::transform_payload detray_converter_transf(
        const Transform3& t, const Transform3JsonConverter::Options& options){
        //nlohmann::json Acts::Transform3JsonConverter::toJson(const Transform3& t, const Transform3JsonConverter::Options& options)

        io::transform_payload p_acts;

        auto translation = t.translation();
        if (translation != Acts::Vector3(0., 0., 0) || options.writeIdentity) {
            std::array<Acts::ActsScalar, 3> tdata = {translation.x(), translation.y(),
                                                    translation.z()};
            p_acts.tr = tdata;
        } else {
            p_acts.tr = {};
        }
        auto rotation =  options.transpose ? t.rotation().transpose() : t.rotation();
        std::array<Acts::ActsScalar, 9> rdata;
        if (rotation != Acts::RotationMatrix3::Identity() || options.writeIdentity) {
            rdata = {
                rotation(0, 0), rotation(0, 1), rotation(0, 2),
                rotation(1, 0), rotation(1, 1), rotation(1, 2),
                rotation(2, 0), rotation(2, 1), rotation(2, 2)};
            p_acts.rot = rdata;
        } else {
            p_acts.rot = {};
        }

        return p_acts;
    }

    /// @return the mask_payload of the surface @param bounds, @param portal
    /// @brief convert the acts bounds to detray mask payload
    static io::mask_payload detray_converter_mask(
        const Acts::SurfaceBounds& bounds, bool portal){
        //Acts::SurfaceBoundsJsonConverter::toJsonDetray

        detray::io::mask_payload mask_pd;
        auto [shape, boundaries] = DetrayJsonHelper::maskFromBounds(bounds, portal);
        mask_pd.shape = static_cast<io::mask_payload::mask_shape>(shape);
        mask_pd.boundaries = static_cast<std::vector<real_io>>(boundaries); //conversion sos??

        //acts/_deps/detray-src/io/include/detray/io/common/geometry_reader.hpp
        //TO DO: use inline single_link_payload convert(const std::size_t idx) 
        detray::io::single_link_payload lnk;
        mask_pd.volume_link = lnk;

        return mask_pd;
    }

    /// @return the surface_payload for portals and volumes by @param Surface acts object
    /// @brief convert the acts surface to detray surface payload and populate the payload
    static io::surface_payload detray_converter_surf(
        const Surface& surface, const Acts::GeometryContext& gctx, const SurfaceJsonConverter::Options& options){
        ///acts/_deps/detray-src/core/include/detray/geometry/detail/surface_descriptor.hpp
        using material_link_payload = io::typed_link_payload<io::material_id>;

        detray::io::surface_payload surf_pd;
        Transform3JsonConverter::Options writtenOption;
    
        surf_pd.transform = detray_converter_transf(surface.transform(gctx), writtenOption);
        surf_pd.source = surface.geometryId().value();
        surf_pd.barcode = std::nullopt;//(long unsigned int)0;
        surf_pd.type = static_cast<detray::surface_id>(options.portal ? 0 : (surface.geometryId().sensitive() > 0 ? 1u : 2u));
        surf_pd.mask = detray_converter_mask(surface.bounds(),options.portal);
        //surf_pd.material =io::typed_link_payload<io::material_id>();
        return surf_pd;
        
    }

    /// construct and @return vector of portals and volumes 
    /// @param gctx, portal, ip, volume, orientedSurfaces, detectorVolumes, option
    /// @brief convert the acts portal to detray surface payload and populate the payload
    static std::vector<io::surface_payload> detray_portals(
        const GeometryContext& gctx, const Experimental::Portal& portal,
        std::size_t ip, const Experimental::DetectorVolume& volume,
        const std::vector<Acts::OrientedSurface>& orientedSurfaces,
        const std::vector<const Experimental::DetectorVolume*>& detectorVolumes,
        const Acts::PortalJsonConverter::Options& option){
        //acts/Plugins/Json/src/PortalJsonConverter.cpp

        std::vector<io::surface_payload> portals {};

        const RegularSurface& surface = portal.surface();
        const auto& volumeLinks = portal.detectorVolumeUpdaters();

        // First assumption for outside link (along direction)
        std::size_t outside = 1u;

        // Find out if you need to take the outside or inside volume
        // for planar surfaces that's easy
        if (surface.type() != Acts::Surface::SurfaceType::Cylinder) {
            // Get the two volume center
            const auto volumeCenter = volume.transform(gctx).translation();
            const auto surfaceCenter = surface.center(gctx);
            const auto surfaceNormal = surface.normal(gctx, surfaceCenter);
            // Get the direction from the volume to the surface, correct link
            const auto volumeToSurface = surfaceCenter - volumeCenter;
            if (volumeToSurface.dot(surfaceNormal) < 0.) {
                outside = 0u;
            }
        } else {
            // This is a cylinder portal, inner cover reverses the normal
            if (ip == 3u) {
            outside = 0u;
            }
        }

        const auto& outsideLink = volumeLinks[outside];
        // Grab the corresponding volume link
        // If it is a single link, we are done
        const auto* instance = outsideLink.instance();
        // Single link cast
        auto singleLink =
            dynamic_cast<const Acts::Experimental::SingleDetectorVolumeImpl*>(
                instance);

        auto [surfaceAdjusted, insidePointer] = orientedSurfaces[ip];

        SurfaceJsonConverter::Options surfaceOptions = option.surfaceOptions;
        surfaceOptions.portal = true;
        // Single link detected - just write it out, we use the oriented surface
        // in order to make sure the size is adjusted
        if (singleLink != nullptr) {
            // Single link can be written out
            std::size_t vLink = findVolume(singleLink->dVolume, detectorVolumes);
            auto portal_pd = detray_converter_surf(*surfaceAdjusted, gctx, surfaceOptions);
            portal_pd.mask.volume_link.link= vLink;
            portals.push_back(portal_pd);
        } else {
            // Multi link detected - 1D
            auto multiLink1D =
                dynamic_cast<const Experimental::BoundVolumesGrid1Impl*>(instance);
            if (multiLink1D != nullptr) {
            // Resolve the multi link 1D
            auto boundaries =
                multiLink1D->indexedUpdater.grid.axes()[0u]->getBinEdges();
            const auto& cast = multiLink1D->indexedUpdater.casts[0u];
            const auto& transform = multiLink1D->indexedUpdater.transform;
            const auto& volumes = multiLink1D->indexedUpdater.extractor.dVolumes;

            // Apply the correction from local to global boundaries
            ActsScalar gCorr = VectorHelpers::cast(transform.translation(), cast);
            std::for_each(boundaries.begin(), boundaries.end(),
                            [&gCorr](ActsScalar& b) { b -= gCorr; });

            // Get the volume indices
            auto surfaceType = surfaceAdjusted->type();
            std::vector<unsigned int> vIndices = {};
            for (const auto& v : volumes) {
                vIndices.push_back(findVolume(v, detectorVolumes));
            }

            // Pick the surface dimension - via poly
            std::array<ActsScalar, 2u> clipRange = {0., 0.};
            std::vector<ActsScalar> boundValues = surfaceAdjusted->bounds().values();
            if (surfaceType == Surface::SurfaceType::Cylinder && cast == binZ) {
                ActsScalar zPosition = surfaceAdjusted->center(gctx).z();
                clipRange = {
                    zPosition - boundValues[CylinderBounds::BoundValues::eHalfLengthZ],
                    zPosition + boundValues[CylinderBounds::BoundValues::eHalfLengthZ]};
            } else if (surfaceType == Surface::SurfaceType::Disc && cast == binR) {
                clipRange = {boundValues[RadialBounds::BoundValues::eMinR],
                            boundValues[RadialBounds::BoundValues::eMaxR]};
            } else {
                throw std::runtime_error(
                    "PortalJsonConverter: surface type not (yet) supported for detray "
                    "conversion, only cylinder and disc are currently supported.");
            }

            // Need to clip the parameter space to the surface dimension
            std::vector<unsigned int> clippedIndices = {};
            std::vector<ActsScalar> clippedBoundaries = {};
            clippedBoundaries.push_back(clipRange[0u]);
            for (const auto [ib, b] : enumerate(boundaries)) {
                if (ib > 0) {
                unsigned int vI = vIndices[ib - 1u];
                ActsScalar highEdge = boundaries[ib];
                if (boundaries[ib - 1] >= clipRange[1u]) {
                    break;
                }
                if (highEdge <= clipRange[0u] ||
                    std::abs(highEdge - clipRange[0u]) < 1e-5) {
                    continue;
                }
                if (highEdge > clipRange[1u]) {
                    highEdge = clipRange[1u];
                }
                clippedIndices.push_back(vI);
                clippedBoundaries.push_back(highEdge);
                }
            }
            // Interpret the clipped information
            //
            // Clipped cylinder case
            if (surfaceType == Surface::SurfaceType::Cylinder) {
                for (auto [ib, b] : enumerate(clippedBoundaries)) {
                if (ib > 0) {
                    // Create sub surfaces
                    std::array<ActsScalar, CylinderBounds::BoundValues::eSize>
                        subBoundValues = {};
                    for (auto [ibv, bv] : enumerate(boundValues)) {
                    subBoundValues[ibv] = bv;
                    }
                    subBoundValues[CylinderBounds::BoundValues::eHalfLengthZ] =
                        0.5 * (b - clippedBoundaries[ib - 1u]);
                    auto subBounds = std::make_shared<CylinderBounds>(subBoundValues);
                    auto subTransform = Transform3::Identity();
                    subTransform.pretranslate(Vector3(
                        0., 0.,
                        clippedBoundaries[ib - 1u] +
                            subBoundValues[CylinderBounds::BoundValues::eHalfLengthZ]));
                    auto subSurface = Surface::makeShared<CylinderSurface>(subTransform, subBounds);
                    
                    auto portal_pd = detray_converter_surf(*subSurface, gctx, surfaceOptions);
                    portal_pd.mask.volume_link.link= clippedIndices[ib - 1u];
                    portals.push_back(portal_pd);
                }
                }
            } else {
                for (auto [ib, b] : enumerate(clippedBoundaries)) {
                    if (ib > 0) {
                        // Create sub surfaces
                        std::array<ActsScalar, RadialBounds::BoundValues::eSize>
                            subBoundValues = {};
                        for (auto [ibv, bv] : enumerate(boundValues)) {
                        subBoundValues[ibv] = bv;
                        }
                        subBoundValues[RadialBounds::BoundValues::eMinR] =
                            clippedBoundaries[ib - 1u];
                        subBoundValues[RadialBounds::BoundValues::eMaxR] = b;
                        auto subBounds = std::make_shared<RadialBounds>(subBoundValues);
                        auto subSurface = Surface::makeShared<DiscSurface>(
                            portal.surface().transform(gctx), subBounds);

                        auto portal_pd = detray_converter_surf(*subSurface, gctx, surfaceOptions);
                        portal_pd.mask.volume_link.link= clippedIndices[ib - 1u];
                        portals.push_back(portal_pd);
                    }
                }
            }

            } else {
            // End of world
            // Write surface with invalid link

                auto portal_pd = detray_converter_surf(*surfaceAdjusted, gctx, surfaceOptions);
                portal_pd.mask.volume_link.link= std::numeric_limits<std::uint_least16_t>::max();
                portals.push_back(portal_pd);
            }
        }
        

        return portals;
    }
    
    /// @return the volume_payload for portals and volumes from acts @param volume, detectorVolumes, gctx
    /// @brief convert the acts volume to detray volume payload and populate the payload
    static io::volume_payload detray_converter_vol(
        const Acts::Experimental::DetectorVolume& volume, 
        const std::vector<const Experimental::DetectorVolume*>& detectorVolumes, 
        const Acts::GeometryContext& gctx){
        //see DetectorVolumeJsonConverter.cpp >>DetectorVolumeJsonConverter::toJsonDetray

        detray::io::volume_payload vol_pd;
        Transform3JsonConverter::Options writtenOption;
        vol_pd.name = volume.name();
        vol_pd.index.link = findVolume(&volume, detectorVolumes);
        std::cout<<vol_pd.name<<std::endl;
        vol_pd.transform = detray_converter_transf(volume.transform(gctx), writtenOption);

        SurfaceJsonConverter::Options surfaceOptions = SurfaceJsonConverter::Options{};

        std::size_t sIndex =0;
        for (const auto surface : volume.surfaces()) {
            io::surface_payload surf_pd = detray_converter_surf(*surface, gctx, surfaceOptions);// acts transf
            surf_pd.index_in_coll= sIndex++;
            surf_pd.mask.volume_link.link= vol_pd.index.link;//link surface' mask to volume
            vol_pd.surfaces.push_back(surf_pd);
        }

        auto orientedSurfaces = volume.volumeBounds().orientedSurfaces(volume.transform(gctx));

        const Acts::PortalJsonConverter::Options options = Acts::PortalJsonConverter::Options{};

        int portals_counter=0;
        for (const auto& [ip, p] : enumerate(volume.portals())) {

            auto portals = (detray_portals(gctx, *p, ip, volume, orientedSurfaces, detectorVolumes, options));
            std::for_each(portals.begin(), portals.end(),
                        [&](auto& portal_pd) {
                            //io::surface_payload portal_pd = detray_converter_portal(*p, gctx);
                            portal_pd.index_in_coll = sIndex++;
                            vol_pd.surfaces.push_back(portal_pd);
                            portals_counter++;
                        });
        }

        return vol_pd;
    }
    
    /// GRID related functions -- in progress
    io::axis_payload axis_converter(const IAxis& ia) {
        ///home/exochell/docker_dir/ACTS_ODD_D/acts/Plugins/Json/src/GridJsonConverter.cpp: nlohmann::json Acts::AxisJsonConverter::toJsonDetray
        io::axis_payload axis_pd;
        axis_pd.bounds =  
            ia.getBoundaryType() == Acts::detail::AxisBoundaryType::Bound ? axis::bounds::e_closed : axis::bounds::e_circular;
        axis_pd.binning = ia.isEquidistant() ? axis::binning::e_regular : axis::binning::e_irregular;
        axis_pd.bins = ia.getNBins();
        if (ia.isEquidistant()) {
            axis_pd.edges = {ia.getBinEdges().front(), ia.getBinEdges().back()};
        } else {
            axis_pd.edges = ia.getBinEdges();
        }

        nlohmann::json output;
        output["axis_pd"]["bounds"] = axis_pd.bounds;
        output["axis_pd"]["binning"] = axis_pd.binning;
        output["axis_pd"]["bins"] = axis_pd.bins;
        output["axis_pd"]["edges"] = axis_pd.edges;
        std::ofstream file("debug_axes.json", std::ios::app);
        file << output.dump(4);

        return axis_pd;
    }

    template <typename grid_type>
    io::grid_payload<std::size_t, io::accel_id> grid_converter(
        const grid_type& grid, bool swapAxis = false) {
        //nlohmann::json toJsonDetray
        //nlohmann::json jGrid;
        // Get the grid axes & potentially swap them
        io::grid_payload<std::size_t, io::accel_id> grid_pd;

        std::array<const Acts::IAxis*, grid_type::DIM> axes = grid.axes();
        if (swapAxis && grid_type::DIM == 2u) {
            std::cout<<"swap axes"<<std::endl;
            std::swap(axes[0u], axes[1u]);
        }
        

        // Fill the axes in the order they are
        for (unsigned int ia = 0u; ia < grid_type::DIM; ++ia) {            
            io::axis_payload axis_pd = axis_converter(*axes[ia]);
            axis_pd.label = static_cast<axis::label>(ia);
            grid_pd.axes.push_back(axis_pd);//push axis to axes
        }

        // 1D connections
        if constexpr (grid_type::DIM == 1u) {
            for (unsigned int ib0 = 1u; ib0 <= axes[0u]->getNBins(); ++ib0) {
            // Lookup bin
                typename grid_type::index_t lbin;
                io::grid_bin_payload<std::size_t>grid_bin_pd; 
                
                lbin[0u] = ib0;
                grid_bin_pd.content = grid.atLocalBins(lbin);
                // Corrected bin for detray
                lbin[0u] = ib0 - 1u;
                grid_bin_pd.loc_index = std::vector<unsigned int>(lbin.begin(), lbin.end());
                grid_pd.bins.push_back(grid_bin_pd);
            }
        }

        // 2D connections
        if constexpr (grid_type::DIM == 2u) {
            for (unsigned int ib0 = 1u; ib0 <= axes[0u]->getNBins(); ++ib0) {
                for (unsigned int ib1 = 1u; ib1 <= axes[1u]->getNBins(); ++ib1) {
                    typename grid_type::index_t lbin;
                    // Lookup bin - respect swap (if it happened) for the lookup
                    lbin[0u] = swapAxis ? ib1 : ib0;
                    lbin[1u] = swapAxis ? ib0 : ib1;

                    io::grid_bin_payload<std::size_t>grid_bin_pd; 

                    nlohmann::json jBin;
                    grid_bin_pd.content = grid.atLocalBins(lbin);
                    // Corrected bin for detray
                    lbin[0u] = ib0 - 1u;
                    lbin[1u] = ib1 - 1u;
                    grid_bin_pd.loc_index = std::vector<unsigned int>(lbin.begin(), lbin.end());
                    grid_pd.bins.push_back(grid_bin_pd);

                }
            }
        }

        std::cout<<"\tgrid_converter"<<std::endl;
        return grid_pd;
    }
    
    //nlohmann::json convertImpl(const index_grid& indexGrid, bool detray = false, bool checkSwap = false) {
    template <typename index_grid>
    io::grid_payload<std::size_t, io::accel_id> convertImpl(const index_grid& indexGrid) {
        
        bool swapAxes = true;

        if constexpr (index_grid::grid_type::DIM == 2u) {
            // Check for axis swap (detray version)
            std::cout<<"swap Axes changed to "<<swapAxes<<std::endl;
            swapAxes = (indexGrid.casts[0u] == binZ && indexGrid.casts[1u] == binPhi);
        }

        io::grid_payload<std::size_t, io::accel_id> grid_pd = grid_converter(indexGrid.grid, swapAxes);

        return grid_pd;
    }

    template <typename instance_type>
    std::optional<io::grid_payload<std::size_t, io::accel_id>> convert(
                const Experimental::SurfaceCandidatesUpdater& delegate,
                bool detray, [[maybe_unused]] const instance_type& refInstance) {
        using GridType =
            typename instance_type::template grid_type<std::vector<std::size_t>>;
        // Defining a Delegate type
        using DelegateType = Experimental::IndexedSurfacesAllPortalsImpl<
            GridType, Experimental::IndexedSurfacesImpl>;
        using SubDelegateType = Experimental::IndexedSurfacesImpl<GridType>;
        
        std::cout<<"convert"<<std::endl;
        // Get the instance
        const auto* instance = delegate.instance();
        auto castedDelegate = dynamic_cast<const DelegateType*>(instance);
        
        if (castedDelegate != nullptr) {
            // Get the surface updator
            io::grid_payload<std::size_t, io::accel_id> grid_pd;
            auto indexedSurfaces = std::get<SubDelegateType>(castedDelegate->updators);
            grid_pd = convertImpl<SubDelegateType>(indexedSurfaces);
            grid_pd.grid_link.type = static_cast<io::accel_id>(DetrayJsonHelper::accelerationLink(indexedSurfaces.casts));
            grid_pd.grid_link.index = std::numeric_limits<std::size_t>::max();
            return grid_pd;
        }

        return std::nullopt;

    }

    template <typename... Args>
    std::vector<io::grid_payload<std::size_t, io::accel_id>> unrollConvert(const Experimental::SurfaceCandidatesUpdater& delegate,
                    bool detray, TypeList<Args...> ) {

        std::cout<<"call convert"<<std::endl;
        std::vector<io::grid_payload<std::size_t, io::accel_id>> grid_pds;

        ((void)(([&]() {
        auto grid_pd = convert(delegate, detray, Args{});
        if (grid_pd.has_value()) {
                grid_pds.push_back(*grid_pd);
            }
        })(), ...));

        return grid_pds;
    }

    static io::detector_grids_payload<std::size_t, io::accel_id> detray_converter_grid(
        const Acts::Experimental::Detector& detector){
    
        io::detector_grids_payload<std::size_t, io::accel_id> grids_pd = io::detector_grids_payload<std::size_t, io::accel_id>();
        auto volumes = detector.volumes();

        for (const auto [iv, volume] : enumerate(volumes)) {

            //Call an equivalent of IndexedSurfacesJsonConverter::toJson
                //check if it is null
            bool detray = true;
            std::cout<<"call unroll"<<std::endl;
            std::vector<io::grid_payload<std::size_t, io::accel_id>> grid_pd = 
                unrollConvert(volume->surfaceCandidatesUpdater(), detray, GridAxisGenerators::PossibleAxes{});
            
            for (auto& grid : grid_pd) {
                detray::io::single_link_payload lnk;
                lnk.link = iv;
                grid.owner_link = lnk;
                grids_pd.grids[iv].push_back(grid);
            }

            //TO DO:: volume link

            //TO DO:: header payload
            
        }

        detray_grids_print(grids_pd);

        return grids_pd;
    }


    /// @return the geo_header_payload from @param detector object of ACTS
    static io::geo_header_payload detray_converter_head(
        const Acts::Experimental::Detector& detector){
        
        detray::io::geo_header_payload header_pd;
        detray::io::common_header_payload header_data_pd;
        

        //TO DO use inline common_header_payload convert(const std::string_view det_name, const std::string_view tag)
        header_data_pd.version = io::detail::get_detray_version();
        header_data_pd.detector = detector.name();
        header_data_pd.tag = "geometry"; 
        header_data_pd.date = io::detail::get_current_date();

        return header_pd;
    }

    /// @brief visit all ACTS detector information, depth-first hierarchically, populate the corresponding payloads and convert to detray detector 
    /// @return detray detector from @param detector and @param gctx of ACTS 
    detector_t detray_tree_converter(
        const Acts::Experimental::Detector& detector,
        const Acts::GeometryContext& gctx, vecmem::memory_resource& mr){
        
        detray::io::detector_payload dp;
        std::cout<<"-----tree converter-------"<<std::endl;
        for (const auto volume : detector.volumes()) {
            dp.volumes.push_back(detray_converter_vol(*volume, detector.volumes(), gctx));
            std::cout<<std::endl;
        }
        typename detector_t::name_map names{};
        detector_builder<default_metadata> det_builder{};

        detray::io::geometry_reader::convert<detector_t>(det_builder, names, dp);
        
        //io/include/detray/io/common/surface_grid_reader.hpp
        detray::io::surface_grid_reader<
                                    typename detector_t::surface_type,
                                    std::integral_constant<std::size_t, 0>,
                                    std::integral_constant<std::size_t, 2>>
                                    ::convert<detector_t>(det_builder, names, detray_converter_grid(detector));
        

        std::cout<<"detector_t"<<std::endl;
        detray::types::print<types::list<detector_t>>();

        detector_t detrayDet(det_builder.build(mr));
        detray::detail::check_consistency(detrayDet);

        return std::move(detrayDet);
    }

}

