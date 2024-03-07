# This file is part of the Acts project.
#
# Copyright (C) 2021 CERN for the benefit of the Acts project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# CMake include(s).
cmake_minimum_required( VERSION 3.11 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building Detray as part of the Acts project" )

# Declare where to get Detray from.
set( ACTS_DETRAY_GIT_REPOSITORY "https://github.com/acts-project/detray.git"
   CACHE STRING "Git repository to take Detray from" )
set( ACTS_DETRAY_GIT_TAG "v0.55.0" CACHE STRING "Version of Detray to build" )
mark_as_advanced( ACTS_DETRAY_GIT_REPOSITORY ACTS_DETRAY_GIT_TAG )
FetchContent_Declare( Detray
   GIT_REPOSITORY "${ACTS_DETRAY_GIT_REPOSITORY}"
   GIT_TAG "${ACTS_DETRAY_GIT_TAG}" )


set( DETRAY_SETUP_ACTSVG FALSE CACHE BOOL
   "Do not set up Actsvg as part of Detray" )
set( DETRAY_BENCHMARKS FALSE CACHE BOOL
   "Disable benchmark tests" )
set( DETRAY_BUILD_TUTORIALS FALSE CACHE BOOL
   "Don't build the tutorial executables of Detray" )
set( DETRAY_BUILD_TESTING FALSE CACHE BOOL
   "Don't Build the (unit) tests of Detray" )
set( DETRAY_CUSTOM_SCALARTYPE float CACHE STRING
   "Use float type in the Detray code" )
set( DETRAY_EIGEN_PLUGIN FALSE CACHE BOOL
   "Don't build Eigen math plugin" )
set( DETRAY_SETUP_ALGEBRA_PLUGINS ON CACHE BOOL
  "Set up the Algebra Plugins target(s)" )
set( DETRAY_SETUP_BENCHMARK FALSE CACHE BOOL
   "Don't set up the Google Benchmark target(s)" )
set( DETRAY_SETUP_COVFIE ON CACHE BOOL
   "Set up the covfie target" )
set( DETRAY_SETUP_GOOGLETEST OFF CACHE BOOL
   "Set up the covfie target" )
set( DETRAY_VC_PLUGIN OFF CACHE BOOL
  "Set up the covfie target" )
set( DETRAY_SVG_DISPLAY ON CACHE BOOL
   "Set up the covfie target" )
set( DETRAY_SETUP_DFELIBS FALSE CACHE BOOL
   "Do not set up Dfelibs as part of Detray" )
set( DETRAY_SETUP_NLOHMANN FALSE CACHE BOOL
   "Do not set up Nlohmann as part of Detray" )
#set( DETRAY_SETUP_VECMEM FALSE CACHE BOOL
#"Do not set up Vecmem as part of Detray" )

#Now set up its build.
FetchContent_MakeAvailable(Detray)
