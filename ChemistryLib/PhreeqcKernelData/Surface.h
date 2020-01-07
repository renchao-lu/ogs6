/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "ThirdParty/iphreeqc/src/src/phreeqcpp/Surface.h"

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
class SurfaceCharge final : private cxxSurfaceCharge
{
public:
    SurfaceCharge();
};

class SurfaceComponent final : private cxxSurfaceComp
{
public:
    SurfaceComponent(std::string const& site_name,
                     double const site_density,
                     double const specific_surface_area,
                     double const mass);
};

class Surface final : private cxxSurface
{
public:
    Surface(std::vector<double> const& surf);
};
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
