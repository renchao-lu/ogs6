/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Surface.h"

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
SurfaceCharge::SurfaceCharge()
{
    int i = 1;
}

SurfaceComponent::SurfaceComponent(std::string const& site_name,
                                   double const site_density,
                                   double const specific_surface_area,
                                   double const mass)
{
    formula = site_name;
    moles = site_density;
}

Surface::Surface(std::vector<double> const& surf)
{
    new_def = true;

}
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
