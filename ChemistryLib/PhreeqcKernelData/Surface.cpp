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
SurfaceCharge::SurfaceCharge(double const specific_surface_area, double const mass)
{
    specific_area = specific_surface_area;
    grams = mass;
}

SurfaceComponent::SurfaceComponent(std::string const& site_name,
                                   double const site_density)
{
    formula = site_name;
    moles = site_density;
}

Surface::Surface(std::vector<SurfaceComponent> surface_components,
                 std::vector<SurfaceCharge> surface_charges_)
{
    new_def = true;
    sites_units = SITES_DENSITY;

    for (auto const& surface_component: surface_components)
    {
        surface_comps.push_back(*surface_component.castToBaseClass());
    }

    for (auto const& surface_charge : surface_charges_)
    {
        surface_charges.push_back(*surface_charge.castToBaseClass());
    }
}
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
