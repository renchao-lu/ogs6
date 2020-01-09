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
    SurfaceCharge(double const specific_surface_area,
                  double const mass);

    cxxSurfaceCharge const* castToBaseClass() const
    {
        return static_cast<cxxSurfaceCharge const*>(this);
    }
};

class SurfaceComponent final : private cxxSurfaceComp
{
public:
    SurfaceComponent(std::string const& site_name,
                     double const site_density);

    cxxSurfaceComp const* castToBaseClass() const
    {
        return static_cast<cxxSurfaceComp const*>(this);
    }
};

class Surface final : private cxxSurface
{
public:
    Surface(std::vector<SurfaceComponent> surface_components,
            std::vector<SurfaceCharge> surface_charges_);

    void setChemicalSystemID(std::size_t const chemical_system_id)
    {
        Set_n_user_both(chemical_system_id);
    }

    void equilibrateWithSolution(bool const flag,
                                 std::size_t const chemical_system_id)
    {
        solution_equilibria = true;
        n_solution = chemical_system_id;
    }

    cxxSurface* castToBaseClass()
    {
        return static_cast<cxxSurface*>(this);
    }

    void sort()
    {
        Sort_comps();
    }
};
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
