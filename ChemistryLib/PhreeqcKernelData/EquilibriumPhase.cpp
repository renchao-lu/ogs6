/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <algorithm>

#include "EquilibriumPhase.h"

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
EquilibriumPhase::EquilibriumPhase(std::string&& phase_name,
                                   double const initial_amount,
                                   double const saturation_index)
{
    name = phase_name;
    moles = initial_amount;
    si = saturation_index;
    si_org = saturation_index;
}

Equilibriums::Equilibriums(std::vector<EquilibriumPhase>& equilibrium_phases)
{
    // new_def = true;
    for (auto& equilibrium_phase : equilibrium_phases)
    {
        auto name = equilibrium_phase.getName();
        pp_assemblage_comps[name] = *equilibrium_phase.castToBaseClass();
    }
}
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
