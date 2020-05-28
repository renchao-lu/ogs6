/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <ostream>

#include "EquilibriumReactant.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
void EquilibriumReactant::print(std::ostream& os,
                                std::size_t const chemical_system_id,
                                double const porosity,
                                double const rho_w) const
{
    // equilibrium reactant
    // volume [m3]: volume_fraction * total_volume
    // mass [mol]: volume_fraction * total_volume / molar_volume
    // after the water mass [kg] is scaled from porosity * total_volume * rho to
    // 1 i.e., scaling factor [-]: 1 / (porosity * total_volume * rho)
    // equilibrium reactant mass after scaling [mol]:
    // volume_fraction / molar_volume / porosity / rho
    // TODO
    os << name << " " << saturation_index << " " << saturation_index << " "
       << (*mass)[chemical_system_id] << "\n";
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
