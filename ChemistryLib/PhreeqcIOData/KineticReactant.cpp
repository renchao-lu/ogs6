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

#include "KineticReactant.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
void KineticReactant::print(std::ostream& os,
                            std::size_t const chemical_system_id,
                            double const porosity,
                            double const density) const
{
    os << name << "\n";

    if (!chemical_formula.empty())
    {
        os << "-formula " << chemical_formula << "\n";
    }

    // kinetic reactant
    // volume [m3]: volume_fraction * total_volume
    // mass [mol]: volume_fraction * total_volume / molar_volume
    // after the water mass [kg] is scaled from porosity * total_volume * rho to
    // 1, i.e., scaling factor [-]: 1 / (porosity * total_volume * rho) kinetic
    // reactant mass after scaling [mol]: volume_fraction / molar_volume /
    // porosity / rho
    auto const mass = (*volume_fraction)[chemical_system_id] / molar_volume /
                      porosity / density;
    os << "-m  " << mass << "\n";

    if (!parameters.empty())
    {
        os << "-parms";
        for (auto const& parameter : parameters)
        {
            os << " " << parameter;
        }
        os << "\n";
    }
}
const ItemType KineticReactant::item_type;
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
