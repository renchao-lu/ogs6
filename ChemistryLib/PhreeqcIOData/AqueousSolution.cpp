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

#include "AqueousSolution.h"
#include "ChemistryLib/Common/ChargeBalance.h"
#include "MaterialLib/PhysicalConstant.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
void AqueousSolution::print(std::ostream& os,
                            std::size_t const chemical_system_id,
                            double const density) const
{
    os << "temp " << temperature << "\n";

    os << "pressure "
       << (*pressure)[chemical_system_id] /
              MaterialLib::PhysicalConstant::AtmosphericPressure
       << "\n";

    switch (charge_balance)
    {
        case ChargeBalance::pH:
            os << "pH " << -std::log10((*pH)[chemical_system_id]) << " charge"
               << "\n";
            os << "pe " << (*pe)[chemical_system_id] << "\n";
            break;
        case ChargeBalance::pe:
            os << "pH " << -std::log10((*pH)[chemical_system_id]) << "\n";
            os << "pe " << (*pe)[chemical_system_id] << " charge"
               << "\n";
            break;
        case ChargeBalance::Unspecified:
            os << "pH " << -std::log10((*pH)[chemical_system_id]) << "\n";
            os << "pe " << (*pe)[chemical_system_id] << "\n";
            break;
    }

    // Molality
    os << "units mol/kgw\n";

    for (auto const& component : components)
    {
        // unit conversion from mol/m3 to mol/kgw
        os << component.name << " "
           << (*component.amount)[chemical_system_id] / density << "\n";
    }

    os << "\n";
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
