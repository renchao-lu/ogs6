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

#include <string>
#include <vector>

#include "ThirdParty/iphreeqc/src/src/phreeqcpp/PPassemblage.h"

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
class EquilibriumPhase final : private cxxPPassemblageComp
{
public:
    EquilibriumPhase(std::string&& phase_name, double const initial_amount,
                     double const saturation_index);

    cxxPPassemblageComp const* castToBaseClass() const
    {
        return static_cast<cxxPPassemblageComp const*>(this);
    }

    std::string getName() const { return Get_name(); }
};

class Equilibriums final : private cxxPPassemblage
{
public:
    explicit Equilibriums(std::vector<EquilibriumPhase>& EquilibriumPhases);
    /*
        bool isKineticReactantDefined() const
        {
            return !Get_kinetics_comps().empty();
        }

        void setChemicalSystemID(std::size_t const chemical_system_id)
        {
            Set_n_user_both(chemical_system_id);
        }

        cxxKinetics const* castToBaseClass() const
        {
            return static_cast<cxxKinetics const*>(this);
        }
        */
};

}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
