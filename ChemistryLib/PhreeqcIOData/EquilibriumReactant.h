/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "MeshLib/PropertyVector.h"
#include "Output.h"

namespace BaseLib
{
class ConfigTree;
}

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct EquilibriumReactant
{
    EquilibriumReactant(std::string&& name_,
                        double initial_mass_,
                        MeshLib::PropertyVector<double>* mass_,
                        MeshLib::PropertyVector<double>* mass_loss_,
                        double saturation_index_)
        : name(std::move(name_)),
          initial_mass(initial_mass_),
          mass(mass_),
          mass_loss(mass_loss_),
          saturation_index(saturation_index_)
    {
    }

    void print(std::ostream& os,
               std::size_t const chemical_system_id,
               double const porosity,
               double const rho_w) const;

    std::string const name;
    double initial_mass;
    MeshLib::PropertyVector<double>* mass;
    MeshLib::PropertyVector<double>* mass_loss;
    double const saturation_index;
    static const ItemType item_type = ItemType::EquilibriumReactant;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
