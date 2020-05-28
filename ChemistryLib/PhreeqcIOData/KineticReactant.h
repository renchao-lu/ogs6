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

#include <boost/optional/optional.hpp>
#include <iosfwd>
#include <string>
#include <vector>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "MeshLib/PropertyVector.h"
#include "Output.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct KineticReactant
{
    KineticReactant(std::string&& name_,
                    std::string&& chemical_formula_,
                    double molar_volume_,
                    double initial_volume_fraction_,
                    MeshLib::PropertyVector<double>* volume_fraction_,
                    MeshLib::PropertyVector<double>* volume_fraction_change_,
                    std::vector<double>&& parameters_,
                    bool const fix_mass_)
        : name(std::move(name_)),
          chemical_formula(std::move(chemical_formula_)),
          molar_volume(molar_volume_),
          initial_volume_fraction(initial_volume_fraction_),
          volume_fraction(volume_fraction_),
          volume_fraction_change(volume_fraction_change_),
          parameters(std::move(parameters_)),
          fix_mass(fix_mass_)
    {
    }

    void print(std::ostream& os,
               std::size_t const chemical_system_id,
               double const porosity,
               double const density) const;

    std::string const name;
    std::string const chemical_formula;
    double molar_volume;
    double initial_volume_fraction;
    MeshLib::PropertyVector<double>* volume_fraction;
    MeshLib::PropertyVector<double>* volume_fraction_change;
    std::vector<double> const parameters;
    bool const fix_mass;
    static const ItemType item_type = ItemType::KineticReactant;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
