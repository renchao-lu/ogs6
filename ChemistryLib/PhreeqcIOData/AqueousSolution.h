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

namespace ChemistryLib
{
enum class ChargeBalance;

namespace PhreeqcIOData
{
struct Component
{
    explicit Component(std::string&& name_) : name(std::move(name_)) {}

    std::string const name;
    std::unique_ptr<GlobalVector> amount;
    static const ItemType item_type = ItemType::Component;
};

struct AqueousSolution
{
    AqueousSolution(double temperature_,
                    MeshLib::PropertyVector<double>* pe_,
                    double pe0_,
                    std::vector<Component>&& components_,
                    ChargeBalance charge_balance_,
                    MeshLib::PropertyVector<double>* water_mass_)
        : temperature(temperature_),
          pe(pe_),
          pe0(pe0_),
          components(std::move(components_)),
          charge_balance(charge_balance_),
          water_mass(water_mass_)
    {
    }

    void print(std::ostream& os, std::size_t const chemical_system_id,
               double const density) const;

    double temperature;
    std::unique_ptr<GlobalVector> pressure;
    std::unique_ptr<GlobalVector> pH;
    MeshLib::PropertyVector<double>* pe;
    double pe0;
    std::vector<Component> components;
    ChargeBalance const charge_balance;
    MeshLib::PropertyVector<double>* water_mass;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
