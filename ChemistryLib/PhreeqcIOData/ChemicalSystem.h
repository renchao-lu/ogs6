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

#include <memory>
#include <vector>

#include "AqueousSolution.h"

namespace MeshLib
{
class Mesh;
}

namespace MaterialPropertyLib
{
class MaterialSpatialDistributionMap;
}

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct KineticReactant;
struct EquilibriumReactant;

struct ChemicalSystem
{
    ChemicalSystem(std::unique_ptr<AqueousSolution>&& aqueous_solution_,
                   std::vector<KineticReactant>&& kinetic_reactants_,
                   std::vector<EquilibriumReactant>&& equilibrium_reactants_)
        : aqueous_solution(std::move(aqueous_solution_)),
          kinetic_reactants(std::move(kinetic_reactants_)),
          equilibrium_reactants(std::move(equilibrium_reactants_))
    {
    }

    void initialize(std::size_t const num_chemical_systems);

    void print(
        std::size_t const mesh_item_id, std::size_t const chemical_system_id,
        std::ostream& os, double const dt,
        MaterialPropertyLib::MaterialSpatialDistributionMap const& media_map,
        double const& porosity);

    std::unique_ptr<AqueousSolution> aqueous_solution;
    std::vector<KineticReactant> kinetic_reactants;
    std::vector<EquilibriumReactant> equilibrium_reactants;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
