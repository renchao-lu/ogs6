/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <boost/optional/optional.hpp>

#include "BaseLib/ConfigTree.h"
#include "CreateEquilibriumReactants.h"
#include "EquilibriumReactant.h"
#include "MeshLib/Mesh.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::vector<EquilibriumReactant> createEquilibriumReactants(
    boost::optional<BaseLib::ConfigTree> const& config,
    MeshLib::Mesh const& mesh)
{
    if (!config)
    {
        return {};
    }

    std::vector<EquilibriumReactant> equilibrium_reactants;
    for (
        auto const& equilibrium_reactant_config :
        //! \ogs_file_param{prj__chemical_system__equilibrium_reactants__phase_component}
        config->getConfigSubtreeList("phase_component"))
    {
        auto name =
            //! \ogs_file_param{prj__chemical_system__equilibrium_reactants__phase_component__name}
            equilibrium_reactant_config.getConfigParameter<std::string>("name");

        double const initial_mass =
            //! \ogs_file_param{prj__chemical_system__equilibrium_reactants__phase_component__initial_amount}
            equilibrium_reactant_config.getConfigParameter<double>(
                "initial_amount");

        auto mass = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh),
            name,
            MeshLib::MeshItemType::IntegrationPoint,
            1);

        auto mass_loss = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh),
            "mass_loss_" + name,
            MeshLib::MeshItemType::IntegrationPoint,
            1);

        double const saturation_index =
            //! \ogs_file_param{prj__chemical_system__equilibrium_reactants__phase_component__saturation_index}
            equilibrium_reactant_config.getConfigParameter<double>(
                "saturation_index");

        equilibrium_reactants.emplace_back(
            std::move(name), initial_mass, mass, mass_loss, saturation_index);
    }

    return equilibrium_reactants;
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
