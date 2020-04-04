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
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "CreateEquilibriumReactants.h"
#include "EquilibriumReactants.h"
#include "MeshLib/Mesh.h"

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
std::unique_ptr<EquilibriumReactants> createEquilibriumReactants(
    boost::optional<BaseLib::ConfigTree> const& config,
    MeshLib::Mesh const& mesh,
    MeshLib::PropertyVector<std::size_t> const& chemical_system_map)
{
    if (!config)
    {
        return nullptr;
    }

    std::vector<PhaseComponent> phase_components;
    for (
        auto const& phase_component_config :
        //! \ogs_file_param{prj__chemical_system__equilibrium_reactants__phase_component}
        config->getConfigSubtreeList("phase_component"))
    {
        auto name =
            //! \ogs_file_param{prj__chemical_system__equilibrium_reactants__phase_component__name}
            phase_component_config.getConfigParameter<std::string>("name");

        double const initial_amount =
            //! \ogs_file_param{prj__chemical_system__equilibrium_reactants__phase_component__initial_amount}
            phase_component_config.getConfigParameter<double>("initial_amount");

        double const saturation_index =
            //! \ogs_file_param{prj__chemical_system__equilibrium_reactants__phase_component__saturation_index}
            phase_component_config.getConfigParameter<double>(
                "saturation_index");

        auto amount = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh),
            name,
            MeshLib::MeshItemType::Node,
            1);

        std::fill(std::begin(*amount),
                  std::end(*amount),
                  std::numeric_limits<double>::quiet_NaN());

        std::for_each(chemical_system_map.begin(),
                      chemical_system_map.end(),
                      [&amount, initial_amount](auto const& global_id) {
                          (*amount)[global_id] = initial_amount;
                      });

        phase_components.emplace_back(
            std::move(name), initial_amount, saturation_index);
    }

    return std::make_unique<EquilibriumReactants>(phase_components);
}
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
