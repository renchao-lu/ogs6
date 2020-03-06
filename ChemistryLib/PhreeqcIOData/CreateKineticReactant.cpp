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
#include "CreateKineticReactant.h"
#include "KineticReactant.h"
#include "MeshLib/Mesh.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::vector<KineticReactant> createKineticReactants(
    boost::optional<BaseLib::ConfigTree> const& config,
    MeshLib::Mesh const& mesh,
    MeshLib::PropertyVector<std::size_t> const& chemical_system_map)
{
    if (!config)
    {
        return {};
    }

    std::vector<KineticReactant> kinetic_reactants;
    for (
        auto const& reactant_config :
        //! \ogs_file_param{prj__chemical_system__kinetic_reactants__kinetic_reactant}
        config->getConfigSubtreeList("kinetic_reactant"))
    {
        //! \ogs_file_param{prj__chemical_system__kinetic_reactants__kinetic_reactant__name}
        auto name = reactant_config.getConfigParameter<std::string>("name");

        auto chemical_formula =
            //! \ogs_file_param{prj__chemical_system__kinetic_reactants__kinetic_reactant__chemical_formula}
            reactant_config.getConfigParameter<std::string>("chemical_formula",
                                                            "");

        double const initial_amount =
            //! \ogs_file_param{prj__chemical_system__kinetic_reactants__kinetic_reactant__initial_amount}
            reactant_config.getConfigParameter<double>("initial_amount");

        auto parameters =
            //! \ogs_file_param{prj__chemical_system__kinetic_reactants__kinetic_reactant__parameters}
            reactant_config.getConfigParameter<std::vector<double>>(
                "parameters", {});

        bool const fix_amount =
            //! \ogs_file_param{prj__chemical_system__kinetic_reactants__kinetic_reactant__fix_amount}
            reactant_config.getConfigParameter<bool>("fix_amount", false);

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

        auto transferred_amount = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh),
            "dm_" + name,
            MeshLib::MeshItemType::Node,
            1);

        std::fill(std::begin(*transferred_amount),
                  std::end(*transferred_amount),
                  std::numeric_limits<double>::quiet_NaN());

        std::for_each(chemical_system_map.begin(),
                      chemical_system_map.end(),
                      [&transferred_amount](auto const& global_id) {
                          (*transferred_amount)[global_id] = 0.;
                      });

        if (chemical_formula.empty() && fix_amount)
        {
            OGS_FATAL(
                "fix_amount can only be used if a chemical_formula has been "
                "defined");
        }

        kinetic_reactants.emplace_back(std::move(name),
                                       std::move(chemical_formula),
                                       amount,
                                       transferred_amount,
                                       std::move(parameters),
                                       fix_amount);
    }

    return kinetic_reactants;
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
