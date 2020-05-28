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
    MeshLib::Mesh const& mesh)
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

        double const molar_volume =
            //! \ogs_file_param{prj__chemical_system__kinetic_reactants__kinetic_reactant__molar_volume}
            reactant_config.getConfigParameter<double>("molar_volume");

        double const initial_volume_fraction =
            //! \ogs_file_param{prj__chemical_system__kinetic_reactants__kinetic_reactant__initial_volume_fraction}
            reactant_config.getConfigParameter<double>(
                "initial_volume_fraction");

        auto volume_fraction = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh),
            name,
            MeshLib::MeshItemType::IntegrationPoint,
            1);

        auto volume_fraction_change = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh),
            "volume_fraction_change_" + name,
            MeshLib::MeshItemType::IntegrationPoint,
            1);

        auto parameters =
            //! \ogs_file_param{prj__chemical_system__kinetic_reactants__kinetic_reactant__parameters}
            reactant_config.getConfigParameter<std::vector<double>>(
                "parameters", {});

        bool const fix_mass =
            //! \ogs_file_param{prj__chemical_system__kinetic_reactants__kinetic_reactant__fix_amount}
            reactant_config.getConfigParameter<bool>("fix_amount", false);

        if (chemical_formula.empty() && fix_mass)
        {
            OGS_FATAL(
                "fix_amount can only be used if a chemical_formula has been "
                "defined");
        }

        kinetic_reactants.emplace_back(std::move(name),
                                       std::move(chemical_formula),
                                       molar_volume,
                                       initial_volume_fraction,
                                       volume_fraction,
                                       volume_fraction_change,
                                       std::move(parameters),
                                       fix_mass);
    }

    return kinetic_reactants;
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
