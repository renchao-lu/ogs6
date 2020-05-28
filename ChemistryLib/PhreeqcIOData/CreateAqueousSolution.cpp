/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateAqueousSolution.h"
#include "AqueousSolution.h"
#include "BaseLib/ConfigTree.h"
#include "ChemistryLib/Common/CreateChargeBalance.h"
#include "CreateSolutionComponent.h"
#include "MeshLib/Mesh.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::unique_ptr<AqueousSolution> createAqueousSolution(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& mesh)
{
    //! \ogs_file_param{prj__chemical_system__solution__temperature}
    auto const temperature = config.getConfigParameter<double>("temperature");

    auto pe = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh),
        "pe",
        MeshLib::MeshItemType::IntegrationPoint,
        1);

    //! \ogs_file_param{prj__chemical_system__solution__pe}
    auto const pe0 = config.getConfigParameter<double>("pe");

    auto components = createSolutionComponents(config);

    auto charge_balance = createChargeBalance(config);

    auto water_mass = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh),
        "water_mass",
        MeshLib::MeshItemType::IntegrationPoint,
        1);

    return std::make_unique<AqueousSolution>(temperature,
                                             pe,
                                             pe0,
                                             std::move(components),
                                             charge_balance,
                                             water_mass);
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
