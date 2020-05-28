/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateChemicalSolverInterface.h"

#include <iphreeqc/src/src/phreeqcpp/cxxKinetics.h>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/FileTools.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MeshLib/Mesh.h"

#include "Common/CreateReactionRate.h"

#include "PhreeqcIO.h"
#include "PhreeqcIOData/ChemicalSystem.h"
#include "PhreeqcIOData/CreateChemicalSystem.h"
#include "PhreeqcIOData/CreateKnobs.h"
#include "PhreeqcIOData/CreateOutput.h"
#include "PhreeqcIOData/CreateSurface.h"
#include "PhreeqcIOData/CreateUserPunch.h"
#include "PhreeqcIOData/Dump.h"
#include "PhreeqcIOData/EquilibriumReactant.h"
#include "PhreeqcIOData/KineticReactant.h"
#include "PhreeqcIOData/Knobs.h"
#include "PhreeqcIOData/ReactionRate.h"
#include "PhreeqcIOData/Surface.h"
#include "PhreeqcIOData/UserPunch.h"

#include "PhreeqcKernel.h"
#include "PhreeqcKernelData/AqueousSolution.h"
#include "PhreeqcKernelData/CreateAqueousSolution.h"
#include "PhreeqcKernelData/CreateEquilibriumReactants.h"
#include "PhreeqcKernelData/CreateKineticReactant.h"
#include "PhreeqcKernelData/ReactionRate.h"

namespace
{
std::string parseDatabasePath(BaseLib::ConfigTree const& config)
{
    // database
    //! \ogs_file_param{prj__chemical_system__database}
    auto const database = config.getConfigParameter<std::string>("database");
    auto path_to_database =
        BaseLib::joinPaths(BaseLib::getProjectDirectory(), database);

    if (!BaseLib::IsFileExisting(path_to_database))
    {
        OGS_FATAL("Not found the specified thermodynamicdatabase: {:s}",
                  path_to_database);
    }

    INFO("Found the specified thermodynamic database: {:s}", path_to_database);

    return path_to_database;
}
}  // namespace

namespace ChemistryLib
{
template <>
std::shared_ptr<ChemicalSolverInterface>
createChemicalSolverInterface<ChemicalSolver::Phreeqc>(
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
    BaseLib::ConfigTree const& config, std::string const& output_directory,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    auto mesh_name =
        //! \ogs_file_param{prj__chemical_system__mesh}
        config.getConfigParameter<std::string>("mesh");

    // Find and extract mesh from the list of meshes.
    auto const& mesh = *BaseLib::findElementOrError(
        std::begin(meshes), std::end(meshes),
        [&mesh_name](auto const& mesh) {
            assert(mesh != nullptr);
            return mesh->getName() == mesh_name;
        },
        "Required mesh with name '" + mesh_name + "' not found.");

    assert(mesh.getID() != 0);
    DBUG("Found mesh '{:s}' with id {:d}.", mesh.getName(), mesh.getID());

    auto path_to_database = parseDatabasePath(config);

    // chemical system
    auto chemical_system =
        PhreeqcIOData::createChemicalSystem(config, *meshes[0]);

    // rates
    auto reaction_rates = createReactionRates<PhreeqcIOData::ReactionRate>(
        //! \ogs_file_param{prj__chemical_system__rates}
        config.getConfigSubtreeOptional("rates"));

    // surface
    auto surface = PhreeqcIOData::createSurface(
        //! \ogs_file_param{prj__chemical_system__surface}
        config.getConfigSubtreeOptional("surface"));

    // dump
    auto const project_file_name = BaseLib::joinPaths(
        output_directory,
        BaseLib::extractBaseNameWithoutExtension(config.getProjectFileName()));

    auto dump = surface.empty()
                    ? nullptr
                    : std::make_unique<PhreeqcIOData::Dump>(project_file_name);

    // knobs
    auto knobs = PhreeqcIOData::createKnobs(
        //! \ogs_file_param{prj__chemical_system__knobs}
        config.getConfigSubtree("knobs"));

    // user punch
    auto user_punch = PhreeqcIOData::createUserPunch(
        //! \ogs_file_param{prj__chemical_system__user_punch}
        config.getConfigSubtreeOptional("user_punch"), *meshes[0]);

    // output
    auto const use_high_precision =
        //! \ogs_file_param{prj__chemical_system__use_high_precision}
        config.getConfigParameter<bool>("use_high_precision", true);
    auto output = PhreeqcIOData::createOutput(
        chemical_system->aqueous_solution->components,
        chemical_system->equilibrium_reactants,
        chemical_system->kinetic_reactants, user_punch, use_high_precision,
        project_file_name);

    auto porosity = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(*meshes[0]),
        "porosity",
        MeshLib::MeshItemType::IntegrationPoint,
        1);

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    return std::make_shared<PhreeqcIOData::PhreeqcIO>(
        std::move(project_file_name), *meshes[mesh.getID()],
        std::move(path_to_database), std::move(chemical_system),
        std::move(reaction_rates), std::move(surface), std::move(user_punch),
        std::move(output), std::move(dump), std::move(knobs), porosity,
        std::move(media_map));
}

template <>
std::shared_ptr<ChemicalSolverInterface>
createChemicalSolverInterface<ChemicalSolver::PhreeqcKernel>(
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
    //    std::vector<std::pair<int, std::string>> const&
    //        process_id_to_component_name_map,
    BaseLib::ConfigTree const& config, std::string const& /*output_directory*/,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    auto mesh = *meshes[0];
    auto path_to_database = parseDatabasePath(config);

    // solution
    auto aqueous_solution = PhreeqcKernelData::createAqueousSolution(
        //! \ogs_file_param{prj__chemical_system__solution}
        config.getConfigSubtree("solution")/*,
        process_id_to_component_name_map*/);

    // kinetic reactants
    auto kinetic_reactants = PhreeqcKernelData::createKineticReactants(
        //! \ogs_file_param{prj__chemical_system__kinetic_reactants}
        config.getConfigSubtreeOptional("kinetic_reactants"), mesh);

    // rates
    auto reaction_rates = createReactionRates<PhreeqcKernelData::ReactionRate>(
        //! \ogs_file_param{prj__chemical_system__rates}
        config.getConfigSubtreeOptional("rates"));

    // equilibrium reactants
    auto equilibrium_reactants = PhreeqcKernelData::createEquilibriumReactants(
        //! \ogs_file_param{prj__chemical_system__equilibrium_reactants}
        config.getConfigSubtreeOptional("equilibrium_reactants"), mesh);

    return std::make_shared<PhreeqcKernelData::PhreeqcKernel>(
        mesh.getNumberOfBaseNodes(), /* process_id_to_component_name_map,*/
        std::move(path_to_database), std::move(aqueous_solution),
        std::move(equilibrium_reactants), std::move(kinetic_reactants),
        std::move(reaction_rates));
}
}  // namespace ChemistryLib
