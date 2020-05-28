/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PhreeqcIO.h"

#include <iphreeqc/src/src/IPhreeqc.h>

#include <boost/algorithm/string.hpp>
#include <cmath>
#include <fstream>

#include "BaseLib/Algorithm.h"
#include "BaseLib/ConfigTreeUtil.h"
#include "MaterialLib/MPL/Medium.h"
#include "MeshLib/Mesh.h"
#include "ParameterLib/SpatialPosition.h"

#include "PhreeqcIOData/AqueousSolution.h"
#include "PhreeqcIOData/ChemicalSystem.h"
#include "PhreeqcIOData/Dump.h"
#include "PhreeqcIOData/EquilibriumReactant.h"
#include "PhreeqcIOData/KineticReactant.h"
#include "PhreeqcIOData/Knobs.h"
#include "PhreeqcIOData/Output.h"
#include "PhreeqcIOData/ReactionRate.h"
#include "PhreeqcIOData/Surface.h"
#include "PhreeqcIOData/UserPunch.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
namespace
{
template <typename DataBlock>
std::ostream& operator<<(std::ostream& os,
                         std::vector<DataBlock> const& data_blocks)
{
    std::copy(data_blocks.begin(), data_blocks.end(),
              std::ostream_iterator<DataBlock>(os));
    return os;
}
}  // namespace

PhreeqcIO::PhreeqcIO(
    std::string const project_file_name,
    MeshLib::Mesh const& mesh,
    std::string&& database,
    std::unique_ptr<ChemicalSystem>&& chemical_system,
    std::vector<ReactionRate>&& reaction_rates,
    std::vector<SurfaceSite>&& surface,
    std::unique_ptr<UserPunch>&& user_punch,
    std::unique_ptr<Output>&& output,
    std::unique_ptr<Dump>&& dump,
    Knobs&& knobs,
    MeshLib::PropertyVector<double>* porosity,
    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>&&
        media_map)
    : _phreeqc_input_file(project_file_name + "_phreeqc.inp"),
      _mesh(mesh),
      _database(std::move(database)),
      _chemical_system(std::move(chemical_system)),
      _reaction_rates(std::move(reaction_rates)),
      _surface(std::move(surface)),
      _user_punch(std::move(user_punch)),
      _output(std::move(output)),
      _dump(std::move(dump)),
      _knobs(std::move(knobs)),
      _porosity(porosity),
      _media_map(std::move(media_map))
{
    // initialize phreeqc instance
    if (CreateIPhreeqc() != phreeqc_instance_id)
    {
        OGS_FATAL(
            "Failed to initialize phreeqc instance, due to lack of memory.");
    }

    // load specified thermodynamic database
    if (LoadDatabase(phreeqc_instance_id, _database.c_str()) != IPQ_OK)
    {
        OGS_FATAL(
            "Failed in loading the specified thermodynamic database file: "
            "{:s}.",
            _database);
    }

    if (SetSelectedOutputFileOn(phreeqc_instance_id, 1) != IPQ_OK)
    {
        OGS_FATAL(
            "Failed to fly the flag for the specified file {:s} where phreeqc "
            "will write output.",
            _output->basic_output_setups.output_file);
    }

    if (_dump)
    {
        // Chemical composition of the aqueous solution of last time step will
        // be written into .dmp file once the second function argument is set to
        // one.
        SetDumpFileOn(phreeqc_instance_id, 1);
    }
}

void PhreeqcIO::initialize()
{
    num_chemical_systems = [&] {
        std::size_t num_chemical_systems = 0;
        for (std::size_t i = 0; i < _chemical_system_index_map.size(); ++i)
            num_chemical_systems += _chemical_system_index_map[i].size();
        return num_chemical_systems;
    }();

    _chemical_system->initialize(num_chemical_systems);
}

void PhreeqcIO::executeInitialCalculation(
    std::vector<GlobalVector> const& int_pt_x)
{
    setAqueousSolutions(int_pt_x);

    writeInputsToFile();

    execute();

    readOutputsFromFile();
}

std::vector<GlobalVector> PhreeqcIO::getIntPtProcessSolutions() const
{
    std::vector<GlobalVector> int_pt_x;

    auto const& aqueous_solution = _chemical_system->aqueous_solution;
    int_pt_x.push_back(*aqueous_solution->pressure);

    auto const& components = aqueous_solution->components;
    std::transform(components.begin(), components.end(),
                   std::back_inserter(int_pt_x),
                   [](auto const& c) { return *c.amount; });

    int_pt_x.push_back(*aqueous_solution->pH);

    return int_pt_x;
}

void PhreeqcIO::doWaterChemistryCalculation(
    std::vector<GlobalVector> const& int_pt_x, double const dt)
{
    setAqueousSolutions(int_pt_x);

    setAqueousSolutionsPrevFromDumpFile();

    writeInputsToFile(dt);

    execute();

    readOutputsFromFile();
}

std::vector<MeshLib::PropertyVector<double>*>
PhreeqcIO::getReactantVolumeFractionChange(MeshLib::Mesh& mesh)
{
    std::vector<MeshLib::PropertyVector<double>*> reactant_volume_change;

    if (!_chemical_system->kinetic_reactants.empty())
    {
        auto const& kinetic_reactants = _chemical_system->kinetic_reactants;
        std::transform(
            kinetic_reactants.begin(), kinetic_reactants.end(),
            std::back_inserter(reactant_volume_change), [&mesh](auto const& r) {
                return mesh.getProperties().template getPropertyVector<double>(
                    "volume_fraction_change_" + r.name,
                    MeshLib::MeshItemType::IntegrationPoint, 1);
            });
    }

    if (!_chemical_system->equilibrium_reactants.empty())
    {
        auto const& equilibrium_reactants =
            _chemical_system->equilibrium_reactants;
        std::transform(
            equilibrium_reactants.begin(),
            equilibrium_reactants.end(),
            std::back_inserter(reactant_volume_change),
            [&mesh](auto const& r) {
                return mesh.getProperties().template getPropertyVector<double>(
                    "volume_fraction_change_" + r.name,
                    MeshLib::MeshItemType::IntegrationPoint, 1);
            });
    }

    return reactant_volume_change;
}

void PhreeqcIO::setAqueousSolutions(std::vector<GlobalVector> const& int_pt_x)
{
    // pressure
    auto& aqueous_solution = _chemical_system->aqueous_solution;
    *aqueous_solution->pressure = int_pt_x.front();

    // components
    auto& components = aqueous_solution->components;
    for (unsigned component_id = 0; component_id < components.size();
         ++component_id)
    {
        *components[component_id].amount = int_pt_x[component_id + 1];
    }

    // pH
    *aqueous_solution->pH = int_pt_x.back();
}

void PhreeqcIO::setAqueousSolutionsPrevFromDumpFile()
{
    if (!_dump)
    {
        return;
    }

    auto const& dump_file = _dump->dump_file;
    std::ifstream in(dump_file);
    if (!in)
    {
        OGS_FATAL("Could not open phreeqc dump file '{:s}'.", dump_file);
    }

    std::size_t const num_chemical_systems = _mesh.getNumberOfBaseNodes();
    _dump->readDumpFile(in, num_chemical_systems);

    if (!in)
    {
        OGS_FATAL("Error when reading phreeqc dump file '{:s}'", dump_file);
    }

    in.close();
}

void PhreeqcIO::writeInputsToFile(double const dt)
{
    DBUG("Writing phreeqc inputs into file '{:s}'.", _phreeqc_input_file);
    std::ofstream out(_phreeqc_input_file, std::ofstream::out);

    if (!out)
    {
        OGS_FATAL("Could not open file '{:s}' for writing phreeqc inputs.",
                  _phreeqc_input_file);
    }

    print(out, dt);

    if (!out)
    {
        OGS_FATAL("Failed in generating phreeqc input file '{:s}'.",
                  _phreeqc_input_file);
    }

    out.close();
}

void PhreeqcIO::print(std::ostream& os, double const dt)
{
    os << _knobs << "\n";

    os << *_output << "\n";

    if (_user_punch)
    {
        os << *_user_punch << "\n";
    }

    if (!_reaction_rates.empty())
    {
        os << "RATES" << "\n";
        os << _reaction_rates << "\n";
    }

    for (std::size_t mesh_item_id = 0;
         mesh_item_id < _chemical_system_index_map.size();
         ++mesh_item_id)
    {
        auto const& local_chemical_system =
            _chemical_system_index_map[mesh_item_id];
        for (std::size_t ip = 0; ip < local_chemical_system.size(); ++ip)
        {
            auto const& chemical_system_id = local_chemical_system[ip];
            _chemical_system->print(mesh_item_id, chemical_system_id, os, dt,
                                    *_media_map,
                                    (*_porosity)[chemical_system_id]);
        }
    }

    if (_dump)
    {
        //        dump->print(os, num_chemical_systems);
    }
}

void PhreeqcIO::execute()
{
    INFO("Phreeqc: Executing chemical calculation.");
    if (RunFile(phreeqc_instance_id, _phreeqc_input_file.c_str()) != IPQ_OK)
    {
        OutputErrorString(phreeqc_instance_id);
        OGS_FATAL(
            "Failed in performing speciation calculation with the generated "
            "phreeqc input file '{:s}'.",
            _phreeqc_input_file);
    }
}

void PhreeqcIO::readOutputsFromFile()
{
    auto const& basic_output_setups = _output->basic_output_setups;
    auto const& phreeqc_result_file = basic_output_setups.output_file;
    DBUG("Reading phreeqc results from file '{:s}'.", phreeqc_result_file);
    std::ifstream in(phreeqc_result_file);

    if (!in)
    {
        OGS_FATAL("Could not open phreeqc result file '{:s}'.",
                  phreeqc_result_file);
    }

    read(in);

    if (!in)
    {
        OGS_FATAL("Error when reading phreeqc result file '{:s}'",
                  phreeqc_result_file);
    }

    in.close();
}

void PhreeqcIO::read(std::istream& in)
{
    // Skip the headline
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::string line;
    auto const& dropped_item_ids = _output->dropped_item_ids;

    int const num_skipped_lines = _surface.empty() ? 1 : 2;

    for (std::size_t mesh_item_id = 0;
         mesh_item_id < _chemical_system_index_map.size();
         ++mesh_item_id)
    {
        auto const& local_chemical_system =
            _chemical_system_index_map[mesh_item_id];
        for (std::size_t ip = 0; ip < local_chemical_system.size(); ++ip)
        {
            auto const& chemical_system_id = local_chemical_system[ip];

            // Skip equilibrium calculation result of initial solution
            for (int i = 0; i < num_skipped_lines; ++i)
            {
                in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }

            // Get calculation result of the solution after the reaction
            if (!std::getline(in, line))
            {
                OGS_FATAL(
                    "Error when reading calculation result of Solution {:d} "
                    "after the reaction.",
                    chemical_system_id + 1);
            }

            std::vector<double> accepted_items;
            std::vector<std::string> items;
            boost::trim_if(line, boost::is_any_of("\t "));
            boost::algorithm::split(items, line, boost::is_any_of("\t "),
                                    boost::token_compress_on);
            for (int item_id = 0; item_id < static_cast<int>(items.size());
                 ++item_id)
            {
                if (std::find(dropped_item_ids.begin(), dropped_item_ids.end(),
                              item_id) == dropped_item_ids.end())
                {
                    double value;
                    try
                    {
                        value = std::stod(items[item_id]);
                    }
                    catch (const std::invalid_argument& e)
                    {
                        OGS_FATAL(
                            "Invalid argument. Could not convert string '{:s}' "
                            "to "
                            "double for chemical system {:d}, column {:d}. "
                            "Exception "
                            "'{:s}' was thrown.",
                            items[item_id], chemical_system_id + 1, item_id,
                            e.what());
                    }
                    catch (const std::out_of_range& e)
                    {
                        OGS_FATAL(
                            "Out of range error. Could not convert string "
                            "'{:s}' "
                            "to "
                            "double for chemical system {:d}, column {:d}. "
                            "Exception "
                            "'{:s}' was thrown.",
                            items[item_id], chemical_system_id + 1, item_id,
                            e.what());
                    }
                    accepted_items.push_back(value);
                }
            }
            assert(accepted_items.size() == _output->accepted_items.size());

            ParameterLib::SpatialPosition pos;
            MaterialPropertyLib::VariableArray vars;

            auto const& medium = *_media_map->getMedium(mesh_item_id);
            auto const& phase = medium.phase("AqueousLiquid");

            auto const density =
                phase.property(MaterialPropertyLib::PropertyType::density)
                    .template value<double>(vars, pos, 0, 0);

            for (int item_id = 0;
                 item_id < static_cast<int>(accepted_items.size());
                 ++item_id)
            {
                auto const& accepted_item = _output->accepted_items[item_id];
                auto const& item_name = accepted_item.name;

                auto compare_by_name = [&item_name](auto const& item) {
                    return item.name == item_name;
                };

                switch (accepted_item.item_type)
                {
                    case ItemType::pH:
                    {
                        // Update pH value
                        (*_chemical_system->aqueous_solution
                              ->pH)[chemical_system_id] =
                            std::pow(10, -accepted_items[item_id]);
                        break;
                    }
                    case ItemType::pe:
                    {
                        // Update pe value
                        //                    (*aqueous_solution.pe)[chemical_system_id]
                        //                    =
                        //                        accepted_items[item_id];
                        break;
                    }
                    case ItemType::Component:
                    {
                        // Update component concentrations
                        auto& component = BaseLib::findElementOrError(
                            _chemical_system->aqueous_solution->components
                                .begin(),
                            _chemical_system->aqueous_solution->components
                                .end(),
                            compare_by_name,
                            "Could not find component '" + item_name + "'.");
                        component.amount->set(
                            chemical_system_id,
                            accepted_items[item_id] * density);
                        break;
                    }
                    case ItemType::EquilibriumReactant:
                    {
                        // Update amounts of equilibrium reactant
                        auto& equilibrium_reactant =
                            BaseLib::findElementOrError(
                                _chemical_system->equilibrium_reactants.begin(),
                                _chemical_system->equilibrium_reactants.end(),
                                compare_by_name,
                                "Could not find equilibrium reactant '" +
                                    item_name + "'.");
                        // TODO: unit conversion
                        (*equilibrium_reactant.mass)[chemical_system_id] +=
                            accepted_items[item_id];
                        (*equilibrium_reactant.mass_loss)[chemical_system_id] =
                            accepted_items[item_id];
                        break;
                    }
                    case ItemType::KineticReactant:
                    {
                        // Update amounts of kinetic reactants
                        auto& kinetic_reactant = BaseLib::findElementOrError(
                            _chemical_system->kinetic_reactants.begin(),
                            _chemical_system->kinetic_reactants.end(),
                            compare_by_name,
                            "Could not find kinetic reactant '" + item_name +
                                "'.");
                        (*kinetic_reactant
                              .volume_fraction)[chemical_system_id] +=
                            accepted_items[item_id] *
                            kinetic_reactant.molar_volume *
                            (*_porosity)[chemical_system_id] * density;
                        (*kinetic_reactant
                              .volume_fraction_change)[chemical_system_id] =
                            accepted_items[item_id] *
                            kinetic_reactant.molar_volume *
                            (*_porosity)[chemical_system_id] * density;
                        break;
                    }
                    case ItemType::SecondaryVariable:
                    {
                        assert(_user_punch);
                        // Update values of secondary variables
                        auto& secondary_variable = BaseLib::findElementOrError(
                            _user_punch->secondary_variables.begin(),
                            _user_punch->secondary_variables.end(),
                            compare_by_name,
                            "Could not find secondary variable '" + item_name +
                                "'.");
                        //                    (*secondary_variable.value)[global_id]
                        //                    =
                        //                            accepted_items[item_id];
                        break;
                    }
                }
            }
        }
    }
}

std::vector<std::string> PhreeqcIO::getComponentList() const
{
    std::vector<std::string> component_names;
    auto const& aqueous_solution = _chemical_system->aqueous_solution;
    auto const& components = aqueous_solution->components;
    std::transform(components.begin(), components.end(),
                   std::back_inserter(component_names),
                   [](auto const& c) { return c.name; });

    component_names.push_back("H");

    return component_names;
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
