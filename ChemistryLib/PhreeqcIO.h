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

#include "ChemicalSolverInterface.h"
#include "PhreeqcIOData/Knobs.h"

namespace MeshLib
{
class Mesh;
}

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct ChemicalSystem;
struct ReactionRate;
struct Output;
struct SurfaceSite;
struct Dump;
struct UserPunch;

enum class Status
{
    SettingAqueousSolutions,
    UpdatingProcessSolutions
};

class PhreeqcIO final : public ChemicalSolverInterface
{
public:
    PhreeqcIO(std::string const project_file_name,
              MeshLib::Mesh const& mesh,
              std::string&& database,
              std::unique_ptr<ChemicalSystem>&& chemical_system,
              std::vector<ReactionRate>&& reaction_rates,
              std::vector<SurfaceSite>&& surface,
              std::unique_ptr<UserPunch>&& user_punch,
              std::unique_ptr<Output>&& output,
              std::unique_ptr<Dump>&& dump,
              Knobs&& knobs);

    void executeInitialCalculation(
        std::vector<GlobalVector*>& process_solutions) override;

    void doWaterChemistryCalculation(
        std::vector<GlobalVector*>& process_solutions,
        double const dt) override;

    void setAqueousSolutionsOrUpdateProcessSolutions(
        std::vector<GlobalVector*> const& process_solutions,
        Status const status);

    void writeInputsToFile(double const dt = 0);

    void execute();

    void readOutputsFromFile();

    friend std::ostream& operator<<(std::ostream& os,
                                    PhreeqcIO const& phreeqc_io);

    friend std::istream& operator>>(std::istream& in, PhreeqcIO& phreeqc_io);

    std::vector<std::string> const getComponentList() const override;

    std::string const _phreeqc_input_file;

private:
    PhreeqcIO& operator<<(double const dt)
    {
        _dt = dt;
        return *this;
    }

    void setAqueousSolutionsPrevFromDumpFile();

    MeshLib::Mesh const& _mesh;
    std::string const _database;
    std::unique_ptr<ChemicalSystem> _chemical_system;
    std::vector<ReactionRate> const _reaction_rates;
    std::vector<SurfaceSite> const _surface;
    std::unique_ptr<UserPunch> _user_punch;
    std::unique_ptr<Output> const _output;
    std::unique_ptr<Dump> const _dump;
    Knobs const _knobs;
    double _dt = std::numeric_limits<double>::quiet_NaN();
    const int phreeqc_instance_id = 0;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
