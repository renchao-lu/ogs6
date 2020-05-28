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

#include "ChemicalSolverInterface.h"
#include "PhreeqcIOData/Knobs.h"

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
struct ChemicalSystem;
struct KineticReactant;
struct ReactionRate;
struct Output;
struct SurfaceSite;
struct Dump;
struct UserPunch;

class PhreeqcIO final : public ChemicalSolverInterface
{
public:
    PhreeqcIO(
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
            media_map);

    void initialize() override;

    void executeInitialCalculation(
        std::vector<GlobalVector> const& int_pt_x) override;

    void doWaterChemistryCalculation(std::vector<GlobalVector> const& int_pt_x,
                                     double const dt) override;

    void setAqueousSolutions(
        std::vector<GlobalVector> const& ip_transport_solutions);

    void writeInputsToFile(double const dt = 0);

    void execute();

    void readOutputsFromFile();

    std::vector<GlobalVector> getIntPtProcessSolutions() const override;

    std::vector<MeshLib::PropertyVector<double>*>
    getReactantVolumeFractionChange(MeshLib::Mesh& mesh) override;

    MeshLib::PropertyVector<double>* getPorosity() override
    {
        return _porosity;
    }

    std::vector<std::string> getComponentList() const;

    std::vector<std::vector<GlobalIndexType>>& getChemicalSystemIndexMap()
    {
        return _chemical_system_index_map;
    }

    friend std::istream& operator>>(std::istream& in, PhreeqcIO& phreeqc_io);

    std::string const _phreeqc_input_file;

private:
    void print(std::ostream& os, double const dt);

    void read(std::istream& in);

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
    MeshLib::PropertyVector<double>* _porosity;
    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        _media_map;
    const int phreeqc_instance_id = 0;
    std::vector<std::vector<GlobalIndexType>> _chemical_system_index_map;
    std::size_t num_chemical_systems =
        std::numeric_limits<std::size_t>::quiet_NaN();
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
