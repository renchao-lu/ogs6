/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <map>
#include <vector>

#include "ChemicalSolverInterface.h"
#include "PhreeqcKernelData/EquilibriumPhase.h"
#include "PhreeqcKernelData/KineticReactant.h"
#include "PhreeqcKernelData/Surface.h"

#include "ThirdParty/iphreeqc/src/src/phreeqcpp/Phreeqc.h"

class cxxSolution;
class cxxISolution;

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
class AqueousSolution;
class ReactionRate;

class PhreeqcKernel final : public ChemicalSolverInterface, private Phreeqc
{
public:
    PhreeqcKernel(std::size_t const num_chemical_systems,
                  std::vector<std::pair<int, std::string>> const&
                      process_id_to_component_name_map,
                  std::string const database,
                  AqueousSolution aqueous_solution,
                  std::unique_ptr<Equilibriums>&& equilibrium_phases,
                  std::unique_ptr<Kinetics>&& kinetic_reactants,
                  std::vector<ReactionRate>&& reaction_rates,
                  std::unique_ptr<Surface>&& surface);

    void executeInitialCalculation(
        std::vector<GlobalVector*>& process_solutions) override;

    void doWaterChemistryCalculation(
        std::vector<GlobalVector*>& process_solutions,
        double const dt) override;

    void setAqueousSolutions(
        std::vector<GlobalVector*> const& process_solutions);

    void execute(std::vector<GlobalVector*>& process_solutions);

    void updateNodalProcessSolutions(
        std::vector<GlobalVector*> const& process_solutions,
        std::size_t const node_id);

private:
    void initializePhreeqcGeneralSettings() { do_initialize(); }

    void loadDatabase(std::string const& database);

    void reinitializeRates();

    void getElementsInSpecies(std::string formula, double moles)
    {
        count_elts = 0;
        paren_count = 0;
        char *formula_char = string_duplicate(formula.c_str());
        get_elts_in_species(&formula_char, moles);
    }

    char* getElement(std::string formula)
    {
        char *formula_char = string_duplicate(formula.c_str());
        char *name = string_duplicate(formula.c_str());
        name[0] = '\0';
        int l = formula.size();

        get_elt(&formula_char, name, &l);
        return name;
    }

    cxxNameDouble getElementListNameDouble()
    {
        return elt_list_NameDouble();
    }

    void setConvergenceTolerance()
    {
        convergence_tolerance = 1e-12;

        // knobs
        {
            itmax = 250;
            convergence_tolerance = 1e-6;
            ineq_tol = 1e-20;
            step_size = 5;
            diagonal_scale = 1;
        }
    }

    void configureOutputSettings() { pr.all = false; }

    cxxISolution* getOrCreateInitialAqueousSolution(
        cxxSolution& aqueous_solution);

    bool isHydrogen(char const* element) const
    {
        return strcmp(element, "H") == 0;
    }

    void setTimeStepSize(double const dt);

    std::map<int, struct master*> _process_id_to_master_map;
    std::unique_ptr<cxxISolution const> _initial_aqueous_solution;
    std::unique_ptr<cxxSolution const> _aqueous_solution;
    std::unique_ptr<cxxSurface const> _surface;
    std::vector<ReactionRate> const _reaction_rates;
    bool initial_step = true;
};
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
