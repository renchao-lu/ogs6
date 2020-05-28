/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ChemicalSystem.h"
#include "AqueousSolution.h"

#include "MaterialLib/MPL/Medium.h"
#include "MathLib/LinAlg/MatrixVectorTraits.h"
#include "MathLib/LinAlg/UnifiedMatrixSetters.h"
#include "MeshLib/Mesh.h"
#include "ParameterLib/SpatialPosition.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
void ChemicalSystem::initialize(std::size_t const num_chemical_systems)
{
    aqueous_solution->pressure =
        MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
            num_chemical_systems);

    aqueous_solution->pH =
        MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
            num_chemical_systems);

    aqueous_solution->pe->resize(num_chemical_systems);
    std::fill(aqueous_solution->pe->begin(), aqueous_solution->pe->end(),
              aqueous_solution->pe0);

    auto& components = aqueous_solution->components;
    for (auto& component : components)
    {
        component.amount =
            MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
                num_chemical_systems);
    }

    // kinetic reactants
    if (!kinetic_reactants.empty())
    {
        for (auto& kinetic_reactant : kinetic_reactants)
        {
            // volume fraction
            kinetic_reactant.volume_fraction->resize(num_chemical_systems);
            std::fill(std::begin(*kinetic_reactant.volume_fraction),
                      std::end(*kinetic_reactant.volume_fraction),
                      kinetic_reactant.initial_volume_fraction);

            // volume fraction change
            kinetic_reactant.volume_fraction_change->resize(
                num_chemical_systems);
            std::fill(std::begin(*kinetic_reactant.volume_fraction_change),
                      std::end(*kinetic_reactant.volume_fraction_change),
                      std::numeric_limits<double>::quiet_NaN());
        }
    }
}

void ChemicalSystem::print(
    std::size_t const mesh_item_id, std::size_t const chemical_system_id,
    std::ostream& os, double const dt,
    MaterialPropertyLib::MaterialSpatialDistributionMap const& media_map,
    double const& porosity)
{
    ParameterLib::SpatialPosition pos;
    MaterialPropertyLib::VariableArray vars;

    auto const& medium = *media_map.getMedium(mesh_item_id);
    auto const& phase = medium.phase("AqueousLiquid");

    auto const density =
        phase.property(MaterialPropertyLib::PropertyType::density)
            .template value<double>(vars, pos, 0, dt);

    os << "SOLUTION " << chemical_system_id + 1 << "\n";
    aqueous_solution->print(os, chemical_system_id, density);

    //        if (_dump)
    //        {
    //            auto const& aqueous_solutions_prev =
    //            _dump->aqueous_solutions_prev; if
    //            (!aqueous_solutions_prev.empty())
    //            {
    //                //                    os <<
    //                aqueous_solutions_prev[local_id] <<
    //                //                    "\n\n";
    //            }
    //        }

    os << "USE solution none"
       << "\n";
    os << "END"
       << "\n\n";

    os << "USE solution " << chemical_system_id + 1 << "\n\n";

    if (!equilibrium_reactants.empty())
    {
        os << "EQUILIBRIUM_PHASES " << chemical_system_id + 1 << "\n";
        for (auto const& equilibrium_reactant : equilibrium_reactants)
        {
            equilibrium_reactant.print(os, chemical_system_id, porosity,
                                       density);
        }
        os << "\n";
    }

    if (!kinetic_reactants.empty())
    {
        os << "KINETICS " << chemical_system_id + 1 << "\n";
        for (auto const& kinetic_reactant : kinetic_reactants)
        {
            kinetic_reactant.print(os, chemical_system_id, porosity, density);
        }
        os << "-steps " << dt << "\n"
           << "\n";
    }

    //    if (!_surface.empty())
    //    {
    //        //                os << "SURFACE " << global_id + 1 << "\n";
    //        //                std::size_t aqueous_solution_id =
    //        //                        dump->aqueous_solutions_prev.empty()
    //        //                        ? global_id + 1
    //        //                        : num_chemical_systems + global_id +
    //        //                        1;
    //        //                os << "-equilibrate with solution " <<
    //        //                aqueous_solution_id << "\n"; os <<
    //        //                "-sites_units DENSITY"
    //        //                   << "\n";
    //        //                os << surface << "\n";
    //        //                os << "SAVE solution " << global_id + 1 <<
    //        //                "\n";
    //    }

    os << "END"
       << "\n\n";
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
