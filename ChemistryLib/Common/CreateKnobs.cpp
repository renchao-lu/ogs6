/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateKnobs.h"
#include "Knobs.h"
#include "BaseLib/ConfigTree.h"

namespace ChemistryLib
{
Knobs createKnobs(BaseLib::ConfigTree const& config)
{
    auto const max_iterations =
        //! \ogs_file_param{prj__chemical_system__knobs__max_iter}
        config.getConfigParameter<int>("max_iter");

    auto const relative_convergence_tolerance =
        //! \ogs_file_param{prj__chemical_system__knobs__relative_convergence_tolerance}
        config.getConfigParameter<double>("relative_convergence_tolerance");

    auto const tolerance =
        //! \ogs_file_param{prj__chemical_system__knobs__tolerance}
        config.getConfigParameter<double>("tolerance");

    auto const step_size =
        //! \ogs_file_param{prj__chemical_system__knobs__step_size}
        config.getConfigParameter<int>("step_size");

    auto const scaling =
        //! \ogs_file_param{prj__chemical_system__knobs__scaling}
        config.getConfigParameter<bool>("scaling");

    return {max_iterations, relative_convergence_tolerance, tolerance,
            step_size, scaling};
}
}  // namespace ChemistryLib
