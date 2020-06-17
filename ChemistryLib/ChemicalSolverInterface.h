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

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

namespace ChemistryLib
{
class ChemicalSolverInterface
{
public:
    virtual void executeInitialCalculation(
        std::vector<GlobalVector*>& process_solutions) = 0;

    virtual void doWaterChemistryCalculation(
        std::vector<GlobalVector*>& process_solutions, double const dt) = 0;

    virtual std::vector<std::string> const getComponentList() const
    {
        return {};
    }

    std::vector<std::vector<GlobalIndexType>>& getChemicalSystemIndexMap()
    {
        return _chemical_system_index_map;
    }

    virtual ~ChemicalSolverInterface() = default;

protected:
    std::vector<std::vector<GlobalIndexType>> _chemical_system_index_map;
};
}  // namespace ChemistryLib
