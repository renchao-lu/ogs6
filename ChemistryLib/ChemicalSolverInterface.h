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

#include <vector>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "MeshLib/PropertyVector.h"

namespace MeshLib
{
class Mesh;
}

namespace ChemistryLib
{
class ChemicalSolverInterface
{
public:
    virtual void initialize() {}

    virtual void executeInitialCalculation(
        std::vector<GlobalVector> const& int_pt_x) = 0;

    virtual void doWaterChemistryCalculation(
        std::vector<GlobalVector> const& int_pt_x, double const dt) = 0;

    virtual ~ChemicalSolverInterface() = default;

    virtual std::vector<GlobalVector> getIntPtProcessSolutions() const
    {
        return {};
    }
    virtual std::vector<std::string> getComponentList() const { return {}; }

    virtual std::vector<std::vector<GlobalIndexType>>&
    getChemicalSystemIndexMap() = 0;
};
}  // namespace ChemistryLib
