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

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "ChemicalSolverType.h"

namespace BaseLib
{
class ConfigTree;
}

namespace MeshLib
{
class Mesh;
}

namespace MaterialPropertyLib
{
class Medium;
}

namespace ChemistryLib
{
class ChemicalSolverInterface;

template <ChemicalSolver chemical_solver>
std::shared_ptr<ChemicalSolverInterface> createChemicalSolverInterface(
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
    BaseLib::ConfigTree const& config, std::string const& output_directory,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);
}  // namespace ChemistryLib
