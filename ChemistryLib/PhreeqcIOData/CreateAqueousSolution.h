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

namespace BaseLib
{
class ConfigTree;
}

namespace MeshLib
{
class Mesh;
}

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct AqueousSolution;

std::unique_ptr<AqueousSolution> createAqueousSolution(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& mesh);
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
