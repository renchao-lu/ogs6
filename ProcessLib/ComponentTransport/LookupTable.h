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

namespace ProcessLib
{
namespace ComponentTransport
{
struct LookupTable
{
    LookupTable(std::map<std::string, std::vector<double>> input_parameters_,
                std::map<std::string, std::vector<double>>
                    kd_matrix_)
        : input_parameters(std::move(input_parameters_)),
          kd_matrix(std::move(kd_matrix_))
    {
    }

    std::map<std::string, std::vector<double>> const input_parameters;
    std::map<std::string, std::vector<double>> const kd_matrix;
};
}  // namespace ComponentTransport
}  // namespace ProcessLib
