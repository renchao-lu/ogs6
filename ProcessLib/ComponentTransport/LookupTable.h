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

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "BaseLib/Error.h"

namespace
{
std::pair<double, double> getNeighboringValues(
        std::vector<double> const& vec, double value)
{
    auto const it = std::lower_bound(vec.begin(), vec.end(), value);

    if (it == vec.begin())
    {
        // negative concentration?
        return std::make_pair(*it, *(it+1));
    }
    else if (it == vec.end())
    {
        return std::make_pair(*vec.end(), *(vec.end()-1));
    }
    else
    {
        return std::make_pair(*(it-1), *it);
    }
}

std::size_t getIntersectedIndex(
        std::vector<std::size_t> const& vec1,
        std::vector<std::size_t> const& vec2)
{
    // intersect vectors
    std::vector<std::size_t> vec;
    std::set_intersection(
                vec1.begin(), vec1.end(), vec2.begin(), vec2.end(),
                std::back_inserter(vec));

    return vec[0];
}
}  // namespace

namespace ProcessLib
{
namespace ComponentTransport
{
struct LookupTable
{
    LookupTable(
        std::map<std::string, std::map<double, std::vector<std::size_t>>>&&
            radionuclides_concentrations_,
        std::map<std::string, std::vector<double>>&& kd_matrix_)
        : radionuclides_concentrations(
              std::move(radionuclides_concentrations_)),
          kd_matrix(std::move(kd_matrix_))
    {
    }

    void lookup(std::vector<GlobalVector*> const& x,
                std::vector<std::unique_ptr<GlobalVector>> const& _x_previous_timestep,
                std::vector<std::pair<int, std::string>> const& process_id_to_component_name_map)
    {
        auto const process_id = 2;
        auto const n_nodes = x[process_id]->size();
        for (auto node_id = 0; node_id < n_nodes; ++node_id)
        {
            auto node_value = x[process_id]->get(node_id);
            auto node_value_previous_timestep =
                    _x_previous_timestep[process_id]->get(node_id);
            auto const neighboring_values = getNeighboringValues(uranium, node_value);
            auto const neighboring_values_previous_timestep =
                    getNeighboringValues(uranium, node_value_previous_timestep);

            // look up the table
            auto const base_point_index = getIntersectedIndex(uranium_cur[neighboring_values.first],
                    uranium_prev[neighboring_values_previous_timestep.first]);
            auto& base_value = kd_matrix.at("U(6)")[base_point_index];

            // linear approximation
            auto const point_1_index = getIntersectedIndex(uranium_cur[neighboring_values.second],
                    uranium_prev[neighboring_values_previous_timestep.first]);
            auto& point_1_value = kd_matrix.at("U(6)")[point_1_index];
            auto df_x = (point_1_value - base_value) / (neighboring_values.second
                                                        - neighboring_values.first);

            auto const point_2_index = getIntersectedIndex(uranium_cur[neighboring_values.first],
                    uranium_prev[neighboring_values_previous_timestep.second]);
            auto& point_2_value = kd_matrix.at("U(6)")[point_2_index];
            auto df_y = (point_2_value - base_value) /
                    (neighboring_values_previous_timestep.second
                                  - neighboring_values_previous_timestep.first);

            auto new_value = base_value + df_x * (node_value - neighboring_values.first)
                    + df_y * (node_value_previous_timestep
                              - neighboring_values_previous_timestep.first);

            x[process_id]->set(node_id, new_value);
        }
    }

    std::map<std::string, std::vector<double>> const kd_matrix;
    std::map<std::string, std::map<double, std::vector<std::size_t>>>
        radionuclides_concentrations;
};
}  // namespace ComponentTransport
}  // namespace ProcessLib
