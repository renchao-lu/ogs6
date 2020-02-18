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
        return std::make_pair(*it, *(it+1));
    }
    else if (it == vec.end())
    {
        return std::make_pair(*(vec.rbegin() + 1), *vec.rbegin());
    }
    else
    {
        return std::make_pair(*(it-1), *it);
    }
}

std::size_t getIntersectedIndex(
    std::map<std::string, std::map<double, std::vector<std::size_t>>>& map,
    std::map<std::string, int> const& concentration_field_to_process_id,
    std::vector<double>
        values,
    std::vector<double>
        values_previous_timestep)
{
    std::vector<std::size_t> vec1;
    auto map_values = map[concentration_field_to_process_id.begin()->first];
    std::vector<std::size_t> vec = map_values[values[0]];
    std::vector<std::size_t> temp_vec;

    int i = 0;
    for (auto const& pair : concentration_field_to_process_id)
    {
        map_values = map[pair.first];
        vec1 = map_values[values[i]];

        // intersect vectors
        std::set_intersection(vec1.begin(), vec1.end(), vec.begin(), vec.end(),
                              std::back_inserter(temp_vec));

        std::swap(vec, temp_vec);

        temp_vec.clear();
        ++i;
    }

    i = 0;
    for (auto const& pair : concentration_field_to_process_id)
    {
        map_values = map[pair.first + "_prev"];
        vec1 = map_values[values_previous_timestep[i]];

        // intersect vectors
        std::set_intersection(vec1.begin(), vec1.end(), vec.begin(), vec.end(),
                              std::back_inserter(temp_vec));

        std::swap(vec, temp_vec);

        temp_vec.clear();
        ++i;
    }

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
        std::map<std::string, int>&& concentration_field_to_process_id_,
        std::map<std::string, std::vector<double>>&& concentration_seeds_,
        std::map<std::string, std::vector<double>>&& surface_field_seeds_,
        std::map<std::string, std::map<double, std::vector<std::size_t>>>&&
            radionuclides_concentrations_,
        std::map<std::string, std::vector<double>>&& kd_matrix_)
        : concentration_field_to_process_id(
              std::move(concentration_field_to_process_id_)),
          concentration_seeds(std::move(concentration_seeds_)),
          surface_field_seeds(std::move(surface_field_seeds_)),
          radionuclides_concentrations(
              std::move(radionuclides_concentrations_)),
          kd_matrix(std::move(kd_matrix_))
    {
    }

    void lookup(
        std::vector<GlobalVector*> const& x,
        std::vector<std::unique_ptr<GlobalVector>> const& _x_previous_timestep)
    {
        auto const n_nodes = x[0]->size();
        for (auto node_id = 0; node_id < n_nodes; ++node_id)
        {
            std::vector<double> node_values;
            std::vector<double> node_values_previous_timestep;
            std::vector<std::pair<double, double>> neighboring_values;
            std::vector<std::pair<double, double>>
                neighboring_values_previous_timestep;
            for (auto const& pair : concentration_field_to_process_id)
            {
                auto const node_value = x[pair.second]->get(node_id);
                // how to deal with negative concentration
                auto const neighboring_value = getNeighboringValues(
                    concentration_seeds[pair.first], node_value);
                node_values.push_back(node_value);
                neighboring_values.push_back(neighboring_value);

                auto const node_value_previous_timestep =
                    _x_previous_timestep[pair.second]->get(node_id);
                auto const neighboring_value_previous_timestep =
                    getNeighboringValues(
                        surface_field_seeds[pair.first + "_prev"],
                        node_value_previous_timestep);
                node_values_previous_timestep.push_back(
                    node_value_previous_timestep);
                neighboring_values_previous_timestep.push_back(
                    neighboring_value_previous_timestep);
            }

            // look up the table
            std::vector<double> values;
            std::vector<double> values_previous_timestep;

            for (auto const& neighboring_value : neighboring_values)
            {
                values.push_back(neighboring_value.first);
            }
            for (auto const& neighboring_value_previous_timestep :
                 neighboring_values_previous_timestep)
            {
                values_previous_timestep.push_back(
                    neighboring_value_previous_timestep.first);
            }

            auto const base_point_index = getIntersectedIndex(
                radionuclides_concentrations, concentration_field_to_process_id,
                values, values_previous_timestep);

            for (auto const& pair : concentration_field_to_process_id)
            {
                auto base_value =
                    kd_matrix.at(pair.first + "_new")[base_point_index];
                auto new_value = base_value;

                // linear approximation
                // loop over values
                for (auto i = 0; i < values.size(); ++i)
                {
                    values[i] = neighboring_values[i].second;
                    auto const interpolation_point_index =
                        getIntersectedIndex(radionuclides_concentrations,
                                            concentration_field_to_process_id,
                                            values, values_previous_timestep);
                    auto& interpolation_point_value = kd_matrix.at(
                        pair.first + "_new")[interpolation_point_index];
                    auto df = (interpolation_point_value - base_value) /
                              (neighboring_values[i].second -
                               neighboring_values[i].first);

                    new_value +=
                        df * (node_values[i] - neighboring_values[i].first);
                    values[i] = neighboring_values[i].first;
                }

                // loop over values_previous_timestep
                for (auto i = 0; i < values_previous_timestep.size(); ++i)
                {
                    values_previous_timestep[i] =
                        neighboring_values_previous_timestep[i].second;
                    auto const interpolation_point_index =
                        getIntersectedIndex(radionuclides_concentrations,
                                            concentration_field_to_process_id,
                                            values, values_previous_timestep);
                    auto& interpolation_point_value = kd_matrix.at(
                        pair.first + "_new")[interpolation_point_index];
                    auto df = (interpolation_point_value - base_value) /
                              (neighboring_values_previous_timestep[i].second -
                               neighboring_values_previous_timestep[i].first);

                    new_value +=
                        df * (node_values_previous_timestep[i] -
                              neighboring_values_previous_timestep[i].first);
                    values[i] = neighboring_values[i].first;
                }

                x[pair.second]->set(node_id, new_value);
            }
        }
    }

    std::map<std::string, int> concentration_field_to_process_id;
    std::map<std::string, std::vector<double>> concentration_seeds;
    std::map<std::string, std::vector<double>> surface_field_seeds;
    std::map<std::string, std::map<double, std::vector<std::size_t>>>
        radionuclides_concentrations;
    std::map<std::string, std::vector<double>> const kd_matrix;
};
}  // namespace ComponentTransport
}  // namespace ProcessLib
