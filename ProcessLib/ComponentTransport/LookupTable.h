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

#include "BaseLib/Algorithm.h"
#include "BaseLib/Error.h"
#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

namespace ProcessLib
{
namespace ComponentTransport
{
struct Field
{
    Field(std::string field_, std::vector<double> interpolation_points_,
          std::vector<std::vector<std::size_t>> indices_)
        : field(field_),
          interpolation_points(interpolation_points_),
          indices(indices_)
    {
    }

    std::pair<double, double> getNeighboringInterpolationPoints(double value)
    {
        auto const it = std::lower_bound(interpolation_points.begin(),
                                         interpolation_points.end(), value);

        if (it == interpolation_points.begin())
        {
            return std::make_pair(*it, *(it + 1));
        }
        else if (it == interpolation_points.end())
        {
            return std::make_pair(*(interpolation_points.rbegin() + 1),
                                  *interpolation_points.rbegin());
        }
        else
        {
            return std::make_pair(*(it - 1), *it);
        }
    }

    std::string field;
    std::vector<double> interpolation_points;
    std::vector<std::vector<std::size_t>> indices;
};

struct LookupTable
{
    LookupTable(std::vector<std::pair<std::string, int>>&&
                    concentration_field_to_process_id_,
                std::vector<Field>&& fields_,
                std::map<std::string, std::vector<double>>&& table_)
        : concentration_field_to_process_id(
              std::move(concentration_field_to_process_id_)),
          fields(std::move(fields_)),
          table(std::move(table_))
    {
    }

    void lookup(
        std::vector<GlobalVector*> const& x,
        std::vector<std::unique_ptr<GlobalVector>> const& _x_previous_timestep)
    {
        auto const n_nodes = x[0]->size();
        for (auto node_id = 0; node_id < n_nodes; ++node_id)
        {
            std::vector<double> nodal_x;
            for (auto const& pair : concentration_field_to_process_id) {
                nodal_x.emplace_back(x[pair.second]->get(node_id));
            }
            for (auto const& pair : concentration_field_to_process_id) {
                nodal_x.emplace_back(
                    _x_previous_timestep[pair.second]->get(node_id));
            }

            // how to deal with negative concentration
            std::vector<std::pair<double, double>>
                neighboring_interpolation_points;
            for (auto i = 0; i < nodal_x.size(); ++i)
            {
                auto ip =
                    fields[i].getNeighboringInterpolationPoints(nodal_x[i]);
                neighboring_interpolation_points.push_back(ip);
            }

            // look up the table
            std::vector<double> reference_points;
            for (auto const& neighboring_interpolation_point :
                 neighboring_interpolation_points)
            {
                reference_points.push_back(
                    neighboring_interpolation_point.first);
            }

            auto const base_point_index =
                getIntersectedIndex(fields, reference_points);

            for (auto const& pair : concentration_field_to_process_id)
            {
                auto base_value =
                    table.at(pair.first + "_new")[base_point_index];
                auto new_value = base_value;

                // linear approximation
                for (auto i = 0; i < fields.size(); ++i)
                {
                    reference_points[i] =
                        neighboring_interpolation_points[i].second;
                    auto const interpolation_point_index =
                        getIntersectedIndex(fields, reference_points);
                    auto& interpolation_point_value = table.at(
                        pair.first + "_new")[interpolation_point_index];
                    auto slope = (interpolation_point_value - base_value) /
                                 (neighboring_interpolation_points[i].second -
                                  neighboring_interpolation_points[i].first);

                    new_value +=
                        slope * (nodal_x[i] -
                                 neighboring_interpolation_points[i].first);
                    reference_points[i] =
                        neighboring_interpolation_points[i].first;
                }

                x[pair.second]->set(node_id, new_value);
            }
        }
    }

    std::size_t getIntersectedIndex(
        std::vector<Field> const& fields, std::vector<double> values)
    {
        std::vector<std::size_t> vec1;
        std::vector<std::size_t> temp_vec;
        std::vector<std::size_t> vec = fields[0].indices[BaseLib::findIndex(
            fields[0].interpolation_points, values[0])];

        auto num_variables = values.size();
        for (auto i = 0; i < num_variables; ++i)
        {
            auto index =
                BaseLib::findIndex(fields[i].interpolation_points, values[i]);
            vec1 = fields[i].indices[index];
            // intersect vectors
            std::set_intersection(vec1.begin(), vec1.end(), vec.begin(),
                                  vec.end(), std::back_inserter(temp_vec));

            std::swap(vec, temp_vec);
            temp_vec.clear();
        }

        return vec[0];
    }

    std::vector<std::pair<std::string, int>> concentration_field_to_process_id;
    std::vector<Field> fields;
    std::map<std::string, std::vector<double>> const table;
};
}  // namespace ComponentTransport
}  // namespace ProcessLib
