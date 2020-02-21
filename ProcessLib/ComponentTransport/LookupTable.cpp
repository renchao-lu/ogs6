/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LookupTable.h"
#include "BaseLib/Algorithm.h"

namespace ProcessLib
{
namespace ComponentTransport
{
std::pair<double, double> Field::getNeighboringPoints(double value)
{
    auto const it = std::lower_bound(points.begin(), points.end(), value);

    if (it == points.begin())
    {
        return std::make_pair(*it, *(it + 1));
    }
    else if (it == points.end())
    {
        return std::make_pair(*(points.rbegin() + 1), *points.rbegin());
    }
    else
    {
        return std::make_pair(*(it - 1), *it);
    }
}

void Table::lookup(
    std::vector<GlobalVector*> const& x,
    std::vector<std::unique_ptr<GlobalVector>> const& _x_previous_timestep)
{
    auto const n_nodes = x[0]->size();
    for (auto node_id = 0; node_id < n_nodes; ++node_id)
    {
        std::vector<InterpolationPoint> interpolation_points;
        {
            int i = 0;
            for (auto const& pair : concentration_field_to_process_id)
            {
                auto const process_id = pair.second;
                // check for negative concentration
                auto value = x[process_id]->get(node_id) < 0.
                                 ? 0.
                                 : x[process_id]->get(node_id);
                auto neighboring_points = fields[i].getNeighboringPoints(value);

                interpolation_points.emplace_back(value, neighboring_points);
                i++;
            }

            for (auto const& pair : concentration_field_to_process_id)
            {
                auto const process_id = pair.second;
                auto value =
                    _x_previous_timestep[process_id]->get(node_id) < 0.
                        ? 0.
                        : _x_previous_timestep[process_id]->get(node_id);
                auto neighboring_points = fields[i].getNeighboringPoints(value);
                interpolation_points.emplace_back(value, neighboring_points);
                i++;
            }
        }

        // look up the table
        std::vector<std::size_t> interpolation_point_indices;
        std::vector<double> point_set;
        for (auto const& interpolation_point : interpolation_points)
        {
            point_set.push_back(interpolation_point.neighboring_points.first);
        }
        interpolation_point_indices.push_back(
            getIntersectedIndex(fields, point_set));

        for (auto i = 0; i < static_cast<int>(interpolation_points.size()); ++i)
        {
            point_set[i] = interpolation_points[i].neighboring_points.second;
            interpolation_point_indices.push_back(
                getIntersectedIndex(fields, point_set));
            point_set[i] = interpolation_points[i].neighboring_points.first;
        }

        for (auto const& pair : concentration_field_to_process_id)
        {
            auto base_value =
                table.at(pair.first + "_new")[interpolation_point_indices[0]];
            auto new_value = base_value;

            // linear approximation
            for (auto i = 0; i < static_cast<int>(interpolation_points.size());
                 ++i)
            {
                auto& interpolation_point_value = table.at(
                    pair.first + "_new")[interpolation_point_indices[i + 1]];
                auto slope =
                    (interpolation_point_value - base_value) /
                    (interpolation_points[i].neighboring_points.second -
                     interpolation_points[i].neighboring_points.first);

                new_value +=
                    slope * (interpolation_points[i].value -
                             interpolation_points[i].neighboring_points.first);
            }

            x[pair.second]->set(node_id, new_value);
        }
    }
}

std::size_t Table::getIntersectedIndex(std::vector<Field> const& fields,
                                       std::vector<double> values)
{
    std::vector<std::size_t> vec1;
    std::vector<std::size_t> temp_vec;
    std::vector<std::size_t> vec =
        fields[0].indices[BaseLib::findIndex(fields[0].points, values[0])];

    auto num_variables = values.size();
    for (auto i = 0; i < static_cast<int>(num_variables); ++i)
    {
        auto index = BaseLib::findIndex(fields[i].points, values[i]);
        vec1 = fields[i].indices[index];
        // intersect vectors
        std::set_intersection(vec1.begin(), vec1.end(), vec.begin(), vec.end(),
                              std::back_inserter(temp_vec));

        std::swap(vec, temp_vec);
        temp_vec.clear();
    }

    return vec[0];
}
}  // namespace ComponentTransport
}  // namespace ProcessLib
