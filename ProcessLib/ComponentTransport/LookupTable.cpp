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
std::pair<double, double> Field::getNeighboringDataPoints(double value)
{
    auto const it =
        std::lower_bound(data_points.begin(), data_points.end(), value);

    if (it == data_points.begin())
    {
        return std::make_pair(*it, *(it + 1));
    }
    else if (it == data_points.end())
    {
        return std::make_pair(*(data_points.rbegin() + 1),
                              *data_points.rbegin());
    }
    else
    {
        return std::make_pair(*(it - 1), *it);
    }
}

void LookupTable::lookup(
    std::vector<GlobalVector*> const& x,
    std::vector<std::unique_ptr<GlobalVector>> const& _x_previous_timestep)
{
    auto const n_nodes = x[0]->size();
    for (auto node_id = 0; node_id < n_nodes; ++node_id)
    {
        std::vector<double> nodal_x;
        for (auto const& pair : concentration_field_to_process_id)
        {
            // check for negative concentration
            double value = x[pair.second]->get(node_id) < 0.
                               ? 0.
                               : x[pair.second]->get(node_id);
            nodal_x.push_back(value);
        }
        for (auto const& pair : concentration_field_to_process_id)
        {
            double value =
                _x_previous_timestep[pair.second]->get(node_id) < 0.
                    ? 0.
                    : _x_previous_timestep[pair.second]->get(node_id);
            nodal_x.push_back(value);
        }

        std::vector<std::pair<double, double>> neighboring_data_points;
        for (auto i = 0; i < static_cast<int>(nodal_x.size()); ++i)
        {
            auto ip = fields[i].getNeighboringDataPoints(nodal_x[i]);
            neighboring_data_points.push_back(ip);
        }

        // look up the table
        std::vector<double> reference_points;
        for (auto const& neighboring_data_point : neighboring_data_points)
        {
            reference_points.push_back(neighboring_data_point.first);
        }

        auto const base_point_index =
            getIntersectedIndex(fields, reference_points);

        for (auto const& pair : concentration_field_to_process_id)
        {
            auto base_value = table.at(pair.first + "_new")[base_point_index];
            auto new_value = base_value;

            // linear approximation
            for (auto i = 0; i < static_cast<int>(fields.size()); ++i)
            {
                reference_points[i] = neighboring_data_points[i].second;
                auto const interpolation_point_index =
                    getIntersectedIndex(fields, reference_points);
                auto& interpolation_point_value =
                    table.at(pair.first + "_new")[interpolation_point_index];
                auto slope = (interpolation_point_value - base_value) /
                             (neighboring_data_points[i].second -
                              neighboring_data_points[i].first);

                new_value +=
                    slope * (nodal_x[i] - neighboring_data_points[i].first);
                reference_points[i] = neighboring_data_points[i].first;
            }

            x[pair.second]->set(node_id, new_value);
        }
    }
}

std::size_t LookupTable::getIntersectedIndex(std::vector<Field> const& fields,
                                             std::vector<double> values)
{
    std::vector<std::size_t> vec1;
    std::vector<std::size_t> temp_vec;
    std::vector<std::size_t> vec =
        fields[0].indices[BaseLib::findIndex(fields[0].data_points, values[0])];

    auto num_variables = values.size();
    for (auto i = 0; i < static_cast<int>(num_variables); ++i)
    {
        auto index = BaseLib::findIndex(fields[i].data_points, values[i]);
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
