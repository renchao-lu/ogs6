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

namespace ProcessLib
{
namespace ComponentTransport
{
struct InterpolationPoint
{
    InterpolationPoint(double value_,
                       std::pair<double, double> neighboring_points_)
        : value(value_), neighboring_points(neighboring_points_)
    {
    }

    double value;
    std::pair<double, double> neighboring_points;
};

struct Field
{
    Field(std::string name_, std::vector<double> points_,
          std::vector<std::vector<std::size_t>> indices_)
        : name(name_), points(points_), indices(indices_)
    {
    }

    std::pair<double, double> getNeighboringPoints(double value);

    std::string name;
    std::vector<double> points;
    std::vector<std::vector<std::size_t>> indices;
};

struct Table
{
    Table(std::vector<std::pair<std::string, int>>&&
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
        std::vector<std::unique_ptr<GlobalVector>> const& _x_previous_timestep);

    std::size_t getIntersectedIndex(
        std::vector<Field> const& fields, std::vector<double> values);

    std::vector<std::pair<std::string, int>> concentration_field_to_process_id;
    std::vector<Field> fields;
    std::map<std::string, std::vector<double>> const table;
};
}  // namespace ComponentTransport
}  // namespace ProcessLib
