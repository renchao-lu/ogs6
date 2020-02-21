/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <boost/algorithm/string.hpp>
#include <boost/optional.hpp>
#include <fstream>
#include <map>

#include "BaseLib/Error.h"
#include "BaseLib/Algorithm.h"
#include "CreateLookupTable.h"
#include "LookupTable.h"

namespace
{
std::size_t getSingleValueFromFile(std::ifstream& in)
{
    std::size_t v;
    in >> v;
    in.ignore();

    return v;
}

std::vector<std::size_t> searchList(std::vector<double> const& v,
                                    double const& find_for)
{
    std::vector<std::size_t> indices;
    for (auto i = 0; i < static_cast<int>(v.size()); i++)
    {
        if (v[i] == find_for)
        {
            indices.push_back(i);
        }
    }

    return indices;
}

}  // namespace

namespace ProcessLib
{
namespace ComponentTransport
{
std::unique_ptr<Table> createTable(
    boost::optional<std::string> spreadsheet_file,
    std::vector<std::pair<int, std::string>> const&
        process_id_to_component_name_map)
{
    if (!spreadsheet_file)
        return nullptr;

    std::ifstream in(*spreadsheet_file);
    if (!in)
    {
        OGS_FATAL("Could not open spreadsheet file '%s'.",
                  (*spreadsheet_file).c_str());
    }

    auto const num_fields = getSingleValueFromFile(in);
    auto const num_items = getSingleValueFromFile(in);

    // read field names
    std::string line;
    std::getline(in, line);
    std::vector<std::string> field_names;
    boost::split(field_names, line, boost::is_any_of("\t "));
    assert(field_names.size() == num_fields);

    std::vector<std::string> variable_fields;
    std::vector<std::string> result_fields;
    for (auto const& field : field_names)
    {
        if (field.find("_new") != std::string::npos)
        {
            result_fields.push_back(field);
        }
        else
        {
            variable_fields.push_back(field);
        }
    }

    // read table
    std::map<std::string, std::vector<double>> table;
    for (auto item_id = 0; item_id < static_cast<int>(num_items); ++item_id)
    {
        std::getline(in, line);
        std::vector<std::string> buckets;
        boost::trim_if(line, boost::is_any_of("\t "));
        boost::algorithm::split(buckets, line, boost::is_any_of("\t "),
                                boost::token_compress_on);

        for (auto bucket_id = 0; bucket_id < static_cast<int>(buckets.size());
             ++bucket_id)
        {
            table[field_names[bucket_id]].push_back(
                std::stod(buckets[bucket_id]));
        }
    }

    if (!in)
    {
        OGS_FATAL("Error in reading spreadsheet file '%s'",
                  (*spreadsheet_file).c_str());
    }

    in.close();

    std::vector<std::pair<std::string, int>> concentration_field_to_process_id;
    for (auto const& variable_field : variable_fields)
    {
        auto pair = std::find_if(process_id_to_component_name_map.begin(),
                                 process_id_to_component_name_map.end(),
                                 [&variable_field](auto const& p) {
                                     return p.second == variable_field;
                                 });

        if (pair != process_id_to_component_name_map.end())
        {
            concentration_field_to_process_id.emplace_back(
                std::make_pair(variable_field, pair->first));
        }
    }

    std::vector<Field> fields;
    for (auto const variable_field : variable_fields)
    {
        auto interpolation_points = table[variable_field];
        BaseLib::makeVectorUnique(interpolation_points);

        std::vector<std::vector<std::size_t>> indices_vec;
        for (auto const& ip : interpolation_points)
        {
            auto indices = searchList(table[variable_field], ip);
            indices_vec.push_back(indices);
        }

        fields.emplace_back(variable_field, interpolation_points, indices_vec);
    }

    return std::make_unique<Table>(std::move(concentration_field_to_process_id),
                                   std::move(fields), std::move(table));
}

}  // namespace ComponentTransport
}  // namespace ProcessLib
