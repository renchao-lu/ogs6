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
std::vector<std::size_t> searchList(std::vector<double> const& v,
                                    double const& find_for)
{
    std::vector<std::size_t> indexes;
    for (auto i = 0; i < v.size(); i++)
    {
        if (v[i] == find_for)
        {
            indexes.push_back(i);
        }
    }

    return indexes;
}

}  // namespace

namespace ProcessLib
{
namespace ComponentTransport
{
std::unique_ptr<LookupTable> createLookupTable(
    boost::optional<std::string> lookup_table_file,
    std::vector<std::pair<int, std::string>> const&
        process_id_to_component_name_map)
{
    if (!lookup_table_file)
        return nullptr;

    std::ifstream in(*lookup_table_file);
    if (!in)
    {
        OGS_FATAL("Could not open Kd matrix file '%s'.",
                  (*lookup_table_file).c_str());
    }

//    int num_skipped_lines = 5;
//    for (int i = 0; i < num_skipped_lines; ++i)
//    {
//        in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
//    }

    // input parameters
    std::string line;
//    std::getline(in, line);
//    std::vector<std::string> params;
//    boost::split(params, line, boost::is_any_of(" "));
//    assert(params.size() == 5);

//    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
//    std::getline(in, line);
//    std::istringstream iss(line);
//    double value, step_size;
//    int num_steps;
//    std::size_t num_tuples = 1;
//    std::vector<std::vector<double>> values;
//    while (iss >> value >> step_size >> num_steps)
//    {
//        std::vector<double> values_(num_steps);
//        value -= step_size;
//        std::generate(values_.begin(), values_.end(), [&step_size, &value] {
//            return value += step_size;
//        });
//        values.push_back(values_);
//        num_tuples *= num_steps;
//    }

//    std::map<std::string, std::vector<double>> input_parameters;
//    for (int i = 0; i < static_cast<int>(params.size()); ++i)
//    {
//        auto const& key = params[i];
//        input_parameters[key] = values[i];
//    }

    // radionuclides
    std::getline(in, line);
    std::vector<std::string> fields;
    boost::split(fields, line, boost::is_any_of("\t "));
    assert(fields.size() == 12);

    std::vector<std::string> concentration_fields;
    std::vector<std::string> surface_fields;
    std::vector<std::string> result_fields;
    for (auto const& field : fields)
    {
        if (field.find("_prev") != std::string::npos)
        {
            surface_fields.push_back(field);
        }
        else if (field.find("_new") != std::string::npos)
        {
            result_fields.push_back(field);
        }
        else
        {
            concentration_fields.push_back(field);
        }
    }

    // skip scaling factors by far
//    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::size_t num_items;
    in >> num_items;
    in.ignore();
    // read matrix
    std::map<std::string, std::vector<double>> table;
    for (auto item_id = 0; item_id < num_items; ++item_id)
    {
        std::getline(in, line);
        std::vector<std::string> cells;
        boost::trim_if(line, boost::is_any_of("\t "));
        boost::algorithm::split(cells, line, boost::is_any_of("\t "),
                                boost::token_compress_on);

        for (auto cell_id = 0; cell_id < cells.size(); ++cell_id)
        {
            table[fields[cell_id]].push_back(std::stod(cells[cell_id]));
        }
    }

    if (!in)
    {
        OGS_FATAL("Error when reading Kd matrix file '%s'",
                  (*lookup_table_file).c_str());
    }

    in.close();

    std::map<std::string, std::vector<double>> variables;
    for (auto const& concentration_field : concentration_fields)
    {
        variables.emplace(concentration_field, table[concentration_field]);
    }
    for (auto const& surface_field : surface_fields)
    {
        variables.emplace(surface_field, table[surface_field]);
    }

    std::map<std::string, int> concentration_field_to_process_id;
    for (auto const& concentration_field : concentration_fields)
    {
        auto pair = std::find_if(process_id_to_component_name_map.begin(),
                                 process_id_to_component_name_map.end(),
                                 [&concentration_field](auto const& p) {
                                     return p.second == concentration_field;
                                 });

        if (pair != process_id_to_component_name_map.end())
        {
            concentration_field_to_process_id[concentration_field] =
                pair->first;
        }
    }

    std::map<std::string, std::vector<double>> concentration_seeds;
    for (auto const& concentration_field : concentration_fields)
    {
        concentration_seeds[concentration_field] = table[concentration_field];
        BaseLib::makeVectorUnique(concentration_seeds[concentration_field]);
    }

    std::map<std::string, std::vector<double>> surface_field_seeds;
    for (auto const& surface_field : surface_fields)
    {
        surface_field_seeds[surface_field] = table[surface_field];
        BaseLib::makeVectorUnique(surface_field_seeds[surface_field]);
    }

    std::map<std::string, std::map<double, std::vector<std::size_t>>>
        radionuclides_concentrations;
    for (auto it = variables.begin(); it != variables.end(); it++)
    {
        BaseLib::makeVectorUnique(it->second);

        std::map<double, std::vector<std::size_t>> radionuclide_concentration;
        for (auto const& value : it->second)
        {
            auto matrix_index = searchList(table[it->first], value);

            radionuclide_concentration[value] = matrix_index;
        }

        radionuclides_concentrations[it->first] = radionuclide_concentration;
    }

    return std::make_unique<LookupTable>(
        std::move(concentration_field_to_process_id),
        std::move(concentration_seeds), std::move(surface_field_seeds),
        std::move(radionuclides_concentrations), std::move(table));
}

}  // namespace ComponentTransport
}  // namespace ProcessLib
