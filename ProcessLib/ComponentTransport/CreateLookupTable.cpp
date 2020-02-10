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
#include "CreateLookupTable.h"
#include "LookupTable.h"

namespace ProcessLib
{
namespace ComponentTransport
{
std::unique_ptr<LookupTable> createLookupTable(
    boost::optional<std::string> lookup_table_file)
{
    if (!lookup_table_file)
        return nullptr;

    std::ifstream in(*lookup_table_file);
    if (!in)
    {
        OGS_FATAL("Could not open Kd matrix file '%s'.",
                  (*lookup_table_file).c_str());
    }

    int num_skipped_lines = 5;
    for (int i = 0; i < num_skipped_lines; ++i)
    {
        in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    // input parameters
    std::string line;
    std::getline(in, line);
    std::vector<std::string> params;
    boost::split(params, line, boost::is_any_of("\t "));
    assert(params.size() == 5);

    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::getline(in, line);
    std::istringstream iss(line);
    double value, step_size;
    int num_steps;
    std::size_t num_tuples = 1;
    std::vector<std::vector<double>> values;
    while (iss >> value >> step_size >> num_steps)
    {
        std::vector<double> values_(num_steps);
        value -= step_size;
        std::generate(values_.begin(), values_.end(), [&step_size, &value] {
            return value += step_size;
        });
        values.push_back(values_);
        num_tuples *= num_steps;
    }

    std::map<std::string, std::vector<double>> input_parameters;
    for (int i = 0; i < static_cast<int>(params.size()); ++i)
    {
        auto const& key = params[i];
        input_parameters[key] = values[i];
    }

    // radionuclides
    std::getline(in, line);
    std::vector<std::string> radionuclides;
    boost::split(radionuclides, line, boost::is_any_of("\t "));
    assert(radionuclides.size() == 8);

    // skip scaling factors by far
    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // read matrix
    std::map<std::string, std::vector<double>> kd_matrix;
    for (std::size_t tuple_id = 0; tuple_id < num_tuples; ++tuple_id)
    {
        std::getline(in, line);
        std::vector<std::string> items;
        boost::trim_if(line, boost::is_any_of("\t "));
        boost::algorithm::split(items, line, boost::is_any_of("\t "),
                                boost::token_compress_on);

        for (int item_id = 0; item_id < static_cast<int>(items.size());
             ++item_id)
        {
            kd_matrix[radionuclides[item_id]].push_back(
                std::stod(items[item_id]));
        }
    }

    if (!in)
    {
        OGS_FATAL("Error when reading Kd matrix file '%s'",
                  (*lookup_table_file).c_str());
    }

    in.close();

    return std::make_unique<LookupTable>(input_parameters, kd_matrix);
}

}  // namespace ComponentTransport
}  // namespace ProcessLib
