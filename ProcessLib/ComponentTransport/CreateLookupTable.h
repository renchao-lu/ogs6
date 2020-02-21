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

#include <boost/optional/optional_fwd.hpp>
#include <memory>
#include <string>
#include <vector>

namespace ProcessLib
{
namespace ComponentTransport
{
struct Table;

std::unique_ptr<Table> createTable(
    boost::optional<std::string> lookup_table_file,
    std::vector<std::pair<int, std::string>> const&
        process_id_to_component_name_map);

}  // namespace ComponentTransport
}  // namespace ProcessLib
