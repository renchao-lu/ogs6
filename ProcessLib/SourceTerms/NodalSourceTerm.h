/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/Parameter/Parameter.h"

namespace ProcessLib
{
class NodalSourceTerm final
{
public:
    NodalSourceTerm(const NumLib::LocalToGlobalIndexMap& dof_table,
                    std::size_t const mesh_id, std::size_t const node_id,
                    const int variable_id, const int component_id,
                    Parameter<double> const& parameter);

    void integrateNodalSourceTerm(
        const double t,
        GlobalVector& b) const;

private:
    NumLib::LocalToGlobalIndexMap const& _dof_table;
    std::size_t const _mesh_id;
    std::size_t const _node_id;
    int const _variable_id;
    int const _component_id;
    Parameter<double> const& _parameter;
};

}   // namespace ProcessLib
