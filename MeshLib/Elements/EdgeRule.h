/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MathLib/Vector3.h"
#include "MeshLib/Node.h"

namespace MeshLib
{

class Element;

class EdgeRule
{
public:
    /// Constant: Dimension of this mesh element
    static const unsigned dimension = 1u;

    /// Constant: The number of faces
    static const unsigned n_faces = 0;

    /// Constant: The number of edges
    static const unsigned n_edges = 1;

    /// Returns the i-th face of the element.
    static const Element* getFace(const Element* /*e*/, unsigned /*i*/) { return nullptr; }

    /**
    * Checks if the node order of an element is correct by testing surface normals.
    * For 1D elements this always returns true.
    */
    static bool testElementNodeOrder(const Element* /*e*/) { return true; }

    /// Returns a normal vector perpendicular to the edge.
    /// inward_normal_vector the inward normal vector normal to a surface that
    /// the edge pertains to.
    static MathLib::Vector3 getNormalVector(
        const Element* e, MathLib::Vector3 inward_normal_vector)
    {
        Node* const* const _nodes = e->getNodes();
        MathLib::Vector3 directed_edge = {*_nodes[1], *_nodes[0]};
        return MathLib::crossProduct(inward_normal_vector, directed_edge);
    }
}; /* class */

}  // namespace MeshLib
