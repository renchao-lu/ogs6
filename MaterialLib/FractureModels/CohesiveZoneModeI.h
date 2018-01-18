/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>
#include <utility>

#include "ProcessLib/Parameter/Parameter.h"

#include "FractureModelBase.h"

namespace MaterialLib
{
namespace Fracture
{
namespace CohesiveZoneModeI
{
template <int DisplacementDim>
class CohesiveZoneModeI final : public FractureModelBase<DisplacementDim>
{
public:
    /// Variables specific to the material model
    struct MaterialProperties
    {
        using P = ProcessLib::Parameter<double>;
        using X = ProcessLib::SpatialPosition;

    public:
        MaterialProperties(P const& normal_stiffness_,
                           P const& shear_stiffness_,
                           P const& fracture_toughness,
                           P const& peak_normal_traction,
                           double const residual_stiffness)
            : normal_stiffness(normal_stiffness_),
              shear_stiffness(shear_stiffness_),
              _fracture_toughness(fracture_toughness),
              _peak_normal_traction(peak_normal_traction),
              _residual_stiffness(residual_stiffness)
        {
        }

        /// Assuming initially stress-free state.
        double fracture_opening_at_peak_traction(double const t,
                                                 X const& x) const
        {
            return _peak_normal_traction(t, x)[0] / normal_stiffness(t, x)[0];
        }

        /// Assuming initially stress-free state.
        double fracture_opening_at_residual_traction(double const t,
                                                     X const& x) const
        {
            return 2 * _fracture_toughness(t, x)[0] /
                   _peak_normal_traction(t, x)[0];
        }

    public:
        /// Normal stiffness given in units of stress.
        P const& normal_stiffness;
        /// Shear stiffness given in units of stress.
        P const& shear_stiffness;

    private:
        /// Fracture toughness/critical energy release rate given in of stress
        /// times lengths.
        P const& _fracture_toughness;
        /// Peak normal traction given in units of stress.
        P const& _peak_normal_traction;
        /// Residual stiffness given in units of stress.
        double const _residual_stiffness;
    };

    struct MaterialStateVariables
        : public FractureModelBase<DisplacementDim>::MaterialStateVariables
    {
        void pushBackState() override {}
    };

    std::unique_ptr<
        typename FractureModelBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() override
    {
        return std::unique_ptr<typename FractureModelBase<
            DisplacementDim>::MaterialStateVariables>{
            new MaterialStateVariables};
    }

public:
    explicit CohesiveZoneModeI(double const penalty_aperture_cutoff,
                               bool const tension_cutoff,
                               MaterialProperties material_properties)
        : _penalty_aperture_cutoff(penalty_aperture_cutoff),
          _tension_cutoff(tension_cutoff),
          _mp(std::move(material_properties))
    {
    }

    /**
     * Computation of the constitutive relation for the Mohr-Coulomb model.
     *
     * @param t           current time
     * @param x           current position in space
     * @param aperture0   initial fracture's aperture
     * @param sigma0      initial stress
     * @param w_prev      fracture displacement at previous time step
     * @param w           fracture displacement at current time step
     * @param sigma_prev  stress at previous time step
     * @param sigma       stress at current time step
     * @param Kep         tangent matrix for stress and fracture displacements
     * @param material_state_variables   material state variables
     */
    void computeConstitutiveRelation(
        double const t,
        ProcessLib::SpatialPosition const& x,
        double const aperture0,
        Eigen::Ref<Eigen::VectorXd const>
            sigma0,
        Eigen::Ref<Eigen::VectorXd const>
            w_prev,
        Eigen::Ref<Eigen::VectorXd const>
            w,
        Eigen::Ref<Eigen::VectorXd const>
            sigma_prev,
        Eigen::Ref<Eigen::VectorXd>
            sigma,
        Eigen::Ref<Eigen::MatrixXd>
            Kep,
        typename FractureModelBase<DisplacementDim>::MaterialStateVariables&
            material_state_variables) override;

private:
    /// Compressive normal displacements above this value will not enter the
    /// computation of the normal stiffness modulus of the fracture.
    /// \note Setting this to the initial aperture value allows negative
    /// apertures.
    double const _penalty_aperture_cutoff;

    /// If set no resistance to open the fracture over the initial aperture is
    /// opposed.
    bool const _tension_cutoff;

    MaterialProperties _mp;
};

}  // namespace CohesiveZoneModeI
}  // namespace Fracture
}  // namespace MaterialLib

namespace MaterialLib
{
namespace Fracture
{
namespace CohesiveZoneModeI
{
extern template class CohesiveZoneModeI<2>;
extern template class CohesiveZoneModeI<3>;
}  // namespace CohesiveZoneModeI
}  // namespace Fracture
}  // namespace MaterialLib
