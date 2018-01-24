/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <iostream>

#include "CohesiveZoneModeI.h"
#include "LogPenalty.h"

#include "BaseLib/Error.h"
#include "MathLib/MathTools.h"

namespace MaterialLib
{
namespace Fracture
{
namespace CohesiveZoneModeI
{
namespace
{
struct MaterialPropertyValues
{
    double Kn = 0.0;
    double Ks = 0.0;
    double w_np = 0.0;
    double w_nf = 0.0;

    template <typename MaterialProperties>
    MaterialPropertyValues(MaterialProperties const& mp, double const t,
                           ProcessLib::SpatialPosition const& x,
                           double const aperture0)
    {
        Kn = mp.normal_stiffness(t, x)[0];
        Ks = mp.shear_stiffness(t, x)[0];
        w_np = mp.fracture_opening_at_peak_traction(t, x);
        w_nf = mp.fracture_opening_at_residual_traction(t, x);
    }
};

double computeDamage(double const damage_prev,
                     double const w_n,
                     double const w_np,
                     double const w_nf)
{
    return std::min(
        1.0,
        std::max(damage_prev, std::max(0.0, (w_n - w_np)) / (w_nf - w_np)));
}

}  // namespace

template <int DisplacementDim>
void CohesiveZoneModeI<DisplacementDim>::computeConstitutiveRelation(
    double const t,
    ProcessLib::SpatialPosition const& x,
    double const aperture0,
    Eigen::Ref<Eigen::VectorXd const>
        sigma0,
    Eigen::Ref<Eigen::VectorXd const>
    /*w_prev*/,
    Eigen::Ref<Eigen::VectorXd const>
        w,
    Eigen::Ref<Eigen::VectorXd const>
    /*sigma_prev*/,
    Eigen::Ref<Eigen::VectorXd>
        sigma,
    Eigen::Ref<Eigen::MatrixXd>
        C,
    typename FractureModelBase<DisplacementDim>::MaterialStateVariables&
        material_state_variables)
{
    assert(dynamic_cast<StateVariables<DisplacementDim> const*>(
               &material_state_variables) != nullptr);

    StateVariables<DisplacementDim> state =
        static_cast<StateVariables<DisplacementDim> const&>(
            material_state_variables);
    state.setInitialConditions();

    auto const mp = MaterialPropertyValues(_mp, t, x, aperture0);

    C.setZero();

    const int index_ns = DisplacementDim - 1;
    double const w_n = w[index_ns];
    for (int i = 0; i < index_ns; i++)
        C(i, i) = mp.Ks;

    sigma.noalias() = C * w;

    double const aperture = w_n + aperture0;

    state.damage = computeDamage(state.damage_prev, w_n, mp.w_np, mp.w_nf);

    C(index_ns, index_ns) = mp.Kn * mp.w_np * (1 - state.damage) / (mp.w_np + (mp.w_nf - mp.w_np) * state.damage);
    /***
    C(index_ns, index_ns) =
        mp.Kn * logPenaltyDerivative(aperture0, aperture, _penalty_aperture_cutoff);
    ***/

    sigma.coeffRef(index_ns) =
        C(index_ns, index_ns) * w_n ;//* logPenalty(aperture0, aperture, _penalty_aperture_cutoff);

    if (w_n < 0)
    {
        return;  /// Undamaged stiffness used in compression.
    }

    /*
    C = C * 1e-10;
    sigma = sigma * 0;
    return;

    */

    state.damage = computeDamage(state.damage_prev, w_n, mp.w_np, mp.w_nf);
    std::cerr << "w_nf / w_np " << mp.w_nf << "/" << mp.w_np << "\n";

    if (state.damage > state.damage_prev)
    {
        // If damage is increasing, provide extension to consistent tangent.

        Eigen::Matrix<double, DisplacementDim, 1> dd_dw =
            Eigen::Matrix<double, DisplacementDim, 1>::Zero();
        dd_dw[index_ns] = 1 / (mp.w_nf - mp.w_np);
        std::cerr << "dd_dw " << dd_dw << "\n";

        C(index_ns, index_ns) += mp.Kres;
        // TODO (naumov) noalias
        // C.noalias() = C * (1 - state.damage) - C * w * (dd_dw).transpose();
        C = C * (1 - state.damage) - C * w * (dd_dw).transpose();
    }
    else
    {
        // Degrade stiffness tensor in tension.
        C.noalias() = C * (1 - state.damage);
    }

    sigma.noalias() = sigma * (1 - state.damage);

    std::cerr << "damage " << state.damage << "\n";
    std::cerr << "exit sigma " << sigma << "\n";
    // TODO (nagel) Initial stress not considered, yet.
    // sigma.noalias() += sigma0;
}

template class CohesiveZoneModeI<2>;
template class CohesiveZoneModeI<3>;

}  // namespace CohesiveZoneModeI
}  // namespace Fracture
}  // namespace MaterialLib
