/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ComputeStressDamageBaseParakeet.h"
#include "Function.h"

template<>
InputParameters validParams<ComputeStressDamageBaseParakeet>()
{
  InputParameters params = validParams<Material>();
  params.addParam<std::string>("base_name", "Optional parameter that allows the user to define multiple mechanics material systems on the same block, i.e. for multiple phases");
  params.addParam<bool>("store_stress_old", false, "Parameter which indicates whether the old stress state, required for the HHT time integration scheme and Rayleigh damping, needs to be stored");
  return params;
}

ComputeStressDamageBaseParakeet::ComputeStressDamageBaseParakeet(const InputParameters & parameters) :
    DerivativeMaterialInterface<Material>(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "" ),
    _B_matrix(getMaterialPropertyByName<RankTwoTensor>(_base_name + "B_matrix")),
    _F_matrix(getMaterialPropertyByName<RankTwoTensor>(_base_name + "F_matrix")),
    _stress(declareProperty<RankTwoTensor>(_base_name + "stress")),
    _Jacobian_mult(declareProperty<RankFourTensor>(_base_name + "Jacobian_mult")),
    _store_stress_old(getParam<bool>("store_stress_old")),
    _damage(getMaterialPropertyByName<Real>(_base_name+"damage")),
    _grad_damage(getMaterialPropertyByName<RealGradient>(_base_name + "grad_damage")),
    _ddamage_dt(getMaterialPropertyByName<Real>(_base_name + "ddamage_dt")),
   _Rd_term1(declareProperty<Real>(_base_name + "Rd_term1")),
   _Rd_term2(declareProperty<RealGradient>(_base_name + "Rd_term2")),
   _Rd_term3(declareProperty<Real>(_base_name + "Rd_term3")),
   _Kdd_term1(declareProperty<Real>(_base_name + "Kdd_term1")),
   _Kdd_term2(declareProperty<Real>(_base_name + "Kdd_term2")),
   _Kdd_term3(declareProperty<Real>(_base_name + "Kdd_term3")),
   _Kdi_term1(declareProperty<RankTwoTensor>(_base_name + "Kdi_term1")),
   _Kid_term1(declareProperty<RankTwoTensor>(_base_name + "Kid_term1")),
   _state_var_psi_max(declareProperty<Real>(_base_name+"state_var_psi_max")),
   _state_var_psi_max_old(declarePropertyOld<Real>(_base_name+"state_var_psi_max")) 
{

  if (_store_stress_old)
  {
    declarePropertyOld<RankTwoTensor>(_base_name + "stress");
    declarePropertyOlder<RankTwoTensor>(_base_name + "stress");
  }

}

void
ComputeStressDamageBaseParakeet::initQpStatefulProperties()
{
  _stress[_qp].zero();
  _state_var_psi_max[_qp] = 0.0;
}

void
ComputeStressDamageBaseParakeet::computeQpProperties()
{
  computeQpStress();
}
