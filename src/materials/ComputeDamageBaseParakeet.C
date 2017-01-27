/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ComputeDamageBaseParakeet.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<ComputeDamageBaseParakeet>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredCoupledVar("damageInpVar", "order parameter for damage");
  params.addParam<std::string>("base_name", "Optional parameter that allows the user to define multiple mechanics material systems on the same block, i.e. for multiple phases");
  return params;
}

ComputeDamageBaseParakeet::ComputeDamageBaseParakeet(const InputParameters & parameters) :
    DerivativeMaterialInterface<Material>(parameters),
    _order_parameter(coupledValue("damageInpVar")),
    _grad_order_parameter(coupledGradient("damageInpVar")),
    _deta_dt(coupledDot("damageInpVar")),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "" ),
    _damage(declareProperty<Real>(_base_name + "damage")),
    _grad_damage(declareProperty<RealGradient>(_base_name + "grad_damage")),
    _ddamage_dt(declareProperty<Real>(_base_name + "ddamage_dt"))
{

  // fetch coupled variables and gradients (as stateful properties if necessary)
  //  _order_parameter = coupledValue("damageInpVar");
  //  if (_stateful_damage)
  //     _damage_old = coupledValueOld("damageInpVar");
  //  else
  //     _damage_old = _zero;

}

void
ComputeDamageBaseParakeet::initQpStatefulProperties()
{
    _damage[_qp] = 0.0;
    _grad_damage[_qp].zero();
    _ddamage_dt[_qp] = 0.0;
    //_state_var_psi_max[_qp] = 0.0;
    //_state_var_psi_max_old[_qp] = 0.0;
}

//void
//ComputeDamageBaseParakeet::computeQpProperties()
//{
//  computeQpDamageEvol();
//}
