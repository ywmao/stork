/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ComputeDamageParakeet.h"
#include "Assembly.h"

// libmesh includes
#include "libmesh/quadrature.h"
// this is for JxW
//#include "libmesh/fe_map.h"
//#include "libmesh/fe.h"

template<>
InputParameters validParams<ComputeDamageParakeet>()
{
  InputParameters params = validParams<ComputeDamageBaseParakeet>();
  params.addClassDescription("Compute a damage for order parameter.");

  return params;
}

ComputeDamageParakeet::ComputeDamageParakeet(const InputParameters & parameters) :
    ComputeDamageBaseParakeet(parameters)
{
}

void
ComputeDamageParakeet::computeProperties()
{

  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
     _damage[_qp] = _order_parameter[_qp];
     _grad_damage[_qp] = _grad_order_parameter[_qp];
     _ddamage_dt[_qp] = _deta_dt[_qp];
  }

}

