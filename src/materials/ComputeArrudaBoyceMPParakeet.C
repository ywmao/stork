/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ComputeArrudaBoyceMPParakeet.h"

template<>
InputParameters validParams<ComputeArrudaBoyceMPParakeet>()
{
  InputParameters params = validParams<ComputeArrudaBoyceMPBaseParakeet>();
  params.addClassDescription("Compute the material properties for AB model.");
  params.addParam<Real>("G0", "Initial shear modulus for the material.");
  params.addParam<Real>("lmd_L", "Locking stretch for the material.");
  params.addParam<Real>("kbulk", "Bulk modulus for the material.");
  return params;
}

ComputeArrudaBoyceMPParakeet::ComputeArrudaBoyceMPParakeet(const InputParameters & parameters) :
    ComputeArrudaBoyceMPBaseParakeet(parameters),
    _G0R_set( parameters.isParamValid("G0") ),
    _lmd_LR_set( parameters.isParamValid("lmd_L") ),
    _kbulkR_set( parameters.isParamValid("kbulk") ),
    _G0R( _G0R_set ? getParam<Real>("G0") : -1 ),
    _lmd_LR( _lmd_LR_set ? getParam<Real>("lmd_L") : -1 ),
    _kbulkR( _kbulkR_set ?  getParam<Real>("kbulk") : -1 )
{
}

void
ComputeArrudaBoyceMPParakeet::computeQpMaterialProperties()
{
  //Assign elasticity tensor at a given quad point
  _G0[_qp] = _G0R;
  _lmd_L[_qp] = _lmd_LR;
  _kbulk[_qp] = _kbulkR;
}

