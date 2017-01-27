/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ComputeArrudaBoyceDamageMPParakeet.h"

template<>
InputParameters validParams<ComputeArrudaBoyceDamageMPParakeet>()
{
  InputParameters params = validParams<ComputeArrudaBoyceDamageMPBaseParakeet>();
  params.addClassDescription("Compute the material properties for AB model with damage.");
  params.addParam<Real>("G0", "Initial shear modulus for the material.");
  params.addParam<Real>("lmd_L", "Locking stretch for the material.");
  params.addParam<Real>("kbulk", "Bulk modulus for the material.");
  params.addParam<Real>("PsiCrit","Critical energy to break.");
  params.addParam<Real>("ell","Length scale in fracture.");
  params.addParam<Real>("BetaCoeff","Time derivative coeff in fracture process.");
  return params;
}

ComputeArrudaBoyceDamageMPParakeet::ComputeArrudaBoyceDamageMPParakeet(const InputParameters & parameters) :
    ComputeArrudaBoyceDamageMPBaseParakeet(parameters),
    _G0R_set( parameters.isParamValid("G0") ),
    _lmd_LR_set( parameters.isParamValid("lmd_L") ),
    _kbulkR_set( parameters.isParamValid("kbulk") ),
    _PsiCritR_set( parameters.isParamValid("PsiCrit") ),
    _ellR_set( parameters.isParamValid("ell") ),
    _BetaCoeffR_set( parameters.isParamValid("BetaCoeff") ),
    _G0R( _G0R_set ? getParam<Real>("G0") : -1 ),
    _lmd_LR( _lmd_LR_set ? getParam<Real>("lmd_L") : -1 ),
    _kbulkR( _kbulkR_set ?  getParam<Real>("kbulk") : -1 ),
    _PsiCritR( _PsiCritR_set ?  getParam<Real>("PsiCrit") : -1 ),
    _ellR( _ellR_set ?  getParam<Real>("ell") : -1 ),
    _BetaCoeffR( _BetaCoeffR_set ?  getParam<Real>("BetaCoeff") : -1 )
{
}

void
ComputeArrudaBoyceDamageMPParakeet::computeQpMaterialProperties()
{
  //Assign elasticity tensor at a given quad point
  _G0[_qp] = _G0R;
  _lmd_L[_qp] = _lmd_LR;
  _kbulk[_qp] = _kbulkR;
  _PsiCrit[_qp] = _PsiCritR;
  _ell[_qp] = _ellR;
  _BetaCoeff[_qp] = _BetaCoeffR;
}

