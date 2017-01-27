/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "DamageParakeet.h"
#include "Material.h"
#include "MooseMesh.h"
#include "ElasticityTensorToolsParakeet.h"
#include "libmesh/quadrature.h"

template<>
InputParameters validParams<DamageParakeet>()
{
  InputParameters params = validParams<ALEKernelParakeet>();
  params.addClassDescription("damage kernel for the Cartesian coordinate system");
  //params.addRequiredCoupledVar("damageInpVar", "The string of damage suitable for the problem statement");
  params.addRequiredCoupledVar("displacements", "The coupled displacement DOF.");
  params.addParam<std::string>("base_name", "Material property base name");
  return params;
}


DamageParakeet::DamageParakeet(const InputParameters & parameters) :
    ALEKernelParakeet(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _Rd_term1(getMaterialPropertyByName<Real>(_base_name + "Rd_term1")),
    _Rd_term2(getMaterialPropertyByName<RealGradient>(_base_name + "Rd_term2")),
    _Rd_term3(getMaterialPropertyByName<Real>(_base_name + "Rd_term3")),
    _Kdd_term1(getMaterialPropertyByName<Real>(_base_name + "Kdd_term1")),
    _Kdd_term2(getMaterialPropertyByName<Real>(_base_name + "Kdd_term3")),
    _Kdd_term3(getMaterialPropertyByName<Real>(_base_name + "Kdd_term2")),
    _Kdi_term1(getMaterialPropertyByName<RankTwoTensor>(_base_name + "Kdi_term1")),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp)
{

    // _ndisp & _disp_var are used for the offJacobian calculation.
    for(unsigned int i=0;i<_ndisp;++i)
       _disp_var[i] = coupled("displacements",i);
}

Real
DamageParakeet::computeQpResidual()
{
  //return _stress[_qp].row(_component) * _grad_test[_i][_qp];
  // more details about this residual see example 3.
  //     http://mooseframework.org/wiki/MooseExamples/Example_03/
  //return  -(_BetaCoeff[_qp]*_ddamage_dt[_qp]+_dpsi_deta[_qp])*_test[_i][_qp]-_xip[_qp]*_grad_test[_i][_qp];
  //return (A1+B1);
  //return 0;
  return _Rd_term1[_qp]*_test[_i][_qp]+_Rd_term2[_qp]*_grad_test[_i][_qp]+_Rd_term3[_qp]*_test[_i][_qp]; 
}

void
DamageParakeet::computeJacobian()
{
    Kernel::computeJacobian();
}

void
DamageParakeet::computeOffDiagJacobian(unsigned int jvar)
{
    Kernel::computeOffDiagJacobian(jvar);
}

Real
DamageParakeet::computeQpJacobian()
{
    return _test[_i][_qp]*_Kdd_term1[_qp]*_phi[_j][_qp]+_grad_test[_i][_qp]*_Kdd_term2[_qp]*_grad_phi[_j][_qp]+_test[_i][_qp]*_Kdd_term3[_qp]*_phi[_j][_qp];
}

Real
DamageParakeet::computeQpOffDiagJacobian(unsigned int jvar)
{
    for(unsigned int i=0;i<_ndisp;++i)
       if(jvar == _disp_var[i])
       {
              RealGradient AHL;
              AHL.zero();
              for(unsigned beta=0;beta<3;beta++)
                   AHL(i) = AHL(i)+_Kdi_term1[_qp](i,beta)*_grad_phi[_j][_qp](beta);
           return _test[_qp][_i]*AHL(i);
       }
    return 0;
}
