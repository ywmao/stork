/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "StressDivergenceDamageParakeet.h"
#include "Material.h"
#include "MooseMesh.h"
#include "ElasticityTensorToolsParakeet.h"
#include "libmesh/quadrature.h"

template<>
InputParameters validParams<StressDivergenceDamageParakeet>()
{
  InputParameters params = validParams<ALEKernelParakeet>();
  params.addClassDescription("Stress divergence kernel for the Cartesian coordinate system");
  params.addRequiredParam<unsigned int>("component", "An integer corresponding to the direction the variable this kernel acts in. (0 for x, 1 for y, 2 for z)");
  // shenghua thought this one shall be commented.
  params.addRequiredCoupledVar("displacements", "The string of displacements suitable for the problem statement");
  params.addRequiredCoupledVar("damageInpVar", "The damage quantity.");
  params.addParam<std::string>("base_name", "Material property base name");
  params.set<bool>("use_displaced_mesh") = false;
  params.addParam<bool>("use_finite_deform_jacobian", false, "Jacobian for corotational finite strain");
  return params;
}


StressDivergenceDamageParakeet::StressDivergenceDamageParakeet(const InputParameters & parameters) :
    ALEKernelParakeet(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _use_finite_deform_jacobian(getParam<bool>("use_finite_deform_jacobian")),
    _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
    _Jacobian_mult(getMaterialPropertyByName<RankFourTensor>(_base_name + "Jacobian_mult")),
    _Kid_term1(getMaterialPropertyByName<RankTwoTensor>(_base_name + "Kid_term1")),
    _component(getParam<unsigned int>("component")),
    _ndisp(coupledComponents("displacements")),
    _disp_var(3),
    _damage_var(coupled("damageInpVar"))
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);

  //_damage_var = coupled("damageInpVar");
 
  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_ndisp != _mesh.dimension())
    mooseError("The number of displacement variables supplied must match the mesh dimension.");
}

Real
StressDivergenceDamageParakeet::computeQpResidual()
{
  return _stress[_qp].row(_component) * _grad_test[_i][_qp];
}

void
StressDivergenceDamageParakeet::computeJacobian()
{
    Kernel::computeJacobian();
}

void
StressDivergenceDamageParakeet::computeOffDiagJacobian(unsigned int jvar)
{
    Kernel::computeOffDiagJacobian(jvar);
}

Real
StressDivergenceDamageParakeet::computeQpJacobian()
{
  return ElasticityTensorToolsParakeet::elasticJacobian(_Jacobian_mult[_qp], _component, _component, _grad_test[_i][_qp], _grad_phi[_j][_qp]);
}

Real
StressDivergenceDamageParakeet::computeQpOffDiagJacobian(unsigned int jvar)
{

   unsigned int coupled_component = 0;
   bool active(false);
   for(unsigned int i=0;i< _ndisp; ++i)
      if(jvar == _disp_var[i])
      {  
         coupled_component = i;
         active = true;
      }
   if(active)
   {        
          return ElasticityTensorToolsParakeet::elasticJacobian(_Jacobian_mult[_qp], _component, coupled_component,
                                          _grad_test[_i][_qp], _grad_phi[_j][_qp]);
   }

   if(jvar == _damage_var)
   {
      Real AA = 0;
      for(unsigned k=0;k<3;k++)
         AA = AA+ _Kid_term1[_qp](_component,k)*_phi[_i][_qp]*_grad_test[_j][_qp](k);
      return AA;
   }
   return 0;
}
