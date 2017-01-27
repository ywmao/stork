/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ComputeStrainBaseParakeet.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<ComputeStrainBaseParakeet>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredCoupledVar("displacements", "The displacements appropriate for the simulation geometry and coordinate system");
  params.addParam<std::string>("base_name", "Optional parameter that allows the user to define multiple mechanics material systems on the same block, i.e. for multiple phases");
  return params;
}

ComputeStrainBaseParakeet::ComputeStrainBaseParakeet(const InputParameters & parameters) :
    DerivativeMaterialInterface<Material>(parameters),
    _ndisp(coupledComponents("displacements")),
    _disp(3),
    _grad_disp(3),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "" ),
    _B_matrix(declareProperty<RankTwoTensor>(_base_name + "B_matrix")),
    _F_matrix(declareProperty<RankTwoTensor>(_base_name + "F_matrix")),
    _stateful_displacements(_fe_problem.isTransient()),
    _stateful_deformation_gradient(getParam<bool>("stateful_deformation_gradient") && _fe_problem.isTransient()),
    _grad_disp_old(3),
    _deformation_gradient(declareProperty<RankTwoTensor>(_base_name + "deformation_gradient")),
    _deformation_gradient_old(_stateful_deformation_gradient ? &declarePropertyOld<RankTwoTensor>(_base_name + "deformation_gradient") : NULL)
{


  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_ndisp != _mesh.dimension())
    mooseError("The number of variables supplied in 'displacements' must match the mesh dimension.");

  // fetch coupled variables and gradients (as stateful properties if necessary)
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    _disp[i] = &coupledValue("displacements", i);
    _grad_disp[i] = &coupledGradient("displacements", i);
    if (_stateful_displacements)
       _grad_disp_old[i] = &coupledGradientOld("displacements" ,i);
    else
       _grad_disp_old[i] = &_grad_zero;
  }

  // set unused dimensions to zero
  for (unsigned i = _ndisp; i < 3; ++i)
  {
    _disp[i] = &_zero;
    _grad_disp[i] = &_grad_zero;
    _grad_disp_old[i] = &_grad_zero;
  }
}

void
ComputeStrainBaseParakeet::initQpStatefulProperties()
{
    _deformation_gradient[_qp].zero();
    _deformation_gradient[_qp].addIa(1.0);
    _F_matrix[_qp].zero();
    _F_matrix[_qp].addIa(1.0);
    _B_matrix[_qp].zero();
    _B_matrix[_qp].addIa(1.0);

}
