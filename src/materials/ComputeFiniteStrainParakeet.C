/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ComputeFiniteStrainParakeet.h"
#include "Assembly.h"

// libmesh includes
#include "libmesh/quadrature.h"
// this is for JxW
//#include "libmesh/fe_map.h"
//#include "libmesh/fe.h"

template<>
InputParameters validParams<ComputeFiniteStrainParakeet>()
{
  InputParameters params = validParams<ComputeStrainBaseParakeet>();
  params.addClassDescription("Compute a strain measure for finite strains.");
  params.set<bool>("stateful_deformation_gradient") = true;

  return params;
}

ComputeFiniteStrainParakeet::ComputeFiniteStrainParakeet(const InputParameters & parameters) :
    ComputeStrainBaseParakeet(parameters)
{
}

void
ComputeFiniteStrainParakeet::computeProperties()
{

  for (_qp=0;_qp<_qrule->n_points();++_qp)
  {
    // Deformation gradient
    RankTwoTensor A((*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp], (*_grad_disp[2])[_qp]); //Deformation gradient
    RankTwoTensor Fbar((*_grad_disp_old[0])[_qp], (*_grad_disp_old[1])[_qp], (*_grad_disp_old[2])[_qp]); //Old Deformation gradient

    _deformation_gradient[_qp] = A;
    _deformation_gradient[_qp].addIa(1.0);//Gauss point deformation gradient

  }

  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    //computeQpStrain();
     _F_matrix[_qp] = _deformation_gradient[_qp];
     _B_matrix[_qp] = _deformation_gradient[_qp]*_deformation_gradient[_qp].transpose();
  }
}

