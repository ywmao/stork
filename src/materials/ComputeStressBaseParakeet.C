/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ComputeStressBaseParakeet.h"
#include "Function.h"

template<>
InputParameters validParams<ComputeStressBaseParakeet>()
{
  InputParameters params = validParams<Material>();
  params.addParam<std::string>("base_name", "Optional parameter that allows the user to define multiple mechanics material systems on the same block, i.e. for multiple phases");
  params.addParam<bool>("store_stress_old", false, "Parameter which indicates whether the old stress state, required for the HHT time integration scheme and Rayleigh damping, needs to be stored");
  return params;
}

ComputeStressBaseParakeet::ComputeStressBaseParakeet(const InputParameters & parameters) :
    DerivativeMaterialInterface<Material>(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "" ),
    _B_matrix(getMaterialPropertyByName<RankTwoTensor>(_base_name + "B_matrix")),
    _F_matrix(getMaterialPropertyByName<RankTwoTensor>(_base_name + "F_matrix")),
    _stress(declareProperty<RankTwoTensor>(_base_name + "stress")),
    _Jacobian_mult(declareProperty<RankFourTensor>(_base_name + "Jacobian_mult")),
    _store_stress_old(getParam<bool>("store_stress_old"))
{
  //Declares old stress and older stress if the parameter _store_stress_old is true. This parameter can be set from the input file using any of the child classes of ComputeStressBase.

  if (_store_stress_old)
  {
    declarePropertyOld<RankTwoTensor>(_base_name + "stress");
    declarePropertyOlder<RankTwoTensor>(_base_name + "stress");
  }

}

void
ComputeStressBaseParakeet::initQpStatefulProperties()
{
  _stress[_qp].zero();
}

void
ComputeStressBaseParakeet::computeQpProperties()
{
  computeQpStress();
}
