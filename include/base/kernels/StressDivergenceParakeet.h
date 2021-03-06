/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef STRESSDIVERGENCEPARAKEET_H
#define STRESSDIVERGENCEPARAKEET_H

#include "ALEKernelParakeet.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

//Forward Declarations
class StressDivergenceParakeet;
class RankTwoTensor;
class RankFourTensor;

template<>
InputParameters validParams<StressDivergenceParakeet>();

/**
 * StressDivergenceTensorsParakeet mostly copies from StressDivergence.  There are small changes to use
 * RankFourTensor and RankTwoTensors instead of SymmElasticityTensors and SymmTensors.  This is done
 * to allow for more mathematical transparancy.
 */
class StressDivergenceParakeet : public ALEKernelParakeet
{
public:
  StressDivergenceParakeet(const InputParameters & parameters);

protected:

  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  virtual void computeJacobian();
  virtual void computeOffDiagJacobian(unsigned int jvar);

  //virtual void computeFiniteDeformJacobian();

  std::string _base_name;
  bool _use_finite_deform_jacobian;

  const MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RankFourTensor> & _Jacobian_mult;

  const unsigned int _component;

  /// Coupled displacement variables
  unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;

  //const bool _temp_coupled;

  //const unsigned int _temp_var;
};

#endif //STRESSDIVERGENCEPARAKEET_H
