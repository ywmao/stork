/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTESTRESSBASEPARAKEET_H
#define COMPUTESTRESSBASEPARAKEET_H

// all these header files are defined in moose framework
#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "DerivativeMaterialInterface.h"

// for our own purpose we add two classes.
#include "RankTwoMatrixParakeet.h"
#include "RankThreeMatrixParakeet.h"


/**
 * ComputeStressBaseParakeet is the base class for stress tensors
 */
class ComputeStressBaseParakeet : public DerivativeMaterialInterface<Material>
{
public:
  ComputeStressBaseParakeet(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();
  virtual void computeQpStress() = 0;

  std::string _base_name;

  const MaterialProperty<RankTwoTensor> & _B_matrix;
  const MaterialProperty<RankTwoTensor> & _F_matrix;

  MaterialProperty<RankTwoTensor> & _stress;

  /// derivative of stress w.r.t. strain (_dstress_dstrain)
  MaterialProperty<RankFourTensor> & _Jacobian_mult;

  /// Parameter which decides whether to store old stress.
  /// This is required for HHT time integration and Rayleigh damping
  const bool _store_stress_old;
};

#endif //COMPUTESTRESSBASEPARAKEET_H
