/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ARRUDABOYCEPARAKEET_H
#define ARRUDABOYCEPARAKEET_H

#include "ComputeStressBaseParakeet.h"

/**
 * ArrudaBoyceParakeet computes the stress using AB model with three materials parameters
 * as input
 */
class ArrudaBoyceParakeet : public ComputeStressBaseParakeet
{
public:
  ArrudaBoyceParakeet(const InputParameters & parameters);

protected:
  virtual void computeQpStress();

  MaterialProperty<RankTwoTensor> & _stress_old;
  const MaterialProperty<RankTwoTensor> & _B_matrix;
  const MaterialProperty<RankTwoTensor> & _F_matrix;

  const MaterialProperty<Real> & _G0;
  const MaterialProperty<Real> & _lmd_L;
  const MaterialProperty<Real> & _kbulk;
};

#endif //ARRUDABOYCEPARAKEET_H
