/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef RANKTWOAUXPARAKEET_H
#define RANKTWOAUXPARAKEET_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

class RankTwoAuxParakeet;

/**
 * RankTwoAuxParakeet is designed to take the data in the RankTwoTensor material
 * property, for example stress or strain, and output the value for the
 * supplied indices.
 */

template<>
InputParameters validParams<RankTwoAuxParakeet>();

class RankTwoAuxParakeet : public AuxKernel
{
public:
  RankTwoAuxParakeet(const InputParameters & parameters);

protected:
  virtual Real computeValue();

private:
  const MaterialProperty<RankTwoTensor> & _tensor;
  const unsigned int _i;
  const unsigned int _j;
};

#endif //RANKTWOAUXPARAKEET_H
