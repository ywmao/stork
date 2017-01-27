/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTEFINITESTRAINPARAKEET_H
#define COMPUTEFINITESTRAINPARAKEET_H

#include "ComputeStrainBaseParakeet.h"

/**
 * ComputeFiniteStrainParakeet defines related strain/stretch measure for Fe/Fp, for finite strains.
 */

class ComputeFiniteStrainParakeet : public ComputeStrainBaseParakeet
{
public:
  ComputeFiniteStrainParakeet(const InputParameters & parameters);

  virtual void computeProperties();

protected:

};

#endif //COMPUTEFINITESTRAINPARAKEET_H
