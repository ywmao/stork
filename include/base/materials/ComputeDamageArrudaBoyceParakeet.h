/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTEDAMAGEARRUDABOYCEPARAKEET_H
#define COMPUTEDAMAGEARRUDABOYCEPARAKEET_H

#include "ComputeDamageParakeet.h"

/**
 * ComputeDamageArrudaBoyceParakeet defines related damage.
 */

class ComputeDamageArrudaBoyceParakeet : public ComputeDamageParakeet
{
public:
  ComputeDamageArrudaBoyceParakeet(const InputParameters & parameters);

  virtual void computeProperties();

protected:

};

#endif //COMPUTEDAMAGEPARAKEET_H
