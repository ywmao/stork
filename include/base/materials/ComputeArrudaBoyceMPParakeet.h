/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTEARRUDABOYCEMPPARAKEET_H
#define COMPUTEARRUDABOYCEMPPARAKEET_H

#include "ComputeArrudaBoyceMPBaseParakeet.h"

/**
 * ComputeArrudaBoyceMPParakeet defines an elasticity tensor material for
 * AB model.
 */
class ComputeArrudaBoyceMPParakeet : public ComputeArrudaBoyceMPBaseParakeet
{
public:
  ComputeArrudaBoyceMPParakeet(const InputParameters & parameters);

protected:
  virtual void computeQpMaterialProperties();

  /// Elastic constants
  bool _G0R_set;
  bool _lmd_LR_set;
  bool _kbulkR_set;

  Real _G0R;
  Real _lmd_LR;
  Real _kbulkR;
 //// R here for "Read from input file"

};

#endif //COMPUTEARRUDABOYCEMPPARAKEET_H
