/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTEARRUDABOYCEDAMAGEMPPARAKEET_H
#define COMPUTEARRUDABOYCEDAMAGEMPPARAKEET_H

#include "ComputeArrudaBoyceDamageMPBaseParakeet.h"

/**
 * ComputeArrudaBoyceDamageMPParakeet defines an elasticity tensor material for
 * AB model.
 */
class ComputeArrudaBoyceDamageMPParakeet : public ComputeArrudaBoyceDamageMPBaseParakeet
{
public:
  ComputeArrudaBoyceDamageMPParakeet(const InputParameters & parameters);

protected:
  virtual void computeQpMaterialProperties();

  /// Elastic constants
  bool _G0R_set;
  bool _lmd_LR_set;
  bool _kbulkR_set;
  bool _PsiCritR_set;
  bool _ellR_set;
  bool _BetaCoeffR_set;

  Real _G0R;
  Real _lmd_LR;
  Real _kbulkR;
  Real _PsiCritR;
  Real _ellR;
  Real _BetaCoeffR;
 //// R here for "Read from input file"

};

#endif //COMPUTEARRUDABOYCEDAMAGEMPPARAKEET_H
