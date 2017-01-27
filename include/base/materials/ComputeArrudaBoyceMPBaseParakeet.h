/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTEARRUDABOYCEMPBASEPARAKEET_H
#define COMPUTEARRUDABOYCEMPBASEPARAKEET_H

#include "Material.h"
#include "RankFourTensor.h"

/**
 * ComputeArrudaBoyceMPBaseParakeet the base class for computing Arruda Boyce model material parameter.
 */

class ComputeArrudaBoyceMPBaseParakeet : public DerivativeMaterialInterface<Material>
{
public:
  ComputeArrudaBoyceMPBaseParakeet(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void computeQpMaterialProperties() = 0;

  std::string _base_name;
  std::string _G0_name;
  std::string _lmd_L_name;
  std::string _kbulk_name;

  MaterialProperty<Real> & _G0;
  MaterialProperty<Real> & _lmd_L;
  MaterialProperty<Real> & _kbulk;

};

#endif //COMPUTEARRUDABOYCEMPBASEPARAKEET_H
