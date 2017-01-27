/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTEARRUDABOYCEDAMAGEMPBASEPARAKEET_H
#define COMPUTEARRUDABOYCEDAMAGEMPBASEPARAKEET_H

#include "Material.h"
#include "RankFourTensor.h"

/**
 * ComputeArrudaBoyceDamageMPBaseParakeet the base class for computing Arruda Boyce model material parameter. Here we will use the damage part.
 */

class ComputeArrudaBoyceDamageMPBaseParakeet : public DerivativeMaterialInterface<Material>
{
public:
  ComputeArrudaBoyceDamageMPBaseParakeet(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void computeQpMaterialProperties() = 0;

  std::string _base_name;
  std::string _G0_name;
  std::string _lmd_L_name;
  std::string _kbulk_name;
  std::string _PsiCrit_name;
  std::string _ell_name;
  std::string _BetaCoeff_name;

  MaterialProperty<Real> & _G0;
  MaterialProperty<Real> & _lmd_L;
  MaterialProperty<Real> & _kbulk;
  MaterialProperty<Real> & _PsiCrit;
  MaterialProperty<Real> & _ell;
  MaterialProperty<Real> & _BetaCoeff;

};

#endif //COMPUTEARRUDABOYCEDAMAGEMPBASEPARAKEET_H
