/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ARRUDABOYCEDAMAGEPARAKEET_H
#define ARRUDABOYCEDAMAGEPARAKEET_H

#include "ComputeStressDamageBaseParakeet.h"

/**
 * ArrudaBoyceDamageParakeet computes the stress 
 * using AB model with damage parameter; with three materials parameters
 * as input
 */

/**
 * for the details about the state variables, see solid mechanics/materials/abaqus.
 *  Another possible source is in tensor mechanics/mechanics/FiniteStrainPlasticMaterial
 */
/**
 * Actually more direct one is the Example 09. But take care look at the updated one.
 */


class ArrudaBoyceDamageParakeet : public ComputeStressDamageBaseParakeet
{
public:
  ArrudaBoyceDamageParakeet(const InputParameters & parameters);

protected:
  virtual void computeQpStress();

  //virtual void initQpStatefulProperties();

  MaterialProperty<RankTwoTensor> & _stress_old;
  //const MaterialProperty<RankTwoTensor> & _B_matrix;
  //const MaterialProperty<RankTwoTensor> & _F_matrix;
  //const MaterialProperty<Real> & _damage;
  //const MaterialProperty<RealGradient> & _grad_damage;
  //const MaterialProperty<Real> & _ddamage_dt;
  const Real & _dt_ins;
  const MaterialProperty<Real> & _G0;
  const MaterialProperty<Real> & _lmd_L;
  const MaterialProperty<Real> & _kbulk;
  const MaterialProperty<Real> & _PsiCrit;
  const MaterialProperty<Real> & _ell;
  const MaterialProperty<Real> & _BetaCoeff;
  //MaterialProperty<Real> & _state_var_psi_max;
  //MaterialProperty<Real> & _state_var_psi_max_old;
};

#endif //ARRUDABOYCEDAMAGEPARAKEET_H
