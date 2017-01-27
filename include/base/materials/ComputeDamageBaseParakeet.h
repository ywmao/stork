/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTEDAMAGEBASEPARAKEET_H
#define COMPUTEDAMAGEBASEPARAKEET_H

// these header files are all defined in moose framework
#include "Material.h"
#include "DerivativeMaterialInterface.h"

/**
 * ComputeDamageBaseParakeet is the base class for damage quantity (or quantities)
 *   and its (their) gradient in Parakeet application
 * only coupling the order parameter into the system.
 */
/*
 * This class can be a model from OrderParameterFunctionMaterial.h in
 *   Phase field modulus in moose
 */
class ComputeDamageBaseParakeet : public DerivativeMaterialInterface<Material>
{
public:
  ComputeDamageBaseParakeet(const InputParameters & parameters);
  virtual ~ComputeDamageBaseParakeet() {}

protected:
  virtual void initQpStatefulProperties();
  //virtual void computeProperties();
  //virtual void computeQpProperties();
  //virtual void computeQpDamageEvol() = 0;

  /// Coupled damage variable
  const VariableValue &_order_parameter;
  const VariableGradient &_grad_order_parameter;
  const VariableValue &_deta_dt;

  std::string _base_name;

 // following items are modifieyd in Parakeet since
 // in this class we mainly use deformation gradient based mechanics.
 // volumetric locking may always set to true since we mainly want to
 // understand the mechanics of soft materials.

  //const bool _stateful_damage;

  MaterialProperty<Real> & _damage;
  MaterialProperty<RealGradient> & _grad_damage;
  MaterialProperty<Real> & _ddamage_dt;

  //MaterialProperty<Real> & _state_var_psi_max;
  //MaterialProperty<Real> & _state_var_psi_max_old;

  // driving force for damage order parameter.
  //  and its jacobian.
  //MaterialProperty<Real> & _dpsi_deta;
  //MaterialProperty<Real> & _dpsi2_deta2;

};

#endif //COMPUTEDAMAGEBASEPARAKEET_H
