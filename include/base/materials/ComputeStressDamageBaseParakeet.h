/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTESTRESSDAMAGEBASEPARAKEET_H
#define COMPUTESTRESSDAMAGEBASEPARAKEET_H

// all these header files are defined in moose framework
#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "DerivativeMaterialInterface.h"

// for our own purpose we add two classes.
#include "RankTwoMatrixParakeet.h"
#include "RankThreeMatrixParakeet.h"


/**
 * ComputeStressDamageBaseParakeet is the base class for stress tensors
 *  Including the effect of damage.
 */
class ComputeStressDamageBaseParakeet : public DerivativeMaterialInterface<Material>
{
public:
  ComputeStressDamageBaseParakeet(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();
  virtual void computeQpStress() = 0;

  std::string _base_name;

  const MaterialProperty<RankTwoTensor> & _B_matrix;
  const MaterialProperty<RankTwoTensor> & _F_matrix;
  
  MaterialProperty<RankTwoTensor> & _stress;

  /// derivative of stress w.r.t. strain (_dstress_dstrain)
  MaterialProperty<RankFourTensor> & _Jacobian_mult;

  /// Parameter which decides whether to store old stress.
  /// This is required for HHT time integration and Rayleigh damping
  const bool _store_stress_old;

  // coupling damage parameter; and the jacobian related terms.
  const MaterialProperty<Real> & _damage;
  const MaterialProperty<RealGradient> & _grad_damage;
  const MaterialProperty<Real> & _ddamage_dt;
  MaterialProperty<Real> & _Rd_term1;
  MaterialProperty<RealGradient> & _Rd_term2;
  MaterialProperty<Real> & _Rd_term3;
  MaterialProperty<Real> & _Kdd_term1;
  MaterialProperty<Real> & _Kdd_term2;
  MaterialProperty<Real> & _Kdd_term3;
  MaterialProperty<RankTwoTensor> & _Kdi_term1;
  MaterialProperty<RankTwoTensor> & _Kid_term1;

  MaterialProperty<Real> & _state_var_psi_max;
  MaterialProperty<Real> & _state_var_psi_max_old;

};

#endif //COMPUTESTRESSDAMAGEBASEPARAKEET_H
