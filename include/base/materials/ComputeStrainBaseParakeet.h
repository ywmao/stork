/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTESTRAINBASEPARAKEET_H
#define COMPUTESTRAINBASEPARAKEET_H

// these header files are all defined in moose framework
#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "DerivativeMaterialInterface.h"

/**
 * ComputeStrainBaseParakeet is the base class for strain tensors in Parakeet application
 * Since we mainly interested in the finite strain, so we only return F.
 */
class ComputeStrainBaseParakeet : public DerivativeMaterialInterface<Material>
{
public:
  ComputeStrainBaseParakeet(const InputParameters & parameters);
  virtual ~ComputeStrainBaseParakeet() {}

protected:
  virtual void initQpStatefulProperties();

  /// Coupled displacement variables
  unsigned int _ndisp;
  std::vector<const VariableValue *> _disp;
  std::vector<const VariableGradient *> _grad_disp;

  std::string _base_name;

 // following items are modifieyd in Parakeet since
 // in this class we mainly use deformation gradient based mechanics.
 // volumetric locking may always set to true since we mainly want to
 // understand the mechanics of soft materials.
  MaterialProperty<RankTwoTensor> & _B_matrix; 
  MaterialProperty<RankTwoTensor> & _F_matrix; 

  const bool _stateful_displacements;
  const bool _stateful_deformation_gradient;

  std::vector<const VariableGradient *> _grad_disp_old;

  MaterialProperty<RankTwoTensor> & _deformation_gradient;
  const MaterialProperty<RankTwoTensor> * _deformation_gradient_old;

};

#endif //COMPUTESTRAINBASEPARAKEET_H
