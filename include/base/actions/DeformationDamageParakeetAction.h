/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef DEFORMATIONDAMAGEPARAKEETACTION_H
#define DEFORMATIONDAMAGEPARAKEETACTION_H

#include "Action.h"

class DeformationDamageParakeetAction;

template<>
InputParameters validParams<DeformationDamageParakeetAction>();

class DeformationDamageParakeetAction : public Action
{
public:
  DeformationDamageParakeetAction(const InputParameters & params);

  virtual void act();

protected:
  virtual std::string getKernelType();
  virtual InputParameters getParameters(std::string type);

  std::vector<NonlinearVariableName> _displacements;
  unsigned int _ndisp;
  std::vector<VariableName> _coupled_displacements;

  std::vector<NonlinearVariableName> _damageInpVar;
  unsigned int _ndamage;
  std::vector<VariableName> _coupled_damageInpVar;

  //std::vector<AuxVariableName> _save_in;
  std::vector<AuxVariableName> _diag_save_in;

  // see the PressureAction.h in TensorMechanics modulus for details.
  // this is supposed to save the displacement and force.
  std::vector<std::vector<AuxVariableName> > _save_in_vars;
  std::vector<bool> _has_save_in_vars;

  Moose::CoordinateSystemType _coord_system;
};

#endif //DEFORMATIONDAMAGEPARAKEETACTION_H
