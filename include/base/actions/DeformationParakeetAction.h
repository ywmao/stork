/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef DEFORMATIONPARAKEETACTION_H
#define DEFORMATIONPARAKEETACTION_H

#include "Action.h"

class DeformationParakeetAction;

template<>
InputParameters validParams<DeformationParakeetAction>();

class DeformationParakeetAction : public Action
{
public:
  DeformationParakeetAction(const InputParameters & params);

  virtual void act();

protected:
  virtual std::string getKernelType();
  virtual InputParameters getParameters(std::string type);

  std::vector<NonlinearVariableName> _displacements;
  unsigned int _ndisp;

  std::vector<VariableName> _coupled_displacements;

  // this needs to be commented out.
  //std::vector<AuxVariableName> _save_in;
  std::vector<AuxVariableName> _diag_save_in;

  // see the PressureAction.h in TensorMechanics modulus for details.
  // this is supposed to save the displacement and force.
  std::vector<std::vector<AuxVariableName> > _save_in_vars;
  std::vector<bool> _has_save_in_vars;

  Moose::CoordinateSystemType _coord_system;

};

#endif //DEFORMATIONPARAKEETACTION_H
