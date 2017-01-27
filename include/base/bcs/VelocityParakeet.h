/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef VELOCITYPARAKEET_H
#define VELOCITYPARAKEET_H

#include "PresetNodalBC.h"


class VelocityParakeet : public PresetNodalBC
{
public:
  VelocityParakeet(const InputParameters & parameters);

protected:
  virtual Real computeQpValue();

  const VariableValue & _u_old;
  const Real _velocity;
  Function & _function;
};

template<>
InputParameters validParams<VelocityParakeet>();

#endif /* VELOCITYPARAKEET_H */
