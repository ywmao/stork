/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ALEKERNELPARAKEET_H
#define ALEKERNELPARAKEET_H

#include "Kernel.h"
#include "Assembly.h"

class ALEKernelParakeet;

template<>
InputParameters validParams<ALEKernelParakeet>();

class ALEKernelParakeet : public Kernel
{
public:
  ALEKernelParakeet(const InputParameters & parameters);

protected:
  virtual void computeJacobian();
  virtual void computeOffDiagJacobian(unsigned int jvar);

  /// undisplaced problem
  Assembly & _assembly_undisplaced;

  /// Reference to this Kernel's undisplaced MooseVariable object
  MooseVariable & _var_undisplaced;

  ///@{ Shape and test functions on the undisplaced mesh
  const VariablePhiGradient & _grad_phi_undisplaced;
  const VariableTestGradient & _grad_test_undisplaced;
  ///@}
};

#endif //ALEKERNELPARAKEET_H
