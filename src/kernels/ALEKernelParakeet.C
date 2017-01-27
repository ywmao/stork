/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ALEKernelParakeet.h"

template<>
InputParameters validParams<ALEKernelParakeet>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Sets up derivetives with respect to initial configuration");
  return params;
}

ALEKernelParakeet::ALEKernelParakeet(const InputParameters & parameters) :
    Kernel(parameters),
    _assembly_undisplaced(_fe_problem.assembly(_tid)),
    _var_undisplaced(_fe_problem.getVariable(_tid, parameters.get<NonlinearVariableName>("variable"))),
    _grad_phi_undisplaced(_assembly_undisplaced.gradPhi()),
    _grad_test_undisplaced(_var_undisplaced.gradPhi())
{
}

void
ALEKernelParakeet::computeJacobian()
{
  _fe_problem.prepareShapes(_var.number(), _tid);
  Kernel::computeJacobian();
}

void
ALEKernelParakeet::computeOffDiagJacobian(unsigned int jvar)
{
  _fe_problem.prepareShapes(jvar, _tid);
  Kernel::computeOffDiagJacobian(jvar);
}
