#include "ParakeetApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

//add include files for pure deformation problem.
#include "DeformationParakeetAction.h"

#include "DisplacementParakeet.h"
#include "VelocityParakeet.h"
#include "PressureParakeet.h"

#include "StressDivergenceParakeet.h"

#include "ComputeFiniteStrainParakeet.h"
#include "ComputeArrudaBoyceMPParakeet.h"
#include "ArrudaBoyceParakeet.h"

//add include files for deformation+damage problem.
// strain part is the same.
#include "ComputeDamageParakeet.h"
#include "ComputeArrudaBoyceDamageMPParakeet.h"
#include "ArrudaBoyceDamagePart1Parakeet.h"
#include "ArrudaBoyceDamagePart2Parakeet.h"
#include "DamageParakeet.h"
//#include "DamageParakeetAction.h"

template<>
InputParameters validParams<ParakeetApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

ParakeetApp::ParakeetApp(InputParameters parameters) :
    MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  ParakeetApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  ParakeetApp::associateSyntax(_syntax, _action_factory);
}

ParakeetApp::~ParakeetApp()
{
}

// External entry point for dynamic application loading
extern "C" void ParakeetApp__registerApps() { ParakeetApp::registerApps(); }
void
ParakeetApp::registerApps()
{
  registerApp(ParakeetApp);
}

// External entry point for dynamic object registration
extern "C" void ParakeetApp__registerObjects(Factory & factory) { ParakeetApp::registerObjects(factory); }
void
ParakeetApp::registerObjects(Factory & factory)
{

      //kernels
      registerKernel(StressDivergenceParakeet);

      //materials
      registerMaterial(ComputeFiniteStrainParakeet);
      registerMaterial(ComputeArrudaBoyceMPParakeet);
      registerMaterial(ArrudaBoyceParakeet);

      //BCs
      registerBoundaryCondition(VelocityParakeet);
      registerBoundaryCondition(PressureParakeet);
      registerBoundaryCondition(DisplacementParakeet);

      // for damage 
      registerKernel(DamageParakeet);
      registerMaterial(ComputeDamageParakeet);
      registerMaterial(ComputeArrudaBoyceDamageMPParakeet);
      registerMaterial(ArrudaBoyceDamagePart1Parakeet);
      registerMaterial(ArrudaBoyceDamagePart2Parakeet);

}

// External entry point for dynamic syntax association
extern "C" void ParakeetApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { ParakeetApp::associateSyntax(syntax, action_factory); }
void
ParakeetApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{

     syntax.registerActionSyntax("DeformationParakeetAction", "Kernels/DeformationParakeet");
     registerAction(DeformationParakeetAction, "add_kernel");
     //syntax.registerActionSyntax("DamageParakeetAction", "Kernels/DamageParakeet");
     //registerAction(DamageParakeetAction, "add_kernel");

}
