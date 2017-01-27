#include "ComputeArrudaBoyceDamageMPBaseParakeet.h"
//#include "Function.h"

template<>
InputParameters validParams<ComputeArrudaBoyceDamageMPBaseParakeet>()
{
  InputParameters params = validParams<Material>();
  params.addParam<std::string>("base_name", "Optional parameter that allows the user to define multiple mechanics material systems on the same block, i.e. for multiple phases");
  return params;
}

ComputeArrudaBoyceDamageMPBaseParakeet::ComputeArrudaBoyceDamageMPBaseParakeet(const InputParameters & parameters) :
    DerivativeMaterialInterface<Material>(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "" ),
    _G0_name(_base_name + "G0"),
    _lmd_L_name(_base_name + "lmd_L"),
    _kbulk_name(_base_name + "kbulk"),
    _PsiCrit_name(_base_name+"PsiCrit"),
    _ell_name(_base_name+"ell"),
    _BetaCoeff_name(_base_name+"BetaCoeff"),
    _G0(declareProperty<Real>(_G0_name)),
    _lmd_L(declareProperty<Real>(_lmd_L_name)),
    _kbulk(declareProperty<Real>(_kbulk_name)),
    _PsiCrit(declareProperty<Real>(_PsiCrit_name)),
    _ell(declareProperty<Real>(_ell_name)),
    _BetaCoeff(declareProperty<Real>(_BetaCoeff_name))
{
}

void
ComputeArrudaBoyceDamageMPBaseParakeet::computeQpProperties()
{
  computeQpMaterialProperties();
}
