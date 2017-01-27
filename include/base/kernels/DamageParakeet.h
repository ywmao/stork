/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef DAMAGEPARAKEET_H
#define DAMAGEPARAKEET_H

#include "ALEKernelParakeet.h"
#include "RankTwoTensor.h"

//Forward Declarations
class DamageParakeet;
class RankTwoTensor;

template<>
InputParameters validParams<DamageParakeet>();

/**
 * DamageParakeet mostly copies from StressDivergence.
 */
class DamageParakeet : public ALEKernelParakeet
{
public:
  DamageParakeet(const InputParameters & parameters);

protected:

  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  virtual void computeJacobian();
  virtual void computeOffDiagJacobian(unsigned int jvar);

  std::string _base_name;

  const  MaterialProperty<Real> & _Rd_term1;
  const MaterialProperty<RealGradient> & _Rd_term2;
  const MaterialProperty<Real> & _Rd_term3;
  const MaterialProperty<Real> & _Kdd_term1;
  const MaterialProperty<Real> & _Kdd_term2;
  const MaterialProperty<Real> & _Kdd_term3;
  const MaterialProperty<RankTwoTensor> & _Kdi_term1;

  //material parameters
  //const MaterialProperty<Real> & _G0;
  //const MaterialProperty<Real> & _lmd_L;
  //const MaterialProperty<Real> & _kbulk;
  //const MaterialProperty<Real> & _PsiCrit;
  //const MaterialProperty<Real> & _ell;
  //const MaterialProperty<Real> & _BetaCoeff;

  //coupled displacement var.
  unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;
};

#endif //DAMAGEPARAKEET_H
