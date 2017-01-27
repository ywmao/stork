/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ArrudaBoyceDamageParakeet.h"

template<>
InputParameters validParams<ArrudaBoyceDamageParakeet>()
{
  InputParameters params = validParams<ComputeStressDamageBaseParakeet>();
  params.addClassDescription("Compute stress with damage using arruda-boyce model");
  return params;
}

ArrudaBoyceDamageParakeet::ArrudaBoyceDamageParakeet(const InputParameters & parameters) :
    ComputeStressDamageBaseParakeet(parameters),
    _stress_old(declarePropertyOld<RankTwoTensor>(_base_name + "stress")),
    _dt_ins(_fe_problem.dt()),
    _G0(getMaterialPropertyByName<Real>(_base_name + "G0")),
    _lmd_L(getMaterialPropertyByName<Real>(_base_name + "lmd_L")),
    _kbulk(getMaterialPropertyByName<Real>(_base_name + "kbulk")),
    _PsiCrit(getMaterialPropertyByName<Real>(_base_name + "PsiCrit")),
    _ell(getMaterialPropertyByName<Real>(_base_name + "ell")),
    _BetaCoeff(getMaterialPropertyByName<Real>(_base_name + "BetaCoeff"))
{
}

void
ArrudaBoyceDamageParakeet::computeQpStress()
{
  // In this function we will calculate.
  // TR (_stress[_qp]); _Jacobian_mult[_qp];
  // Rd_term (1,2,3); Kdd_term (1,2,3)
  // Kdi_term1; Kid_term1; 


  RankTwoTensor F_tau = _F_matrix[_qp];
  RankTwoTensor B_tau = _B_matrix[_qp];

  Real detB = B_tau.det();
  Real trB  = B_tau.trace();
  Real detF = sqrt(detB);

  RankTwoTensor Finv = F_tau.inverse();

  RankTwoTensor Iden;
  Iden.zero();
  Iden.addIa(1.0);

  RankTwoTensor B_dis = B_tau/detB;
  Real lmd_bar = sqrt(B_dis.trace()/3.0);

  Real rel_sth = lmd_bar/_lmd_L[_qp];
  Real beta = rel_sth*(3.0-rel_sth*rel_sth)/(1.0-rel_sth*rel_sth);
  Real G_bar = _G0[_qp]*(1./(3.0*rel_sth))*beta;

  // calculate the piola stress.
  RankTwoTensor DJDF;
  for(unsigned int i=0;i<3;i++)
     for(unsigned int j=0;j<3;j++)
           DJDF(i,j) = detF*Finv(j,i);

  RankTwoTensor DlmdbarDF;
  for(unsigned int i=0;i<3;i++)
     for(unsigned int j=0;j<3;j++)
          DlmdbarDF(i,j) = (1/(3.0*lmd_bar))*pow(detF,-0.667)*(F_tau(i,j)-0.333*trB*Finv(j,i));


  // calculate the _Kid_term1 and the stress.
  Real gfunc = (1.0-_damage[_qp])*(1.0-_damage[_qp])+0.0001;
  Real dgfunc = -2.0*(1-_damage[_qp]);
  for(unsigned int i=0;i<3;i++)
     for(unsigned int j=0;j<3;j++)
     {
          if(detF>=1)
          {
              _stress[_qp](i,j) = gfunc*(_G0[_qp]*_lmd_L[_qp]*beta*DlmdbarDF(i,j)+_kbulk[_qp]*log(detF)*Finv(j,i));
              _Kid_term1[_qp](i,j) = dgfunc*(_G0[_qp]*_lmd_L[_qp]*beta*DlmdbarDF(i,j)+_kbulk[_qp]*log(detF)*Finv(j,i));
          }else{
              _stress[_qp](i,j) = gfunc*(_G0[_qp]*_lmd_L[_qp]*beta*DlmdbarDF(i,j))+_kbulk[_qp]*log(detF)*Finv(j,i);
              _Kid_term1[_qp](i,j) = dgfunc*(_G0[_qp]*_lmd_L[_qp]*beta*DlmdbarDF(i,j));
          }

     }


  //calculate the Jacobian Kii related term.

 RankFourTensor DFinvDF;
 //DFinvDF(i,j,k,l) = DFinv(i,j)DF(k,l)
 for(unsigned int i=0;i<3;i++)
    for(unsigned int j=0;j<3;j++)
       for(unsigned int k=0;k<3;k++)
          for(unsigned int l=0;l<3;l++)
               DFinvDF(i,j,k,l) = -Finv(l,j)*Finv(i,k);

 RankFourTensor D2lmdbarDF2;
 for(unsigned int i=0;i<3;i++)
    for(unsigned int j=0;j<3;j++)
       for(unsigned int k=0;k<3;k++)
          for(unsigned int l=0;l<3;l++)
          {
             Real A1 = -DlmdbarDF(i,j)*DlmdbarDF(k,l)/lmd_bar;
             Real A2 = -0.667*Finv(l,k)*DlmdbarDF(i,j);
             Real A3 = Iden(i,k)*Iden(j,l)-0.667*Finv(j,i)*F_tau(k,l);
                  A3 = A3-0.333*trB*DFinvDF(j,i,k,l);
                  A3 = A3*pow(detF,-0.666)/(3*lmd_bar);
             D2lmdbarDF2(i,j,k,l) = A1+A2+A3;
          }

 Real RS = (rel_sth*rel_sth);
 Real DbetaDlmdbar;
 DbetaDlmdbar = beta/lmd_bar+((4*RS)/((1-RS)*(1-RS)))/_lmd_L[_qp];

  RankFourTensor dTRdF;
  for(unsigned int i=0;i<3;i++)
     for(unsigned int j=0;j<3;j++)
        for(unsigned int l=0;l<3;l++)
           for(unsigned int k=0;k<3;k++)
           {
                Real A1 = _kbulk[_qp]*Finv(j,i)*Finv(l,k);
                Real A2 = -_kbulk[_qp]*log(detF)*Finv(l,i)*Finv(j,k);
                Real A3 = _G0[_qp]*_lmd_L[_qp]*DbetaDlmdbar*DlmdbarDF(k,l)*DlmdbarDF(i,j);
                Real A4 = (_G0[_qp]*_lmd_L[_qp]*beta)*D2lmdbarDF2(i,j,k,l);
                if(detF>=1)
                   dTRdF(i,j,k,l) = gfunc*(A1+A2+A3+A4);
                else
                   dTRdF(i,j,k,l) = (A1+A2)+gfunc*(A3+A4);
           }

  // apply dTdF to jacobian.
  _Jacobian_mult[_qp] = dTRdF;

  //calculte the Rd related terms.
  // calculate _Rd_term1 and _Rd_term2
  _Rd_term1[_qp] = _BetaCoeff[_qp]*_ddamage_dt[_qp];
  _Rd_term2[_qp] = _PsiCrit[_qp]*_ell[_qp]*_ell[_qp]*_grad_damage[_qp];

  // calculate _Rd_term3;
  Real AH = rel_sth*beta+log(beta/sinh(beta));
  Real psi_dev_tau = _G0[_qp]*_lmd_L[_qp]*_lmd_L[_qp]*AH;
  Real psi_vol_tau = 0.5*_kbulk[_qp]*(log(detF))*(log(detF));
  Real psi_tau = 0;
  Real strain_energy;
  if(detF>=1)
  {
      psi_tau = (psi_dev_tau+psi_vol_tau);
      strain_energy =  gfunc*(psi_dev_tau+psi_vol_tau);
  }else{
      psi_tau = psi_dev_tau;
      strain_energy = gfunc*psi_dev_tau+psi_vol_tau;
  }

  Real psi_t = _state_var_psi_max_old[_qp];
  psi_tau = (psi_tau>psi_t)?psi_tau:psi_t;
  // update the state variables.
  _state_var_psi_max[_qp] = psi_tau;

  // numerical approximation.
  //   use the old strain energy to decouple equation.
  Real psi = psi_tau;
  // avoid the matrix to be singular.
  psi = psi +0.00000001*_PsiCrit[_qp];
  Real H_tau = ((psi-_PsiCrit[_qp])>0.0)?(psi-_PsiCrit[_qp]):0.0; 
  _Rd_term3[_qp] = 2*_damage[_qp]*_PsiCrit[_qp]-2*(1.0-_damage[_qp])*H_tau;

 // calculate Kdd related term
 _Kdd_term1[_qp] = _BetaCoeff[_qp]/_dt_ins;
 _Kdd_term2[_qp] = _PsiCrit[_qp]*_ell[_qp]*_ell[_qp];
 _Kdd_term3[_qp] = 2*_PsiCrit[_qp]+2*H_tau; 

 // calculate Kdi related term
 _Kdi_term1[_qp].zero();
}
