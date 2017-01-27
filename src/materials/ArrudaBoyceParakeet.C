/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ArrudaBoyceParakeet.h"

template<>
InputParameters validParams<ArrudaBoyceParakeet>()
{
  InputParameters params = validParams<ComputeStressBaseParakeet>();
  params.addClassDescription("Compute stress using arruda-boyce model");
  return params;
}

ArrudaBoyceParakeet::ArrudaBoyceParakeet(const InputParameters & parameters) :
    ComputeStressBaseParakeet(parameters),
    _stress_old(declarePropertyOld<RankTwoTensor>(_base_name + "stress")),
    _B_matrix(getMaterialPropertyByName<RankTwoTensor>(_base_name + "B_matrix")),
    _F_matrix(getMaterialPropertyByName<RankTwoTensor>(_base_name + "F_matrix")),
    _G0(getMaterialPropertyByName<Real>(_base_name + "G0")),
    _lmd_L(getMaterialPropertyByName<Real>(_base_name + "lmd_L")),
    _kbulk(getMaterialPropertyByName<Real>(_base_name + "kbulk"))
{
}

void
ArrudaBoyceParakeet::computeQpStress()
{

  //Take care that we use the undisplaced mesh
  // as the gradient variables. Then we need to use
  // T_R instead of T. So _stress = Piola stress.

  // Also the Jacobian_mult = d Piola / d F.

  // for the case for using displaced mesh.
  //   please see the code in kangaroo application.

  // compute the stress (Piola).
  //RankTwoTensor Bm = static_cast<RankTwoTensor>(_B_matrix);
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

  // This stress is the Piola stress.
  for(unsigned int i=0;i<3;i++)
     for(unsigned int j=0;j<3;j++) 
          _stress[_qp](i,j) = _G0[_qp]*_lmd_L[_qp]*beta*DlmdbarDF(i,j)+_kbulk[_qp]*log(detF)*Finv(j,i);

  // push the Piola->Cauchy.
//  RankTwoTensor cauchy_stress;
//  cauchy_stress.zero();
//  cauchy_stress = _stress[_qp]*(_F_matrix[_qp].transpose())/detF;
//  _stress[_qp] = cauchy_stress;
 
  // compute the _jacobian_mult

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
		dTRdF(i,j,k,l) = (A1+A2+A3+A4);
           }
  
  // apply dTdF to jacobian.
  _Jacobian_mult[_qp] = dTRdF;

  // push from material to spatial.
 // RankFourTensor AA;
 // AA.zero();
 //   for(unsigned int i=0;i<3;i++)
 //    for(unsigned int j=0;j<3;j++)
 //       for(unsigned int k=0;k<3;k++)
 //          for(unsigned int l=0;l<3;l++)
 //             for(unsigned int m=0;m<3;m++)
 //                for(unsigned int n=0;n<3;n++)
 //                {
 //                   AA(i,j,k,l) = AA(i,j,k,l)+dTRdF(i,m,k,n)*_F_matrix[_qp](l,n)*_F_matrix[_qp](j,m)/detF;
 //                }
 // _Jacobian_mult[_qp] = AA;
}
