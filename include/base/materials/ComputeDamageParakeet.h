/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTEDAMAGEPARAKEET_H
#define COMPUTEDAMAGEPARAKEET_H

#include "ComputeDamageBaseParakeet.h"

/**
 * ComputeDamageParakeet defines related damage.
 *  Transfer the data from _order_parameter into _damage
 */

class ComputeDamageParakeet : public ComputeDamageBaseParakeet
{
public:
  ComputeDamageParakeet(const InputParameters & parameters);

  virtual void computeProperties();

protected:
  //virtual void computeQpProperties();
  //virtual void computeQpDamageEvol();

  //const MaterialProperty<RankTwoTensor> & _B_matrix;
  //const MaterialProperty<RankTwoTensor> & _F_matrix;
  //const MaterialPro
};

#endif //COMPUTEDAMAGEPARAKEET_H
