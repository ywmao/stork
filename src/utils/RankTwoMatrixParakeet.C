/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "RankTwoMatrixParakeet.h"
#include "MooseException.h"
#include "MatrixTools.h"
#include "MaterialProperty.h"

// Any other includes here
#include "libmesh/utility.h"
#include <ostream>

template<>
void mooseSetToZero<RankTwoMatrixParakeet>(RankTwoMatrixParakeet & v)
{
  v.zero();
}

template<>
void
dataStore(std::ostream & stream, RankTwoMatrixParakeet & rft, void * context)
{
  dataStore(stream, rft._vals, context);
}

template<>
void
dataLoad(std::istream & stream, RankTwoMatrixParakeet & rft, void * context)
{
  dataLoad(stream, rft._vals, context);
}

RankTwoMatrixParakeet::RankTwoMatrixParakeet()
{
  //mooseAssert(N == 3, "RankThreeMatrix is currently only tested for 3 dimensions.");

  for (unsigned int i = 0; i < N+1; ++i)
    for (unsigned int j = 0; j < N+1; ++j)
          _vals[i][j] = 0.0;
}

Real &
RankTwoMatrixParakeet::operator()(unsigned int i, unsigned int j)
{
  return _vals[i][j];
}

Real
RankTwoMatrixParakeet::operator()(unsigned int i, unsigned int j) const
{
  return _vals[i][j];
}

void
RankTwoMatrixParakeet::zero()
{
  for (unsigned int i = 0; i < N+1; ++i)
    for (unsigned int j = 0; j < N+1; ++j)
          _vals[i][j] = 0.0;
}

RankTwoMatrixParakeet &
RankTwoMatrixParakeet::operator=(const RankTwoMatrixParakeet & a)
{
  for (unsigned int i = 0; i < N+1; ++i)
    for (unsigned int j = 0; j < N+1; ++j)
          _vals[i][j] = a(i,j);

  return *this;
}

void
RankTwoMatrixParakeet::print(std::ostream & stm) const
{
  const RankTwoMatrixParakeet & a = *this;

  for (unsigned int i = 0; i < N+1; ++i)
    for (unsigned int j = 0; j < N+1; ++j)
    {
      stm << "i = " << i << " j = " << j << '\n';
          stm << std::setw(15) << a(i,j) << " ";

        stm << '\n';
      }
}
