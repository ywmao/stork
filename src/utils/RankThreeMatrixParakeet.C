/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "RankThreeMatrixParakeet.h"
#include "MooseException.h"
#include "MatrixTools.h"
#include "MaterialProperty.h"

// Any other includes here
#include "libmesh/utility.h"
#include <ostream>

template<>
void mooseSetToZero<RankThreeMatrixParakeet>(RankThreeMatrixParakeet & v)
{
  v.zero();
}

template<>
void
dataStore(std::ostream & stream, RankThreeMatrixParakeet & rft, void * context)
{
  dataStore(stream, rft._vals, context);
}

template<>
void
dataLoad(std::istream & stream, RankThreeMatrixParakeet & rft, void * context)
{
  dataLoad(stream, rft._vals, context);
}

RankThreeMatrixParakeet::RankThreeMatrixParakeet()
{
  //mooseAssert(N == 3, "RankThreeMatrix is currently only tested for 3 dimensions.");

  for (unsigned int i = 0; i < N; ++i)
    for (unsigned int j = 0; j < N+1; ++j)
      for (unsigned int k = 0; k < NB; ++k)
          _vals[i][j][k] = 0.0;
}

Real &
RankThreeMatrixParakeet::operator()(unsigned int i, unsigned int j, unsigned int k)
{
  return _vals[i][j][k];
}

Real
RankThreeMatrixParakeet::operator()(unsigned int i, unsigned int j, unsigned int k) const
{
  return _vals[i][j][k];
}

void
RankThreeMatrixParakeet::zero()
{
  for (unsigned int i = 0; i < N; ++i)
    for (unsigned int j = 0; j < N+1; ++j)
      for (unsigned int k = 0; k < NB; ++k)
          _vals[i][j][k] = 0.0;
}

RankThreeMatrixParakeet &
RankThreeMatrixParakeet::operator=(const RankThreeMatrixParakeet & a)
{
  for (unsigned int i = 0; i < N; ++i)
    for (unsigned int j = 0; j < N+1; ++j)
      for (unsigned int k = 0; k < NB; ++k)
          _vals[i][j][k] = a(i,j,k);

  return *this;
}

void
RankThreeMatrixParakeet::print(std::ostream & stm) const
{
  const RankThreeMatrixParakeet & a = *this;

  for (unsigned int i = 0; i < N; ++i)
    for (unsigned int j = 0; j < N+1; ++j)
    {
      stm << "i = " << i << " j = " << j << '\n';
      for (unsigned int k = 0; k < NB; ++k)
      {
          stm << std::setw(15) << a(i,j,k) << " ";

        stm << '\n';
      }
    }
}
