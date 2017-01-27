/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef RANKTWOMATRIXPARAKEET_H
#define RANKTWOMATRIXPARAKEET_H

// Forward declaration of RankTwoMatrixParakeet
class RankTwoMatrixParakeet;

// MOOSE includes
#include "Moose.h"
#include "PermutationTensor.h"
#include "MooseEnum.h"
#include "DerivativeMaterialInterface.h"

// libMesh includes
#include "libmesh/tensor_value.h"
#include "libmesh/libmesh.h"
#include "libmesh/vector_value.h"

/**
 * Helper function template specialization to set an object to zero.
 * Needed by DerivativeMaterialInterface
 */
template<>
void mooseSetToZero<RankTwoMatrixParakeet>(RankTwoMatrixParakeet & v);

/**
 * RankTwoMatrixParakeet is designed to handle any 3-dimensional two order matrix, p_matrix.
 *
 * It is designed to allow for maximum clarity of the mathematics and ease of use.
 * Original class authors: Yunwei Mao
 *
 * P(NDim+1,NDim+1) : Maximum Q(4,4); 4*4 = 16 entries.
 */
class RankTwoMatrixParakeet
{
public:

  /// Default constructor; fills to zero
  RankTwoMatrixParakeet();

  /// Gets the value for the index specified.
  /// i = 0,1,2; j=0,1,2,3; k=0,1,...,26
  Real & operator()(unsigned int i, unsigned int j);

  /**
   * Gets the value for the index specified.
   * used for const
   /// i = 0,1,2,3; j=0,1,2,3;
   */
  Real operator()(unsigned int i, unsigned int j) const;

  /// Zeros out the tensor.
  void zero();

  /// Print the rank two matrix
  void print(std::ostream & stm = Moose::out) const;

  /// copies values from a into this tensor
  RankTwoMatrixParakeet & operator= (const RankTwoMatrixParakeet & a);


protected:
  /// Dimensionality of rank-three matrix
  static const unsigned int N = LIBMESH_DIM;

  /// The values of the rank-three matrix
  Real _vals[N+1][N+1];


  template<class T>
  friend void dataStore(std::ostream &, T &, void *);

  template<class T>
  friend void dataLoad(std::istream &, T &, void *);
};

template<>
void dataStore(std::ostream &, RankTwoMatrixParakeet &, void *);

template<>
void dataLoad(std::istream &, RankTwoMatrixParakeet &, void *);

#endif //RANKTWOMATRIXPARAKEET_H
