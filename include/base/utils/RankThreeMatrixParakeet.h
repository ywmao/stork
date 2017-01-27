/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef RANKTHREEMATRIXPARAKEET_H
#define RANKTHREEMATRIXPARAKEET_H

// Forward declaration of RankThreeMatrixParakeet
class RankThreeMatrixParakeet;

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
void mooseSetToZero<RankThreeMatrixParakeet>(RankThreeMatrixParakeet & v);

/**
 * RankThreeMatrixParakeet is designed to handle any 3-dimensional three order matrix, Q_matrix.
 *
 * It is designed to allow for maximum clarity of the mathematics and ease of use.
 * Original class authors: Yunwei Mao
 *
 * Q(NDim,NDim+1,Nnode) : Maximum Q(3,4,27); 12*27 = 324 entries.
 */
class RankThreeMatrixParakeet
{
public:

  /// Default constructor; fills to zero
  RankThreeMatrixParakeet();

  /// Gets the value for the index specified.
  /// i = 0,1,2; j=0,1,2,3; k=0,1,...,26
  Real & operator()(unsigned int i, unsigned int j, unsigned int k);

  /**
   * Gets the value for the index specified.
   * used for const
   /// i = 0,1,2; j=0,1,2,3; k=0,1,...,26
   */
  Real operator()(unsigned int i, unsigned int j, unsigned int k) const;

  /// Zeros out the tensor.
  void zero();

  /// Print the rank three matrix
  void print(std::ostream & stm = Moose::out) const;

  /// copies values from a into this tensor
  RankThreeMatrixParakeet & operator= (const RankThreeMatrixParakeet & a);


protected:
  /// Dimensionality of rank-three matrix
  static const unsigned int N = LIBMESH_DIM;
  static const unsigned int NB = 27;

  /// The values of the rank-three matrix
  Real _vals[N][N+1][NB];


  template<class T>
  friend void dataStore(std::ostream &, T &, void *);

  template<class T>
  friend void dataLoad(std::istream &, T &, void *);
};

template<>
void dataStore(std::ostream &, RankThreeMatrixParakeet &, void *);

template<>
void dataLoad(std::istream &, RankThreeMatrixParakeet &, void *);

#endif //RANKTHREEMATRIXPARAKEET_H
