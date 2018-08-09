
#ifndef __itkDerivativeStruct_h
#define __itkDerivativeStruct_h

#include "itkFiniteDifferenceFunction.h"
#include "vnl/vnl_matrix_fixed.h"

// A derivative data type for this class of equations. Used to store
// values that are needed in calculating the time step and other intermediate
// products such as derivatives that may be used by virtual functions called
// from ComputeUpdate. Caching these values here allows the ComputeUpdate
// function to be const and thread safe.
template <class TImageType> struct DerivativeStruct
{
  static const unsigned int ImageDimension =
      itk::FiniteDifferenceFunction<TImageType>::ImageDimension;

  // Hessian matrix
  vnl_matrix_fixed<double, ImageDimension, ImageDimension> m_dxy;

  // diffusion tensor first derivative matrix
  vnl_matrix_fixed<double, ImageDimension, ImageDimension> m_DTdxy;

  // Array of first derivatives
  double m_dx[ImageDimension];

  double m_GradMagSqr;
};

#endif