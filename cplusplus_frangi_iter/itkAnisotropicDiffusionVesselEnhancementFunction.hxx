#ifndef __itkAnisotropicDiffusionVesselEnhancementFunction_hxx
#define __itkAnisotropicDiffusionVesselEnhancementFunction_hxx

#include "itkAnisotropicDiffusionVesselEnhancementFunction.h"

#include "vnl/algo/vnl_symmetric_eigensystem.h"

template <class TImageType>
AnisotropicDiffusionVesselEnhancementFunction<
    TImageType>::AnisotropicDiffusionVesselEnhancementFunction()
{

  RadiusType r;
  r.Fill(1);

  this->SetRadius(r);

  NeighborhoodType itNeighbor;
  itNeighbor.SetRadius(r);

  m_Center = itNeighbor.Size() / 2;

  for (unsigned int i = 0; i < ImageDimension; ++i)
  {
    m_xStride[i] = itNeighbor.GetStride(i);
  }
}

template <class TImageType>
typename AnisotropicDiffusionVesselEnhancementFunction<TImageType>::TimeStepType
AnisotropicDiffusionVesselEnhancementFunction<
    TImageType>::ComputeGlobalTimeStep(void* derivateData) const
{
  // Returns the time step supplied by the user. We don't need
  // to use the global data supplied since we are returning a fixed value.
  return this->GetTimeStep();
}

template <class TImageType>
typename AnisotropicDiffusionVesselEnhancementFunction<TImageType>::PixelType
AnisotropicDiffusionVesselEnhancementFunction<TImageType>::ComputeUpdate(
    const NeighborhoodType& itNeighbor, void* derivateData,
    const FloatOffsetType& offset)
{
  return static_cast<PixelType>(0.0);
}

// =============================================================================
// Compute the diffusion equation to optimize the anisotropy of the tensor
// inside the vessel.
// m_dx -> Intensity first derivative
// m_dxy -> Intensity second derivative
// m_DTdxy -> Diffusion tensor first derivative
// =============================================================================
template <class TImageType>
typename AnisotropicDiffusionVesselEnhancementFunction<TImageType>::PixelType
AnisotropicDiffusionVesselEnhancementFunction<TImageType>::ComputeUpdate(
    const NeighborhoodType& itNeighbor,
    const DiffusionTensorNeighborhoodType& itTensorNeighbor,
    DerivativeStructType* dd)
{
  const double ZERO = itk::NumericTraits<double>::Zero;
  const double centerValue = itNeighbor.GetCenterPixel();

  dd->m_GradMagSqr = 1.0e-6;

  for (unsigned int i = 0; i < ImageDimension; ++i)
  {
    const auto positionA = m_Center + m_xStride[i];
    const auto positionB = m_Center - m_xStride[i];

    dd->m_dx[i] =
        (itNeighbor.GetPixel(positionA) - itNeighbor.GetPixel(positionB)) / 2.0;

    dd->m_dxy[i][i] = itNeighbor.GetPixel(positionA) +
                      itNeighbor.GetPixel(positionB) - 2.0 * centerValue;

    for (unsigned int j = i + 1; j < ImageDimension; ++j)
    {
      const auto positionAa = positionB - m_xStride[j];
      const auto positionBa = positionB + m_xStride[j];
      const auto positionCa = positionA - m_xStride[j];
      const auto positionDa = positionA + m_xStride[j];

      dd->m_dxy[i][j] = dd->m_dxy[j][i] =
          (itNeighbor.GetPixel(positionAa) - itNeighbor.GetPixel(positionBa) -
           itNeighbor.GetPixel(positionCa) + itNeighbor.GetPixel(positionDa)) /
          4.0;
    }

    const TensorPixelType positionATensorValue =
        itTensorNeighbor.GetPixel(positionA);
    const TensorPixelType positionBTensorValue =
        itTensorNeighbor.GetPixel(positionB);

    for (unsigned int j = 0; j < ImageDimension; ++j)
    {
      dd->m_DTdxy[i][j] =
          (positionATensorValue(i, j) - positionBTensorValue(i, j)) / 2.0;
    }
  }

  // Compute the diffusion tensor matrix for the first derivatives
  const TensorPixelType centerTensorValue = itTensorNeighbor.GetCenterPixel();

  const double pdWrtDiffusion1 = dd->m_DTdxy[0][0] * dd->m_dx[0] +
                                 dd->m_DTdxy[0][1] * dd->m_dx[1] +
                                 dd->m_DTdxy[0][2] * dd->m_dx[2];

  const double pdWrtDiffusion2 = dd->m_DTdxy[1][0] * dd->m_dx[0] +
                                 dd->m_DTdxy[1][1] * dd->m_dx[1] +
                                 dd->m_DTdxy[1][2] * dd->m_dx[2];

  const double pdWrtDiffusion3 = dd->m_DTdxy[2][0] * dd->m_dx[0] +
                                 dd->m_DTdxy[2][1] * dd->m_dx[1] +
                                 dd->m_DTdxy[2][2] * dd->m_dx[2];

  const double pdWrtImageIntensity1 =
      centerTensorValue(0, 0) * dd->m_dxy[0][0] +
      centerTensorValue(0, 1) * dd->m_dxy[0][1] +
      centerTensorValue(0, 2) * dd->m_dxy[0][2];

  const double pdWrtImageIntensity2 =
      centerTensorValue(1, 0) * dd->m_dxy[1][0] +
      centerTensorValue(1, 1) * dd->m_dxy[1][1] +
      centerTensorValue(1, 2) * dd->m_dxy[1][2];

  const double pdWrtImageIntensity3 =
      centerTensorValue(2, 0) * dd->m_dxy[2][0] +
      centerTensorValue(2, 1) * dd->m_dxy[2][1] +
      centerTensorValue(2, 2) * dd->m_dxy[2][2];

  const double total = pdWrtDiffusion1 + pdWrtDiffusion2 + pdWrtDiffusion3 +
                       pdWrtImageIntensity1 + pdWrtImageIntensity2 +
                       pdWrtImageIntensity3;

  return static_cast<PixelType>(total);
}

#endif
