#ifndef __itkAnisotropicDiffusionVesselEnhancementFunction_h
#define __itkAnisotropicDiffusionVesselEnhancementFunction_h

#include "itkDerivativeStruct.h"
#include "itkSymmetricSecondRankTensor.h"

#include "itkDiffusionTensor3D.h"
#include "itkFiniteDifferenceFunction.h"
#include "vnl/vnl_matrix_fixed.h"


// \class AnisotropicDiffusionVesselEnhancementFunction
// \brief This class iteratively enhances vessels in an image by solving
// non-linear diffusion equation developed by Manniesing et al.
//
//  Manniesing, R, Viergever, MA, & Niessen, WJ (2006). Vessel Enhancing
//  Diffusion: A Scale Space Representation of Vessel Structures. Medical
//  Image Analysis, 10(6), 815-825./
//
//  Coded by :
//  http://www.insight-journal.org/browse/publication/163
//
template <class TImageType>
class AnisotropicDiffusionVesselEnhancementFunction
    : public itk::FiniteDifferenceFunction<TImageType>
{
public:
  typedef AnisotropicDiffusionVesselEnhancementFunction Self;
  typedef itk::FiniteDifferenceFunction<TImageType> Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  itkNewMacro(Self);

  itkTypeMacro(AnisotropicDiffusionVesselEnhancementFunction,
               itk::FiniteDifferenceFunction);

  static const unsigned int ImageDimension = Superclass::ImageDimension;

  typedef double TimeStepType;
  typedef typename Superclass::ImageType ImageType;
  typedef typename Superclass::PixelType PixelType;
  typedef typename Superclass::RadiusType RadiusType;
  typedef typename Superclass::NeighborhoodType NeighborhoodType;
  typedef typename Superclass::FloatOffsetType FloatOffsetType;

  typedef DerivativeStruct<TImageType> DerivativeStructType;

  typedef itk::Image<itk::DiffusionTensor3D<double>, 3>
      DiffusionTensorImageType;

  // The default boundary condition for finite difference
  // functions that is used unless overridden in the Evaluate() method.
  typedef itk::ZeroFluxNeumannBoundaryCondition<DiffusionTensorImageType>
      DefaultBoundaryConditionType;

  typedef itk::ConstNeighborhoodIterator<DiffusionTensorImageType,
                                         DefaultBoundaryConditionType>
      DiffusionTensorNeighborhoodType;

  typedef itk::SymmetricSecondRankTensor<double> TensorPixelType;

  virtual PixelType
  ComputeUpdate(const NeighborhoodType& neighborhood, void* derivateData,
                const FloatOffsetType& offset = FloatOffsetType(0.0)) override;

  PixelType
  ComputeUpdate(const NeighborhoodType& neighborhood,
                const DiffusionTensorNeighborhoodType& neighborhoodTensor,
                DerivativeStructType* derivateData);

  // Compute the time step for an update given a derivatie data structure.
  virtual TimeStepType ComputeGlobalTimeStep(void* derivateData) const;

  // Returns a pointer to a global data structure that is passed to this
  // object from the solver at each calculation.
  virtual void* GetGlobalDataPointer() const
  {
    return new DerivativeStructType();
  }

  virtual void ReleaseGlobalDataPointer(void* derivateData) const
  {
    delete static_cast<DerivativeStructType*>(derivateData);
  }

  // Set/Get the time step. For this class of anisotropic diffusion filters,
  // the time-step is supplied by the user and remains fixed for all updates.
  void SetTimeStep(const TimeStepType& t) { m_TimeStep = t; }

  const TimeStepType& GetTimeStep() const { return m_TimeStep; }

protected:
  AnisotropicDiffusionVesselEnhancementFunction();

  virtual ~AnisotropicDiffusionVesselEnhancementFunction() {}

  std::slice x_slice[ImageDimension];

  // The offset of the center pixel in the neighborhood.
  unsigned int m_Center;

  // Stride length along the y-dimension.
  unsigned int m_xStride[ImageDimension];

private:
  // purposely not implemented
  AnisotropicDiffusionVesselEnhancementFunction(const Self&);
  void operator=(const Self&); // purposely not implemented

  TimeStepType m_TimeStep;
};

#if ITK_TEMPLATE_TXX
#include "itkAnisotropicDiffusionVesselEnhancementFunction.hxx"
#endif

#endif