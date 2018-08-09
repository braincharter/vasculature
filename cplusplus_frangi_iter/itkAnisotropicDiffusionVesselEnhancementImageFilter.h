#ifndef __itkAnisotropicDiffusionVesselEnhancementImageFilter_h
#define __itkAnisotropicDiffusionVesselEnhancementImageFilter_h

#include "itkAnisotropicDiffusionVesselEnhancementFunction.h"

#include "itkDerivativeStruct.h"
#include "itkMultiScaleHessian.h"
#include "itkSymmetricEigenVectorAnalysisImageFilter.h"

#include "itkDiffusionTensor3D.h"
#include "itkFiniteDifferenceImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkMultiThreader.h"


// \class AnisotropicDiffusionVesselEnhancementImageFilter
// \brief This class create the Tensor D of Manniesing et al.
//
//  Manniesing, R, Viergever, MA, & Niessen, WJ (2006). Vessel Enhancing
//  Diffusion: A Scale Space Representation of Vessel Structures. Medical
//  Image Analysis, 10(6), 815-825./
//
//  Coded by :
//  http://www.insight-journal.org/browse/publication/163

template <class TInputImage, class TOutputImage>
class AnisotropicDiffusionVesselEnhancementImageFilter
    : public itk::FiniteDifferenceImageFilter<TInputImage, TOutputImage>
{
public:
  typedef AnisotropicDiffusionVesselEnhancementImageFilter Self;

  typedef itk::FiniteDifferenceImageFilter<TInputImage, TOutputImage>
      Superclass;

  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;
  itkNewMacro(Self);
  itkTypeMacro(AnisotropicDiffusionVesselEnhancementImageFilter,
               ImageToImageFilter);

  typedef DerivativeStruct<TInputImage> DerivativeStructType;

  typedef typename Superclass::InputImageType InputImageType;
  typedef typename Superclass::OutputImageType OutputImageType;
  typedef typename Superclass::PixelType PixelType;

  typedef itk::Image<itk::DiffusionTensor3D<double>, 3>
      DiffusionTensorImageType;

  // Dimensionality of input and output data is assumed to be the same.
  // It is inherited from the superclass.
  static const unsigned int ImageDimension = Superclass::ImageDimension;

  typedef itk::SymmetricSecondRankTensor<double, ImageDimension>
      TensorPixelType;
  typedef itk::Image<TensorPixelType, ImageDimension> TensorImageType;

  typedef AnisotropicDiffusionVesselEnhancementFunction<InputImageType>
      FiniteDifferenceFunctionType;

  typedef itk::HessianRecursiveGaussianImageFilter<
      InputImageType, TensorImageType> HessianFilterType;

  typedef itk::Image<double, 3> VesselnessOutputImageType;

  typedef MultiScaleHessian<InputImageType, TensorImageType,
                            VesselnessOutputImageType>
      MultiScaleVesselnessFilterType;

  typedef TensorImageType HessianImageType;

  typedef float ScalesPixelType;
  typedef itk::Image<ScalesPixelType, ImageDimension> ScalesImageType;

  typedef itk::Matrix<double, ImageDimension, ImageDimension> MatrixType;

  typedef itk::Image<MatrixType, ImageDimension> OutputMatrixImageType;

  typedef itk::FixedArray<double, ImageDimension> EigenValueArrayType;

  typedef itk::Image<EigenValueArrayType, ImageDimension>
      EigenAnalysisOutputImageType;

  typedef SymmetricEigenVectorAnalysisImageFilter<
      TensorImageType, EigenAnalysisOutputImageType, OutputMatrixImageType>
      EigenVectorMatrixAnalysisFilterType;

  typedef typename Superclass::TimeStepType TimeStepType;

  typedef OutputImageType UpdateBufferType;

  typedef typename FiniteDifferenceFunctionType::DiffusionTensorNeighborhoodType
      DiffusionTensorNeighborhoodType;

  itkSetMacro(TimeStep, double);
  itkSetMacro(Epsilon, double);
  itkSetMacro(WStrength, double);
  itkSetMacro(Sensitivity, double);

  itkGetMacro(TimeStep, double);
  itkGetMacro(Epsilon, double);
  itkGetMacro(WStrength, double);
  itkGetMacro(Sensitivity, double);

  itkSetMacro(NumberOfIterations, unsigned int);
  itkGetMacro(NumberOfIterations, unsigned int);

  itkSetMacro(GenerateIterationFiles, bool);
  itkGetMacro(GenerateIterationFiles, bool);

#ifdef ITK_USE_CONCEPT_CHECKING
  itkConceptMacro(OutputTimesDoubleCheck,
                  (itk::Concept::MultiplyOperator<PixelType, double>));
  itkConceptMacro(OutputAdditiveOperatorsCheck,
                  (itk::Concept::AdditiveOperators<PixelType>));
  itkConceptMacro(
      InputConvertibleToOutputCheck,
      (itk::Concept::Convertible<typename TInputImage::PixelType, PixelType>));
#endif

  void SetSigmaMin(double);
  void SetSigmaMax(double);
  void SetNumberOfSigmaSteps(int);
  void SetBrightBlood(bool);
  void SetFrangiOnly(bool);

  void SetAlpha(double);
  void SetBeta(double);
  void SetC(double);
  void SetScaleObject(bool);
  void SetGenerateScale(bool);
  void SetGenerateHessian(bool);

  double GetSigmaMin();
  double GetSigmaMax();
  int GetNumberOfSigmaSteps();
  bool GetBrightBlood();
  bool GetFrangiOnly();

  double GetAlpha();
  double GetBeta();
  double GetC();
  bool GetScaleObject();
  bool GetGenerateScale();
  bool GetGenerateHessian();

  const HessianImageType* GetHessianOutput() const;
  const ScalesImageType* GetScalesOutput() const;

protected:
  AnisotropicDiffusionVesselEnhancementImageFilter();
  AnisotropicDiffusionVesselEnhancementImageFilter(
      double const& timeStep, const unsigned int nbIteration,
      double const& wStrength, double const& sensitivity, double const& epsilon,
      const bool generateIterationFiles);

  ~AnisotropicDiffusionVesselEnhancementImageFilter() {}

  virtual void GenerateData();

  // A simple method to copy the data from the input to the output. ( Supports
  // "read-only" image adaptors in the case where the input image type converts
  // to a different output image type. )
  virtual void CopyInputToOutput();

  // This method applies changes from the m_UpdateBuffer to the output using
  // the ThreadedApplyUpdate() method and a multithreading mechanism.  "dt" is
  // the time step to use for the update of each pixel.
  virtual void ApplyUpdate(const TimeStepType& dt);

  // Method to allow subclasses to get direct access to the update
  // buffer
  virtual typename UpdateBufferType::Pointer GetUpdateBuffer() { return m_UpdateBuffer; }

  // This method populates an update buffer with changes for each pixel in the
  // output using the ThreadedCalculateChange() method and a multithreading
  // mechanism. Returns value is a time step to be used for the update.
  virtual TimeStepType CalculateChange();

  // This method allocates storage in m_UpdateBuffer.  It is called from
  // Superclass::GenerateData().
  virtual void AllocateUpdateBuffer();

  void AllocateDiffusionTensorImage();
  void UpdateDiffusionTensorImage();

  typedef typename UpdateBufferType::RegionType ThreadRegionType;
  typedef typename DiffusionTensorImageType::RegionType
      ThreadDiffusionImageRegionType;

  //  Does the actual work of updating the output from the UpdateContainer
  //   over an output region supplied by the multithreading mechanism.
  //  \sa ApplyUpdate
  // \sa ApplyUpdateThreaderCallback
  virtual void ThreadedApplyUpdate(
      TimeStepType dt, const ThreadRegionType& regionToProcess,
      const ThreadDiffusionImageRegionType& diffusionRegionToProcess,
      int threadId);

  // Does the actual work of calculating change over a region supplied by
  // the multithreading mechanism.
  // \sa CalculateChange
  // \sa CalculateChangeThreaderCallback
  virtual TimeStepType ThreadedCalculateChange(
      const ThreadRegionType& regionToProcess,
      const ThreadDiffusionImageRegionType& diffusionRegionToProcess,
      int threadId);

  virtual void InitializeIteration();

private:
  // purposely not implemented
  AnisotropicDiffusionVesselEnhancementImageFilter(const Self&);
  void operator=(const Self&); // purposely not implemented

  // Structure for passing information into static callback methods.  Used in
  // the subclasses' threading mechanisms.
  struct DenseFDThreadStruct
  {
    AnisotropicDiffusionVesselEnhancementImageFilter* Filter;
    TimeStepType TimeStep;
    std::vector<TimeStepType> TimeStepList;
    std::vector<bool> ValidTimeStepList;
  };

  // This callback method uses ImageSource::SplitRequestedRegion to acquire an
  // output region that it passes to ThreadedApplyUpdate for processing.
  static ITK_THREAD_RETURN_TYPE ApplyUpdateThreaderCallback(void* arg);

  // This callback method uses SplitUpdateContainer to acquire a region
  // which it then passes to ThreadedCalculateChange for processing.
  static ITK_THREAD_RETURN_TYPE CalculateChangeThreaderCallback(void* arg);

  // The buffer that holds the updates for an iteration of the algorithm.
  typename UpdateBufferType::Pointer m_UpdateBuffer;

  TimeStepType m_TimeStep;
  typename DiffusionTensorImageType::Pointer m_DiffusionTensorImage;
  typename MultiScaleVesselnessFilterType::Pointer m_MultiScaleVesselnessFilter;
  typename HessianFilterType::Pointer m_HessianFilter;

  typename EigenVectorMatrixAnalysisFilterType::Pointer
      m_EigenVectorMatrixAnalysisFilter;

  double m_Epsilon;
  double m_WStrength;
  double m_Sensitivity;

  unsigned int m_NumberOfIterations;

  bool m_GenerateIterationFiles;
};

#if ITK_TEMPLATE_TXX
#include "itkAnisotropicDiffusionVesselEnhancementImageFilter.hxx"
#endif

#endif
