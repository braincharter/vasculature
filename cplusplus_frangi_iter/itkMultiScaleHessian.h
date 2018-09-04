#ifndef __itkMultiScaleHessian_h
#define __itkMultiScaleHessian_h

#include "itkVesselnessMeasurement.h"
#include "itkImageFileWriter.h"

#include "itkImageToImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"


//\class MultiScaleHessian
// \brief A filter to enhance structures using Hessian eigensystem-based
// measures in a multiscale framework
//
// The filter evaluates a Hessian-based enhancement measure, such as vesselness
// or objectness, at different scale levels. The Hessian-based measure is
// computed from the Hessian image at each scale level and the best response is
// selected.
//
// Minimum and maximum sigma value can be set using SetMinSigma and SetMaxSigma
// methods respectively.The number of scale levels is set using
// SetNumberOfSigmaSteps method.Exponentially distributed scale levels are
// computed within the bound set by the minimum and maximum sigma values
//
// The filter computes a second output image (accessed by the GetScalesOutput
// method) containing the scales at which each pixel gave the best response.
//
//  Manniesing, R, Viergever, MA, & Niessen, WJ (2006). Vessel Enhancing
//  Diffusion: A Scale Space Representation of Vessel Structures. Medical
//  Image Analysis, 10(6), 815-825./
//
//  Coded by :
//  http://www.insight-journal.org/browse/publication/163
//

template <typename TInputImage, typename THessianImage, typename TOutputImage>
class MultiScaleHessian
    : public itk::ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  typedef MultiScaleHessian Self;
  typedef itk::ImageToImageFilter<TInputImage, TOutputImage> Superclass;

  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef TInputImage InputImageType;
  typedef TOutputImage OutputImageType;
  typedef THessianImage HessianImageType;

  typedef VesselnessMeasurement<HessianImageType, OutputImageType>
      HessianToMeasureFilterType;

  typedef typename TInputImage::PixelType InputPixelType;
  typedef typename TOutputImage::PixelType OutputPixelType;
  typedef typename TOutputImage::RegionType OutputRegionType;

  static const unsigned int ImageDimension = InputImageType::ImageDimension;

  typedef float ScalesPixelType;
  typedef itk::Image<ScalesPixelType, ImageDimension> ScalesImageType;

  typedef itk::HessianRecursiveGaussianImageFilter<
      InputImageType, HessianImageType> HessianFilterType;

  typedef itk::Image<double, ImageDimension> UpdateBufferType;
  typedef typename UpdateBufferType::ValueType BufferValueType;

  typedef typename Superclass::DataObjectPointer DataObjectPointer;

  itkNewMacro(Self);

  itkTypeMacro(MultiScaleHessian, ImageToImageFilter);

  itkSetMacro(SigmaMinimum, double);
  itkGetConstMacro(SigmaMinimum, double);

  itkSetMacro(SigmaMaximum, double);
  itkGetConstMacro(SigmaMaximum, double);

  itkSetMacro(NumberOfSigmaSteps, unsigned int);
  itkGetConstMacro(NumberOfSigmaSteps, unsigned int);

  itkSetMacro(GenerateScalesOutput, bool);
  itkGetConstMacro(GenerateScalesOutput, bool);
  itkBooleanMacro(GenerateScalesOutput);

  itkSetMacro(GenerateHessianOutput, bool);
  itkGetConstMacro(GenerateHessianOutput, bool);
  itkBooleanMacro(GenerateHessianOutput);

  // Set/Get HessianToMeasureFilter. This will be a filter that takes
  // Hessian input image and produces enhanced output scalar image. The filter
  // must derive from itk::ImageToImage filter
  itkSetObjectMacro(HessianToMeasureFilter, HessianToMeasureFilterType);
  itkGetModifiableObjectMacro(HessianToMeasureFilter,
                              HessianToMeasureFilterType);

  // Methods to turn on/off flag to inform the filter that the Hessian-based
  // measure is non-negative (classical measures like Sato's and Frangi's are),
  // hence it has a minimum at zero. In this case, the update buffer is
  // initialized at zero, and the output scale and Hessian are zero in case the
  // Hessian-based measure returns zero for all scales. Otherwise, the minimum
  // output scale and Hessian are the ones obtained at scale SigmaMinimum. On by
  // default.
  itkSetMacro(NonNegativeHessianBasedMeasure, bool);
  itkGetConstMacro(NonNegativeHessianBasedMeasure, bool);
  itkBooleanMacro(NonNegativeHessianBasedMeasure);

  typedef enum
  {
    EquispacedSigmaSteps = 0,
    LogarithmicSigmaSteps = 1
  } SigmaStepMethodType;

  itkSetMacro(SigmaStepMethod, SigmaStepMethodType);
  itkGetConstMacro(SigmaStepMethod, SigmaStepMethodType);

  void SetSigmaStepMethodToEquispaced();
  void SetSigmaStepMethodToLogarithmic();

  const HessianImageType* GetHessianOutput() const;
  const ScalesImageType* GetScalesOutput() const;
  void EnlargeOutputRequestedRegion(itk::DataObject*);

  typedef itk::ProcessObject::DataObjectPointerArraySizeType
      DataObjectPointerArraySizeType;
  using Superclass::MakeOutput;
  virtual DataObjectPointer MakeOutput(DataObjectPointerArraySizeType idx);

  void SetBrightBlood(bool);
  bool GetBrightBlood();

  void SetFrangiOnly(bool);
  bool GetFrangiOnly();

  void SetScaleObject(bool);
  bool GetScaleObject();

  void SetAlpha(double);
  double GetAlpha();

  void SetBeta(double);
  double GetBeta();

  void SetC(double);
  double GetC();

protected:
  MultiScaleHessian();
  MultiScaleHessian(const bool nonNeg, double const& sigmaMin,
                    double const& sigmaMax, const unsigned int nbSigma,
                    const bool generateScale, const bool generateHessian);
  ~MultiScaleHessian() {}

  void PrintSelf(std::ostream& os, itk::Indent indent) const;
  void GenerateData(void);

private:
  void UpdateMaximumResponse(double sigma);

  double ComputeSigmaValue(int scaleLevel);

  void AllocateUpdateBuffer();

  MultiScaleHessian(const Self&);
  void operator=(const Self&);

  bool m_NonNegativeHessianBasedMeasure;

  double m_SigmaMinimum;
  double m_SigmaMaximum;

  unsigned int m_NumberOfSigmaSteps;
  SigmaStepMethodType m_SigmaStepMethod;

  typename HessianToMeasureFilterType::Pointer m_HessianToMeasureFilter;

  typename HessianFilterType::Pointer m_HessianFilter;

  typename UpdateBufferType::Pointer m_UpdateBuffer;

  bool m_GenerateScalesOutput;
  bool m_GenerateHessianOutput;
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultiScaleHessian.hxx"
#endif

#endif
