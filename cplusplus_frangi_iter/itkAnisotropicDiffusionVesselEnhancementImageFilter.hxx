#ifndef __itkAnisotropicDiffusionVesselEnhancementImageFilter_hxx_
#define __itkAnisotropicDiffusionVesselEnhancementImageFilter_hxx_

#include "itkAnisotropicDiffusionVesselEnhancementImageFilter.h"

#include "itkAnisotropicDiffusionVesselEnhancementFunction.h"
#include "itkDerivativeStruct.h"

#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkNumericTraits.h"
#include "itkVector.h"

#include <list>

template <class TInputImage, class TOutputImage>
AnisotropicDiffusionVesselEnhancementImageFilter<TInputImage, TOutputImage>::
    AnisotropicDiffusionVesselEnhancementImageFilter()
    : AnisotropicDiffusionVesselEnhancementImageFilter<
          TInputImage,
          TOutputImage>::AnisotropicDiffusionVesselEnhancementImageFilter(10e-3,
                                                                          2,
                                                                          25.0,
                                                                          5.0,
                                                                          10e-2,
                                                                          false)
{
} // nbIteration is set to 2 by default, to apply 1 VED iteration and after to 
  // apply a last frangi. (it itkVEDMain, there is also a +1 for that)

template <class TInputImage, class TOutputImage>
AnisotropicDiffusionVesselEnhancementImageFilter<TInputImage, TOutputImage>::
    AnisotropicDiffusionVesselEnhancementImageFilter(
        double const& timeStep, const unsigned int nbIteration,
        double const& wStrength, double const& sensitivity,
        double const& epsilon, const bool generateIterationFiles)
    : m_TimeStep{timeStep}, m_NumberOfIterations{nbIteration},
      m_WStrength{wStrength}, m_Sensitivity{sensitivity}, m_Epsilon{epsilon},
      m_GenerateIterationFiles{generateIterationFiles}
{
  m_UpdateBuffer = UpdateBufferType::New();
  m_DiffusionTensorImage = DiffusionTensorImageType::New();

  this->SetNumberOfIterations(m_NumberOfIterations);


  typename AnisotropicDiffusionVesselEnhancementFunction<
      UpdateBufferType>::Pointer q =
      AnisotropicDiffusionVesselEnhancementFunction<UpdateBufferType>::New();
  this->SetDifferenceFunction(q);

  m_HessianFilter = HessianFilterType::New();
  m_EigenVectorMatrixAnalysisFilter =
      EigenVectorMatrixAnalysisFilterType::New();
  m_EigenVectorMatrixAnalysisFilter->SetDimension(TensorPixelType::Dimension);

  m_MultiScaleVesselnessFilter = MultiScaleVesselnessFilterType::New();
}

template <class TInputImage, class TOutputImage>
void AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::InitializeIteration()
{
  itkDebugMacro(<< "InitializeIteration() called ");

  auto* f = dynamic_cast<
      AnisotropicDiffusionVesselEnhancementFunction<UpdateBufferType>*>(
      this->GetDifferenceFunction().GetPointer());

  if (!f)
  {
    throw itk::ExceptionObject(
        __FILE__, __LINE__,
        "Anisotropic diffusion Vessel Enhancement function is not set.",
        ITK_LOCATION);
  }

  f->SetTimeStep(m_TimeStep);

  // Check the timestep for stability.
  double minSpacing = 1.0;
  if (this->GetUseImageSpacing())
  {
    minSpacing = *(std::min_element(this->GetInput()->GetSpacing().Begin(),
                                    this->GetInput()->GetSpacing().End()));
  }

  const double ratio =
      minSpacing / vcl_pow(2.0, static_cast<double>(ImageDimension) + 1);

  if (m_TimeStep > ratio)
  {
    itkWarningMacro(<< std::endl
                    << "Anisotropic diffusion unstable time step:" << m_TimeStep
                    << std::endl
                    << "Minimum stable time step for this image is " << ratio);
  }

  f->InitializeIteration();

  if (this->GetNumberOfIterations() != 0)
  {
    this->UpdateProgress(static_cast<float>(this->GetElapsedIterations()) /
                         static_cast<float>(this->GetNumberOfIterations()));
  }
  else
  {
    this->UpdateProgress(0);
  }


  this->UpdateDiffusionTensorImage();
}

template <class TInputImage, class TOutputImage>
void AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::SetSigmaMin(double value)
{
  m_MultiScaleVesselnessFilter->SetSigmaMinimum(value);
  this->Modified();
}

template <class TInputImage, class TOutputImage>
void AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::SetSigmaMax(double value)
{
  m_MultiScaleVesselnessFilter->SetSigmaMaximum(value);
  this->Modified();
}

template <class TInputImage, class TOutputImage>
void AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::SetNumberOfSigmaSteps(int value)
{
  m_MultiScaleVesselnessFilter->SetNumberOfSigmaSteps(value);
  this->Modified();
}

template <class TInputImage, class TOutputImage>
void AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::SetBrightBlood(bool value)
{
  m_MultiScaleVesselnessFilter->SetBrightBlood(value);
  this->Modified();
}

template <class TInputImage, class TOutputImage>
void AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::SetFrangiOnly(bool value)
{
  m_MultiScaleVesselnessFilter->SetFrangiOnly(value);
  this->Modified();
}

template <class TInputImage, class TOutputImage>
void AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::SetAlpha(double value)
{
  m_MultiScaleVesselnessFilter->SetAlpha(value);
  this->Modified();
}

template <class TInputImage, class TOutputImage>
void AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::SetBeta(double value)
{
  m_MultiScaleVesselnessFilter->SetBeta(value);
  this->Modified();
}

template <class TInputImage, class TOutputImage>
void AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::SetC(double value)
{
  m_MultiScaleVesselnessFilter->SetC(value);
  this->Modified();
}

template <class TInputImage, class TOutputImage>
void AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::SetScaleObject(bool value)
{
  m_MultiScaleVesselnessFilter->SetScaleObject(value);
  this->Modified();
}

template <class TInputImage, class TOutputImage>
void AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::SetGenerateScale(bool value)
{
  m_MultiScaleVesselnessFilter->SetGenerateScalesOutput(value);
  this->Modified();
}

template <class TInputImage, class TOutputImage>
void AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::SetGenerateHessian(bool value)
{
  m_MultiScaleVesselnessFilter->SetGenerateHessianOutput(value);
  this->Modified();
}

template <class TInputImage, class TOutputImage>
double
AnisotropicDiffusionVesselEnhancementImageFilter<TInputImage,
                                                 TOutputImage>::GetSigmaMin()
{
  return m_MultiScaleVesselnessFilter->GetSigmaMinimum();
}

template <class TInputImage, class TOutputImage>
double
AnisotropicDiffusionVesselEnhancementImageFilter<TInputImage,
                                                 TOutputImage>::GetSigmaMax()
{
  return m_MultiScaleVesselnessFilter->GetSigmaMaximum();
}

template <class TInputImage, class TOutputImage>
int AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::GetNumberOfSigmaSteps()
{
  return m_MultiScaleVesselnessFilter->GetNumberOfSigmaSteps();
}

template <class TInputImage, class TOutputImage>
bool AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::GetBrightBlood()
{
  return m_MultiScaleVesselnessFilter->GetBrightBlood();
}

template <class TInputImage, class TOutputImage>
bool AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::GetFrangiOnly()
{
  return m_MultiScaleVesselnessFilter->GetFrangiOnly();
}

template <class TInputImage, class TOutputImage>
double
AnisotropicDiffusionVesselEnhancementImageFilter<TInputImage,
                                                 TOutputImage>::GetAlpha()
{
  return m_MultiScaleVesselnessFilter->GetAlpha();
}

template <class TInputImage, class TOutputImage>
double AnisotropicDiffusionVesselEnhancementImageFilter<TInputImage,
                                                        TOutputImage>::GetBeta()
{
  return m_MultiScaleVesselnessFilter->GetBeta();
}

template <class TInputImage, class TOutputImage>
double AnisotropicDiffusionVesselEnhancementImageFilter<TInputImage,
                                                        TOutputImage>::GetC()
{
  return m_MultiScaleVesselnessFilter->GetC();
}

template <class TInputImage, class TOutputImage>
bool AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::GetScaleObject()
{
  return m_MultiScaleVesselnessFilter->GetScaleObject();
}

template <class TInputImage, class TOutputImage>
bool AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::GetGenerateScale()
{
  return m_MultiScaleVesselnessFilter->GetGenerateScalesOutput();
}

template <class TInputImage, class TOutputImage>
bool AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::GetGenerateHessian()
{
  return m_MultiScaleVesselnessFilter->GetGenerateHessianOutput();
}

// Get the image containing the Hessian at which each pixel gave the best
// response
template <class TInputImage, class TOutputImage>
const typename AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::HessianImageType*
AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::GetHessianOutput() const
{
  return m_MultiScaleVesselnessFilter->GetHessianOutput();
}

// Get the image containing the scales at which each pixel gave the best
// response
template <class TInputImage, class TOutputImage>
const typename AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::ScalesImageType*
AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::GetScalesOutput() const
{

  return m_MultiScaleVesselnessFilter->GetScalesOutput();
}

template <class TInputImage, class TOutputImage>
void AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::CopyInputToOutput()
{
  typename TInputImage::ConstPointer input = this->GetInput();
  typename TOutputImage::Pointer output = this->GetOutput();

  if (!input || !output)
  {
    itkExceptionMacro(<< "Either input or output is NULL.");
  }

  // Check if we are doing in-place filtering.
  if (this->GetInPlace() && (typeid(TInputImage) == typeid(TOutputImage)))
  {
    typename TInputImage::Pointer const tempPtr =
        dynamic_cast<TInputImage*>(output.GetPointer());
    if (tempPtr.IsNotNull() &&
        tempPtr->GetPixelContainer() == input->GetPixelContainer())
    {
      return;
    }
  }

  itk::ImageRegionConstIterator<TInputImage> itInput(
      input, output->GetRequestedRegion());
  itk::ImageRegionIterator<TOutputImage> itOuput(output,
                                                 output->GetRequestedRegion());

  while (!itOuput.IsAtEnd())
  {
    itOuput.Value() = static_cast<PixelType>(itInput.Get());
    ++itInput;
    ++itOuput;
  }
}

template <class TInputImage, class TOutputImage>
void AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::AllocateUpdateBuffer()
{
  // The update buffer looks just like the output and holds the change in
  // the pixel.
  typename TOutputImage::Pointer output = this->GetOutput();

  m_UpdateBuffer->SetSpacing(output->GetSpacing());
  m_UpdateBuffer->SetOrigin(output->GetOrigin());
  m_UpdateBuffer->SetLargestPossibleRegion(output->GetLargestPossibleRegion());
  m_UpdateBuffer->SetRequestedRegion(output->GetRequestedRegion());
  m_UpdateBuffer->SetBufferedRegion(output->GetBufferedRegion());
  m_UpdateBuffer->Allocate();
}

template <class TInputImage, class TOutputImage>
void AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::AllocateDiffusionTensorImage()
{
  itkDebugMacro(<< "AllocateDiffusionTensorImage() called");

  // The diffusionTensor image has the same size as the output and holds
  // the diffusion tensor matrix for each pixel.
  typename TOutputImage::Pointer output = this->GetOutput();

  m_DiffusionTensorImage->SetSpacing(output->GetSpacing());
  m_DiffusionTensorImage->SetOrigin(output->GetOrigin());
  m_DiffusionTensorImage->SetLargestPossibleRegion(
      output->GetLargestPossibleRegion());
  m_DiffusionTensorImage->SetRequestedRegion(output->GetRequestedRegion());
  m_DiffusionTensorImage->SetBufferedRegion(output->GetBufferedRegion());
  m_DiffusionTensorImage->Allocate();
}

// =============================================================================
// Create the tensor D with the vesselness measure and eigen value/vector from
// the Hessian matrix. Next, the tensor is inserted in the diffusion equation.
// =============================================================================
template <class TInputImage, class TOutputImage>
void AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::UpdateDiffusionTensorImage()
{
  itkDebugMacro(<< "UpdateDiffusionTensorImage() called");

  m_MultiScaleVesselnessFilter->SetInput(this->GetOutput());
  m_MultiScaleVesselnessFilter->Modified();
  m_MultiScaleVesselnessFilter->Update();

  if (this->GetFrangiOnly())
  {
    std::cout << "Frangi vesselness measure has been computed and will be "
                 "saved. Program will exit now based on user argument.";
    return;
  }


  if(this->GetNumberOfIterations() - 1 == this->GetElapsedIterations())
  {
    std::cout << "One more iteration has been added to apply Frangi equation"
                 "on the VED smoothing and save it as output.\n";

    //this->SetOutput(0, m_MultiScaleVesselnessFilter->GetOutput());

    itk::ImageRegionIterator<InputImageType> 
        itUpdate(m_MultiScaleVesselnessFilter->GetOutput(),
        m_MultiScaleVesselnessFilter->GetOutput()->GetRequestedRegion());

    itk::ImageRegionIterator<OutputImageType> 
        itOutput(this->GetOutput(), this->GetOutput()->GetRequestedRegion());

    itUpdate = itUpdate.Begin();
    itOutput = itOutput.Begin();

    while (!itUpdate.IsAtEnd())
    {
      itOutput.Value() = itUpdate.Value();

      ++itOutput;
      ++itUpdate;
    }

    return;
  }

  m_EigenVectorMatrixAnalysisFilter->SetInput(
      m_MultiScaleVesselnessFilter->GetHessianOutput());

  m_EigenVectorMatrixAnalysisFilter->Update();

  typename OutputMatrixImageType::Pointer eigenVectorMatrixOutputImage =
      m_EigenVectorMatrixAnalysisFilter->GetOutput();

  typedef itk::ImageRegionIterator<OutputMatrixImageType>
      EigenVectorMatrixIteratorType;

  EigenVectorMatrixIteratorType itEigen(
      eigenVectorMatrixOutputImage,
      eigenVectorMatrixOutputImage->GetLargestPossibleRegion());
  itEigen.GoToBegin();

  // Obtain the vessleness measure from Frangi equation.
  typename MultiScaleVesselnessFilterType::OutputImageType::Pointer
      MultiScaleHessianOutputImage;
  MultiScaleHessianOutputImage = m_MultiScaleVesselnessFilter->GetOutput();

  typedef typename MultiScaleVesselnessFilterType::OutputImageType
      MultiScaleHessianOutputImageType;
  typedef itk::ImageRegionIterator<MultiScaleHessianOutputImageType>
      MultiScaleHessianIteratorType;
  MultiScaleHessianIteratorType itHessian(
      MultiScaleHessianOutputImage,
      MultiScaleHessianOutputImage->GetLargestPossibleRegion());
  itHessian.GoToBegin();

  const double iS = 1.0 / m_Sensitivity;

  typedef itk::ImageRegionIterator<DiffusionTensorImageType>
      DiffusionTensorIteratorType;

  DiffusionTensorIteratorType itDiff(
      m_DiffusionTensorImage,
      m_DiffusionTensorImage->GetLargestPossibleRegion());
  itDiff.GoToBegin();

  std::cout << "(In UpdateDiffusionTensorImage) Compute D tensor. \n";
    
    
  while (!itDiff.IsAtEnd())
  {

    typename DiffusionTensorImageType::PixelType tensor;
   
    // Generate matrix "Q" with the eigenvectors of the Hessian matrix.
    const MatrixType hessianEigenVectorMatrix = itEigen.Get();
    const MatrixType hessianEigenVectorMatrixTranspose(
        hessianEigenVectorMatrix.GetTranspose());

    // Generate the diagonal matrix with the eigen values.
    MatrixType eigenValueMatrix;
    eigenValueMatrix.SetIdentity();

    const double vesselnessValue = static_cast<double>(itHessian.Get());
    const double powedVesselness = vcl_pow(vesselnessValue, iS);

    const double lambda1 = 1 + (m_WStrength - 1) * powedVesselness;
    const double lambda2 = 1 + (m_Epsilon - 1) * powedVesselness;

    // lambda3 = lambda2, no needs to create lamdba3.
    eigenValueMatrix(0, 0) = lambda1;
    eigenValueMatrix(1, 1) = lambda2;
    eigenValueMatrix(2, 2) = lambda2;

    const MatrixType productMatrix = hessianEigenVectorMatrix *
                                     eigenValueMatrix *
                                     hessianEigenVectorMatrixTranspose;

    for (unsigned int i = 0; i < ImageDimension; ++i)
    {
      for (unsigned int j = 0; j < ImageDimension; ++j)
      {
        tensor(i, j) = productMatrix(i, j);
      }

    }

    itDiff.Set(tensor);
      

    ++itDiff;
    ++itEigen;
    ++itHessian;
  }
}

template <class TInputImage, class TOutputImage>
void AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::ApplyUpdate(const TimeStepType& dt)
{

  itkDebugMacro(<< "ApplyUpdate Invoked with time step size: " << dt);

  DenseFDThreadStruct str;
  str.Filter = this;
  str.TimeStep = dt;
  this->GetMultiThreader()->SetNumberOfThreads(this->GetNumberOfThreads());
  this->GetMultiThreader()->SetSingleMethod(this->ApplyUpdateThreaderCallback,
                                            &str);

  this->GetMultiThreader()->SingleMethodExecute();
}

template <class TInputImage, class TOutputImage>
ITK_THREAD_RETURN_TYPE AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::ApplyUpdateThreaderCallback(void* arg)
{

  const auto threadInfo =
      static_cast<itk::MultiThreader::ThreadInfoStruct*>(arg);
  const int threadId = threadInfo->ThreadID;
  const int threadCount = threadInfo->NumberOfThreads;
  const auto str = static_cast<DenseFDThreadStruct*>(threadInfo->UserData);

  // Execute the actual method with appropriate output region
  // to find out in how many pieces extent can be split into.
  // Using the SplitRequestedRegion method from itk::ImageSource.
  ThreadRegionType splitRegion;
  str->Filter->SplitRequestedRegion(threadId, threadCount, splitRegion);

  ThreadDiffusionImageRegionType splitRegionDiffusionImage;
  const int total = str->Filter->SplitRequestedRegion(
      threadId, threadCount, splitRegionDiffusionImage);

  if (threadId < total)
  {
    str->Filter->ThreadedApplyUpdate(str->TimeStep, splitRegion,
                                     splitRegionDiffusionImage, threadId);
  }

  return ITK_THREAD_RETURN_VALUE;
}

template <class TInputImage, class TOutputImage>
typename AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::TimeStepType
AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::CalculateChange()
{

  itkDebugMacro(<< "CalculateChange called");

  DenseFDThreadStruct str;
  str.Filter = this;
  str.TimeStep = itk::NumericTraits<TimeStepType>::Zero;
  this->GetMultiThreader()->SetNumberOfThreads(this->GetNumberOfThreads());
  this->GetMultiThreader()->SetSingleMethod(
      this->CalculateChangeThreaderCallback, &str);

  // Initialize the list of time step values that will be generated by the
  // various threads. There is one distinct slot for each possible thread,
  // so this data structure is thread-safe.
  const int threadCount = this->GetMultiThreader()->GetNumberOfThreads();
  str.TimeStepList.resize(threadCount);
  str.ValidTimeStepList.assign(threadCount, false);

  this->GetMultiThreader()->SingleMethodExecute();

  return this->ResolveTimeStep(str.TimeStepList, str.ValidTimeStepList);
}

template <class TInputImage, class TOutputImage>
ITK_THREAD_RETURN_TYPE AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::CalculateChangeThreaderCallback(void* arg)
{

  const auto threadInfo =
      static_cast<itk::MultiThreader::ThreadInfoStruct*>(arg);
  const int threadId = threadInfo->ThreadID;
  const int threadCount = threadInfo->NumberOfThreads;
  const auto str = static_cast<DenseFDThreadStruct*>(threadInfo->UserData);

  // Execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  // Using the SplitRequestedRegion method from itk::ImageSource.
  ThreadRegionType splitRegion;
  str->Filter->SplitRequestedRegion(threadId, threadCount, splitRegion);

  ThreadDiffusionImageRegionType splitDiffusionimageRegion;
  const int total = str->Filter->SplitRequestedRegion(
      threadId, threadCount, splitDiffusionimageRegion);

  if (threadId < total)
  {
    str->TimeStepList[threadId] = str->Filter->ThreadedCalculateChange(
        splitRegion, splitDiffusionimageRegion, threadId);
    str->ValidTimeStepList[threadId] = true;
  }

  return ITK_THREAD_RETURN_VALUE;
}

template <class TInputImage, class TOutputImage>
void AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage,
    TOutputImage>::ThreadedApplyUpdate(TimeStepType dt,
                                       const ThreadRegionType& regionToProcess,
                                       const ThreadDiffusionImageRegionType&
                                           diffusionRegionToProcess,
                                       int)
{

  itk::ImageRegionIterator<UpdateBufferType> itUpdate(m_UpdateBuffer,
                                                      regionToProcess);
  itk::ImageRegionIterator<OutputImageType> itOutput(this->GetOutput(),
                                                     regionToProcess);

  itUpdate = itUpdate.Begin();
  itOutput = itOutput.Begin();

  while (!itUpdate.IsAtEnd())
  {
    itOutput.Value() += itUpdate.Value() * dt;

    ++itOutput;
    ++itUpdate;
  }
}

template <class TInputImage, class TOutputImage>
typename AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::TimeStepType
AnisotropicDiffusionVesselEnhancementImageFilter<TInputImage, TOutputImage>::
    ThreadedCalculateChange(
        const ThreadRegionType& regionToProcess,
        const ThreadDiffusionImageRegionType& diffusionRegionToProcess, int)
{
  typedef typename OutputImageType::RegionType RegionType;
  typedef typename OutputImageType::SizeType SizeType;
  typedef typename OutputImageType::SizeValueType SizeValueType;
  typedef typename OutputImageType::IndexType IndexType;
  typedef typename OutputImageType::IndexValueType IndexValueType;
  typedef typename FiniteDifferenceFunctionType::NeighborhoodType
      NeighborhoodIteratorType;
  typedef itk::ImageRegionIterator<UpdateBufferType> UpdateIteratorType;
  typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<
      DiffusionTensorImageType> DiffusionTensorRegionCalculatorType;

  typedef typename DiffusionTensorRegionCalculatorType::FaceListType
      DiffusionTensorRegionListType;
  typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<
      OutputImageType> RegionCalculatorType;
  typedef typename RegionCalculatorType::FaceListType RegionListType;

  typename OutputImageType::Pointer output = this->GetOutput();

  const typename FiniteDifferenceFunctionType::Pointer df = dynamic_cast<
      AnisotropicDiffusionVesselEnhancementFunction<UpdateBufferType>*>(
      this->GetDifferenceFunction().GetPointer());

  const SizeType radius = df->GetRadius();

  RegionCalculatorType regionCalculator;
  RegionListType regionList = regionCalculator(output, regionToProcess, radius);
  typename RegionListType::iterator itRegionList = regionList.begin();

  DiffusionTensorRegionCalculatorType diffusionTensorRegionCalculator;
  DiffusionTensorRegionListType diffusionTensorRegionList =
      diffusionTensorRegionCalculator(m_DiffusionTensorImage,
                                      diffusionRegionToProcess, radius);
  typename DiffusionTensorRegionListType::iterator itDiffTensor =
      diffusionTensorRegionList.begin();

  const DiffusionTensorNeighborhoodType dTN(radius, m_DiffusionTensorImage,
                                            *itDiffTensor);
  DerivativeStructType* derivativeData =
      static_cast<DerivativeStructType*>(df->GetGlobalDataPointer());
  UpdateIteratorType itNeighorNBUpdate(m_UpdateBuffer, *itRegionList);

  // Process the non-boundary region.
  NeighborhoodIteratorType itNeighborNonBoundary(radius, output, *itRegionList);
  itNeighborNonBoundary.GoToBegin();

  while (!itNeighborNonBoundary.IsAtEnd())
  {
    itNeighorNBUpdate.Value() =
        df->ComputeUpdate(itNeighborNonBoundary, dTN, derivativeData);
    ++itNeighborNonBoundary;
    ++itNeighorNBUpdate;
  }

  // Skip the non-boundary region
  ++itRegionList;

  // Process each of the boundary faces.
  for (; itRegionList != regionList.end(); ++itRegionList)
  {
    NeighborhoodIteratorType itNeighborBoundary =
        NeighborhoodIteratorType(radius, output, *itRegionList);
    const DiffusionTensorNeighborhoodType bDD = DiffusionTensorNeighborhoodType(
        radius, m_DiffusionTensorImage, *itDiffTensor);
    UpdateIteratorType itNeighborBUpdate =
        UpdateIteratorType(m_UpdateBuffer, *itRegionList);

    itNeighborBoundary.GoToBegin();
    itNeighborBUpdate.GoToBegin();

    while (!itNeighborBoundary.IsAtEnd())
    {
      itNeighborBUpdate.Value() =
          df->ComputeUpdate(itNeighborBoundary, bDD, derivativeData);
      ++itNeighborBoundary;
      ++itNeighborBUpdate;
    }
    ++itDiffTensor;
  }

  // Ask the finite difference function to compute the time step for
  // this iteration.
  const TimeStepType timeStep = df->ComputeGlobalTimeStep(derivativeData);
  df->ReleaseGlobalDataPointer(derivativeData);

  return timeStep;
}

template <class TInputImage, class TOutputImage>
void AnisotropicDiffusionVesselEnhancementImageFilter<
    TInputImage, TOutputImage>::GenerateData()
{
  itkDebugMacro(<< "GenerateData is called");

  if (!Superclass::GetIsInitialized())
  {
    this->AllocateOutputs();
    this->CopyInputToOutput();
    this->AllocateUpdateBuffer();
    this->AllocateDiffusionTensorImage();
    this->SetStateToInitialized();
    this->SetElapsedIterations(0);

    Superclass::SetNumberOfIterations(m_NumberOfIterations);
  }

   // std::cout << "number of iter : " << this->GetNumberOfIterations() 
   // << std::endl;

  unsigned int iter = 0;
  while (!this->Halt())
  {
    if (iter == this->GetNumberOfIterations() - 1)
    {
      std::cout << "(In AnisotropicFilter) Iteration : Final Frangi "
                   "iteration.\n";
    }
    else    
    {
      std::cout << "(In AnisotropicFilter) Iteration : " << iter + 1 
                << std::endl;
    }

    this->InitializeIteration();

    if ((m_GenerateIterationFiles && iter == 0) || this->GetFrangiOnly())
    {
        const std::string filename = "frangi_only_vesselness_measure.nii.gz";

        //WRITE OUTPUT TO FILE
        typedef itk::Image<float, ImageDimension> floatImageType;
        typedef itk::ImageFileWriter<floatImageType> ImageWriterType;

        typedef itk::CastImageFilter< OutputImageType, floatImageType > CastFilterType;
        typename CastFilterType::Pointer castFilter = CastFilterType::New();
        castFilter->SetInput(m_MultiScaleVesselnessFilter->GetOutput());

        typename ImageWriterType::Pointer writer = ImageWriterType::New();
        writer->SetFileName(filename);
        writer->SetInput(castFilter->GetOutput());
        writer->Update();


        if (this->GetFrangiOnly())
        {
            std::cout << "Only Frangi vesselness will be produced. VED will be "
                        "aborted. \n";
            break;
        }
    }

    if (iter < this->GetNumberOfIterations() - 1)
    {
        std::cout << "Apply update in diffusion.\n";
        this->ApplyUpdate(this->CalculateChange());

        // Ensure to save iteration 1 to N - 1. Because N = output.
        if (m_GenerateIterationFiles && iter < this->GetNumberOfIterations())
        {
            std::stringstream sstm;
            sstm << "ved_iteration_" << (iter + 1) << ".nii.gz";
            const std::string vedIterationFile = sstm.str();

            //WRITE OUTPUT TO FILE
            typedef itk::Image<float, ImageDimension> floatImageType;
            typedef itk::ImageFileWriter<floatImageType> ImageWriterType;

            typedef itk::CastImageFilter< OutputImageType, floatImageType > CastFilterType;
            typename CastFilterType::Pointer castFilter = CastFilterType::New();
            castFilter->SetInput(m_MultiScaleVesselnessFilter->GetOutput());

            typename ImageWriterType::Pointer writer = ImageWriterType::New();
            writer->SetFileName(vedIterationFile);
            writer->SetInput(castFilter->GetOutput());
            writer->Update();
        }
    }   

    this->SetElapsedIterations(++iter);

    this->InvokeEvent(itk::IterationEvent());
    if (this->GetAbortGenerateData())
    {
      this->InvokeEvent(itk::IterationEvent());
      this->ResetPipeline();
      throw itk::ProcessAborted(__FILE__, __LINE__);
    }
  }
}

#endif
