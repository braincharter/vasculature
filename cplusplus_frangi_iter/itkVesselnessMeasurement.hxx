#ifndef __itkVesselnessMeasurement_hxx
#define __itkVesselnessMeasurement_hxx

#include "itkVesselnessMeasurement.h"

#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"
#include "itkSymmetricEigenAnalysis.h"
#include "vnl/vnl_math.h"

#include <algorithm>

#define EPSILON 1e-03

template <typename TInputImage, typename TOutputImage>
VesselnessMeasurement<TInputImage, TOutputImage>::VesselnessMeasurement()
    : VesselnessMeasurement<TInputImage, TOutputImage>::VesselnessMeasurement(
          0.5, 1.0, 10e-6, false, true, false)
{
}

template <typename TInputImage, typename TOutputImage>
VesselnessMeasurement<TInputImage, TOutputImage>::VesselnessMeasurement(
    double const& alpha, double const& beta, double const& c, const bool scale,
    const bool bright, const bool frangi)
    : m_Alpha{alpha}, m_Beta{beta}, m_C{c}, m_ScaleObjectnessMeasure{scale},
      m_BrightObject{bright}, m_FrangiOnly{frangi}
{
}

// =============================================================================
// Threading functions to generate the vesselness measure. This is called from
// multiScaleHessian.
// =============================================================================
template <typename TInputImage, typename TOutputImage>
void VesselnessMeasurement<TInputImage, TOutputImage>::ThreadedGenerateData(
    const OutputImageRegionType& outputRegionForThread,
    itk::ThreadIdType threadId)
{

  typename OutputImageType::Pointer output = this->GetOutput();
  typename InputImageType::ConstPointer input = this->GetInput();

  itk::ProgressReporter progress(this, threadId,
                                 outputRegionForThread.GetNumberOfPixels(),
                                 1000 / this->GetNumberOfThreads());

  // calculator for computation of the eigen values
  typedef itk::SymmetricEigenAnalysis<InputPixelType, EigenValueArrayType>
      CalculatorType;
  CalculatorType eigenCalculator(ImageDimension);

  // walk the region of eigen values and get the vesselness measure
  itk::ImageRegionConstIterator<InputImageType> it(input,
                                                   outputRegionForThread);
  itk::ImageRegionIterator<OutputImageType> oit(output, outputRegionForThread);

  it.GoToBegin();
  while (!it.IsAtEnd())
  {
    EigenValueArrayType eigenValues;
    eigenCalculator.ComputeEigenValues(it.Get(), eigenValues);
    const double S =
        vcl_sqrt(vnl_math_sqr(eigenValues[0]) + vnl_math_sqr(eigenValues[1]) +
                 vnl_math_sqr(eigenValues[2]));

    // use same m_Gamma value as Maxime Descoteaux, which is the max S value.
    m_Gamma = std::max(m_Gamma, S / 2.0);

    ++it;
  }

  oit.GoToBegin();
  it.GoToBegin();

  while (!it.IsAtEnd())
  {
    EigenValueArrayType eigenValues;
      
    eigenCalculator.ComputeEigenValues(it.Get(), eigenValues);

    // The eigenvalues are to be sorted |e1|<=|e2|<=...<=|eN|
    EigenValueArrayType sortedEigenValues = eigenValues;
    std::sort(sortedEigenValues.Begin(), sortedEigenValues.End(),
              AbsLessEqualCompare());

    // Following code will process the Frangi vesselness measure and the v4 from
    // VED.
    const double lambda1 = sortedEigenValues[0];
    const double lambda2 = sortedEigenValues[1];
    const double lambda3 = sortedEigenValues[2];

    if (m_BrightObject)
    {
      // Doing bright extraction then, If blood is dark, skip.
      if (lambda2 >= 0.0 || lambda3 >= 0.0 || vnl_math_abs(lambda2) < EPSILON ||
          vnl_math_abs(lambda3) < EPSILON)
      {
        oit.Set(itk::NumericTraits<OutputPixelType>::Zero);
        ++it;
        ++oit;
        progress.CompletedPixel();
        continue;
      }
    }
    else
    {
      // Doing dark extraction then, If blood is bright, skip.
      if (lambda2 <= 0.0 || lambda3 <= 0.0 || vnl_math_abs(lambda2) < EPSILON ||
          vnl_math_abs(lambda3) < EPSILON)
      {
        oit.Set(itk::NumericTraits<OutputPixelType>::Zero);
        ++it;
        ++oit;
        progress.CompletedPixel();
        continue;
      }
    }

    const double lambda1Abs = vnl_math_abs(lambda1);
    const double lambda2Abs = vnl_math_abs(lambda2);
    const double lambda3Abs = vnl_math_abs(lambda3);

    const double lambda1Sqr = vnl_math_sqr(lambda1);
    const double lambda2Sqr = vnl_math_sqr(lambda2);
    const double lambda3Sqr = vnl_math_sqr(lambda3);

    const double alphaSqr = vnl_math_sqr(m_Alpha);
    const double betaSqr = vnl_math_sqr(m_Beta);
    const double gammaSqr = vnl_math_sqr(m_Gamma);

    const double A = lambda2Abs / lambda3Abs;
    const double B = lambda1Abs / vcl_sqrt(vnl_math_abs(lambda2 * lambda3));
    const double C = lambda1Sqr + lambda2Sqr + lambda3Sqr;

    const double vesMeasure1 =
        1 - vcl_exp(-1.0 * (vnl_math_sqr(A) / (2.0 * alphaSqr)));

    const double vesMeasure2 =
        vcl_exp(-1.0 * (vnl_math_sqr(B) / (2.0 * betaSqr)));

    const double vesMeasure3 =
        1 - vcl_exp(-1.0 * (C / (2.0 * gammaSqr)));

    double vesselnessMeasure = vesMeasure1 * vesMeasure2 * vesMeasure3;

    // If the user wants the pure frangi vesselness measure, skip this.
    if (!m_FrangiOnly)
    {
      const double vesMeasure4 =
          vcl_exp(-1.0 * (2.0 * vnl_math_sqr(m_C)) / (lambda2Abs * lambda3Sqr));
      vesselnessMeasure = vesselnessMeasure * vesMeasure4;
    }

    if (m_ScaleObjectnessMeasure)
    {
      oit.Set(static_cast<OutputPixelType>(lambda3Abs * vesselnessMeasure));
    }
    else
    {
      oit.Set(static_cast<OutputPixelType>(vesselnessMeasure));
    }

    ++it;
    ++oit;
    progress.CompletedPixel();
  }
}

template <typename TInputImage, typename TOutputImage>
void VesselnessMeasurement<TInputImage, TOutputImage>::PrintSelf(
    std::ostream& os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Alpha: " << m_Alpha << std::endl;
  os << indent << "Beta: " << m_Beta << std::endl;
  os << indent << "C: " << m_C << std::endl;
  os << indent << "ScaleObjectnessMeasure: " << m_ScaleObjectnessMeasure
     << std::endl;
  os << indent << "BrightObject: " << m_BrightObject << std::endl;
  os << indent << "FrangiOnly: " << m_FrangiOnly << std::endl;
}

#endif