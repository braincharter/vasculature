#ifndef __itkVesselnessMeasurement_h
#define __itkVesselnessMeasurement_h

#include "itkSymmetricSecondRankTensor.h"
#include "itkSymmetricEigenVectorAnalysisImageFilter.h"
#include "itkImageToImageFilter.h"

//*\class HessianToObjectnessMeasureImageFilter
//	\brief A filter to enhance M-dimensional objects in N-dimensional images
//
// The objectness measure is a generalization of Frangi's vesselness measure,
// which is based on the analysis of the the Hessian eigen system. The
// filter can enhance blob-like structures (M=0), vessel-like structures (M=1),
// 2D plate-like structures (M=2), hyper-plate-like structures (M=3) in
// N-dimensional images, with M<N. The filter takes an image of a Hessian
// pixels ( SymmetricSecondRankTensorpixels pixels ) and produces an enhanced
// image. The Hessian input image can be produced using
// itk::HessianRecursiveGaussianImageFilter.
//
//	Frangi, AF, Niessen, WJ, Vincken, KL, & Viergever, MA (1998). Multiscale
//	Vessel Enhancement Filtering. In Wells, WM, Colchester, A, & Delp, S,
//  Editors, MICCAI '98 Medical Image Computing and Computer-Assisted
//  Intervention, Lecture Notes in Computer Science, pages 130-137, Springer
//  Verlag, 1998.
//
//  Manniesing, R, Viergever, MA, & Niessen, WJ (2006). Vessel Enhancing
//  Diffusion: A Scale Space Representation of Vessel Structures. Medical
//  Image Analysis, 10(6), 815-825./
//
//  Coded by :
//  http://www.insight-journal.org/browse/publication/163
//

template <typename TInputImage, typename TOutputImage>
class VesselnessMeasurement
    : public itk::ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  typedef VesselnessMeasurement Self;

  typedef itk::ImageToImageFilter<TInputImage, TOutputImage> Superclass;

  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef typename Superclass::InputImageType InputImageType;
  typedef typename Superclass::OutputImageType OutputImageType;
  typedef typename InputImageType::PixelType InputPixelType;
  typedef typename OutputImageType::PixelType OutputPixelType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  static const unsigned int ImageDimension = InputImageType::ImageDimension;

  typedef double EigenValueType;
  typedef itk::FixedArray<EigenValueType, ImageDimension> EigenValueArrayType;
  
  typedef itk::Image<double, OutputImageType::ImageDimension> LambdaImageType;
  
    //ADDED/////////
  static const unsigned int vectorlength = 6; 
  typedef itk::Vector<EigenValueType, vectorlength> 
    EigenVectorArrayType;
  typedef itk::Vector<EigenValueType, vectorlength> 
    VectorType;
  typedef itk::Matrix<EigenValueType, OutputImageType::ImageDimension, OutputImageType::ImageDimension>   
    EigenVectorMatrixType;
  typedef itk::Matrix<EigenValueType, OutputImageType::ImageDimension, OutputImageType::ImageDimension>   
    MatrixType;
  typedef itk::Image<EigenVectorArrayType, OutputImageType::ImageDimension>                
    FlippedHessianArrayImageType;
  ///////////////


  void UpdateMeasure(OutputImageType* newMeasure);

  itkNewMacro(Self);

  itkTypeMacro(VesselnessMeasurement, ImageToImageFilter);

  // Set/Get Alpha, the weight corresponding to R_A (the ratio of the smallest
  // eigenvalue that has to be large to the larger  ones). Smaller values lead
  // to increased sensitivity to the object dimensionality.
  itkSetMacro(Alpha, double);
  itkGetConstMacro(Alpha, double);

  // Set/Get Beta, the weight corresponding to R_B (the ratio of the largest
  // eigenvalue that has to be small to the largerones). Smaller values lead to
  // increased sensitivity to the object dimensionality.
  itkSetMacro(Beta, double);
  itkGetConstMacro(Beta, double);

  // Set/Get Gamma, the weight corresponding to S (the Frobenius norm of the
  // Hessian matrix, or second-order structureness)
  itkSetMacro(Gamma, double);
  itkGetConstMacro(Gamma, double);

  // Set/Get C, the weight corresponding to the last objectness measure
  itkSetMacro(C, double);
  itkGetConstMacro(C, double);

  itkSetMacro(MinLambda, double);
  itkGetConstMacro(MinLambda, double);
  
  // Toggle scaling the objectness measure with the magnitude of the
  // largestabsolute eigenvalue
  itkSetMacro(ScaleObjectnessMeasure, bool);
  itkGetConstMacro(ScaleObjectnessMeasure, bool);
  itkBooleanMacro(ScaleObjectnessMeasure);

  // Enhance bright structures on a dark background if true, the opposite if
  // false.
  itkSetMacro(BrightObject, bool);
  itkGetConstMacro(BrightObject, bool);
  itkBooleanMacro(BrightObject);

  itkSetMacro(FrangiOnly, bool);
  itkGetConstMacro(FrangiOnly, bool);
  itkBooleanMacro(FrangiOnly);
  
  itkSetMacro(FirstPass, bool);
  itkGetConstMacro(FirstPass, bool);
  itkBooleanMacro(FirstPass);
  
  itkSetMacro(LastHighLambda3, typename LambdaImageType::Pointer);
  itkGetConstMacro(LastHighLambda3, typename LambdaImageType::Pointer);
  
  itkSetMacro(FlippedHessian, typename FlippedHessianArrayImageType::Pointer);
  itkGetConstMacro(FlippedHessian, typename FlippedHessianArrayImageType::Pointer);

#ifdef ITK_USE_CONCEPT_CHECKING
  itkConceptMacro(DoubleConvertibleToOutputCheck,
                  (itk::Concept::Convertible<double, OutputPixelType>));
#endif

protected:
  VesselnessMeasurement();
  VesselnessMeasurement(double const& alpha, double const& beta,
                        double const& c, const bool scale, const bool bright,
                        const bool frangi, const bool firstPass, double minL3);

  ~VesselnessMeasurement() {}
  void PrintSelf(std::ostream& os, itk::Indent indent) const;

  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                            itk::ThreadIdType threadId);
  
  

private:
  VesselnessMeasurement(const Self&);

  void operator=(const Self&);

  // functor used to sort the eigenvalues are to be sorted
  // |e1|<=|e2|<=...<=|eN|
  // \class AbsLessEqualCompare
  // \brief Returns ( abs(a) <= abs(b) )
  // \ingroup ITKReview
  struct AbsLessEqualCompare
  {
    bool operator()(EigenValueType a, EigenValueType b)
    {
      return vnl_math_abs(a) <= vnl_math_abs(b);
    }
  };

  typename LambdaImageType::Pointer m_LastHighLambda3;
  typename FlippedHessianArrayImageType::Pointer m_FlippedHessian;

/*
        FORMAT
      myTensor = mIt.Get();
      myMatrix[0][0] = myTensor[0]; D11
      myMatrix[0][1] = myTensor[1]; D12
      myMatrix[0][2] = myTensor[2]; D13
      myMatrix[1][0] = myTensor[1]; D21
      myMatrix[1][1] = myTensor[3]; D22
      myMatrix[1][2] = myTensor[4]; D23
      myMatrix[2][0] = myTensor[2]; D31
      myMatrix[2][1] = myTensor[4]; D32
      myMatrix[2][2] = myTensor[5]; D33

      For MRTrix, they want: (volumes 0-5) D11, D22, D33, D12, D13, D23 
      myTensor[0], myTensor[3], myTensor[5], myTensor[1], myTensor[2], myTensor[4]
*/

  double m_MinLambda;
  double m_Alpha;
  double m_Beta;
  double m_Gamma;
  double m_C;
  bool m_FirstPass;
  bool m_BrightObject;
  bool m_ScaleObjectnessMeasure;
  bool m_FrangiOnly;
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVesselnessMeasurement.hxx"
#endif

#endif