#ifndef __itkSymmetricEigenVectorAnalysisImageFilter_h
#define __itkSymmetricEigenVectorAnalysisImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "itkSymmetricEigenAnalysis.h"

// This functor class invokes the computation of Eigen Analysis for
// every pixel. The input pixel type must provide the API for the [][]
// operator, while the output pixel type must provide the API for the
// [] operator. Input pixel matrices should be symmetric.
//
// The default operation is to order eigen values in ascending order.
// You may also use OrderEigenValuesBy( ) to order eigen values by
// magnitude as is common with use of tensors in vessel extraction.

template <typename TInput, typename TOutput, typename TMatrix>
class SymmetricEigenVectorAnalysisFunction
{
public:
  SymmetricEigenVectorAnalysisFunction() {}
  ~SymmetricEigenVectorAnalysisFunction() {}
  typedef itk::SymmetricEigenAnalysis<TInput, TOutput, TMatrix> CalculatorType;

  inline TMatrix operator()(const TInput& x)
  {
    TOutput eigenValues;
    TMatrix eigenVectorMatrix;
    m_Calculator.ComputeEigenValuesAndVectors(x, eigenValues,
                                              eigenVectorMatrix);
    return eigenVectorMatrix;
  }

  void SetDimension(unsigned int n) { m_Calculator.SetDimension(n); }

  // Typdedefs to order eigen values.
  // OrderByValue:      lambda_1 < lambda_2 < ....
  // OrderByMagnitude:  |lambda_1| < |lambda_2| < .....
  // DoNotOrder:        Default order of eigen values obtained after QL method

  typedef enum
  {
    OrderByValue = 1,
    OrderByMagnitude,
    DoNotOrder
  } EigenValueOrderType;

  // Order eigen values. Default is to OrderByValue:  lambda_1 < lambda_2 < ...
  void OrderEigenValuesBy(EigenValueOrderType order)
  {
    if (order == OrderByMagnitude)
    {
      m_Calculator.SetOrderEigenMagnitudes(true);
    }
    else if (order == DoNotOrder)
    {
      m_Calculator.SetOrderEigenValues(false);
    }
  }

private:
  CalculatorType m_Calculator;
};

// \class SymmetricEigenVectorAnalysisImageFilter //
//\ingroup IntensityImageFilters  Multithreaded  TensorObjects
template <typename TInputImage, typename TOutputImage, typename TOutputMatrix>
class SymmetricEigenVectorAnalysisImageFilter
    : public itk::UnaryFunctorImageFilter<
          TInputImage, TOutputMatrix,
          SymmetricEigenVectorAnalysisFunction<
              typename TInputImage::PixelType, typename TOutputImage::PixelType,
              typename TOutputMatrix::PixelType>>
{
public:
  typedef SymmetricEigenVectorAnalysisImageFilter Self;
  typedef itk::UnaryFunctorImageFilter<
      TInputImage, TOutputMatrix,
      SymmetricEigenVectorAnalysisFunction<
          typename TInputImage::PixelType, typename TOutputImage::PixelType,
          typename TOutputMatrix::PixelType>> Superclass;

  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef typename Superclass::OutputImageType OutputImageType;
  typedef typename TOutputImage::PixelType OutputPixelType;
  typedef typename TInputImage::PixelType InputPixelType;
  typedef typename Superclass::FunctorType FunctorType;

  // Typdedefs to order eigen values.
  // OrderByValue:      lambda_1 < lambda_2 < ....
  // OrderByMagnitude:  |lambda_1| < |lambda_2| < .....
  // DoNotOrder:        Default order of eigen values obtained after QL method

  typedef typename FunctorType::EigenValueOrderType EigenValueOrderType;

  // Order eigen values. Default is to OrderByValue:  lambda_1 < lambda_2 < ....
  void OrderEigenValuesBy(EigenValueOrderType order)
  {
    this->GetFunctor().OrderEigenValuesBy(order);
  }

  itkNewMacro(Self);

  // Set the dimension of the tensor. (For example the
  // SymmetricSecondRankTensor is a pxp matrix)
  void SetDimension(unsigned int p) { this->GetFunctor().SetDimension(p); }

protected:
  SymmetricEigenVectorAnalysisImageFilter(){};
  virtual ~SymmetricEigenVectorAnalysisImageFilter(){};

private:
  SymmetricEigenVectorAnalysisImageFilter(
      const Self&);            // purposely not implemented
  void operator=(const Self&); // purposely not implemented
};

#endif
