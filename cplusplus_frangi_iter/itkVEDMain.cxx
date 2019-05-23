#if defined(_MSC_VER)
#pragma warning(disable : 4786)
#endif

#include "itkAnisotropicDiffusionVesselEnhancementImageFilter.h"
#include "itkSymmetricEigenVectorAnalysisImageFilter.h"

#include "boost/program_options.hpp"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkLiThresholdImageFilter.h"
#include "itkHuangThresholdImageFilter.h"
#include "itkIntermodesThresholdImageFilter.h"
#include "itkIsoDataThresholdImageFilter.h"
#include "itkKittlerIllingworthThresholdImageFilter.h"
#include "itkMaximumEntropyThresholdImageFilter.h"
#include "itkMomentsThresholdImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkRenyiEntropyThresholdImageFilter.h"
#include "itkShanbhagThresholdImageFilter.h"
#include "itkTriangleThresholdImageFilter.h"
#include "itkYenThresholdImageFilter.h"
#include "itkStatisticsImageFilter.h"

#include "itkCurvatureFlowImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkAnisotropicDiffusionImageFilter.h"
#include "itkDenseFiniteDifferenceImageFilter.h"

bool process_command_line(int argc, char** argv,
                          boost::program_options::variables_map& vm)
{
  try
  {
    boost::program_options::options_description program(
        "Program allowed options\n");
    program.add_options()("help,h", "produce help message.");

    boost::program_options::options_description requiredVariable("Required\n");
    requiredVariable.add_options()(
        "input,i", boost::program_options::value<std::string>()->required(),
        "the input file name.")(
        "output,o", boost::program_options::value<std::string>()->required(),
        "the output file name.");

    boost::program_options::options_description multipleHessianVariable(
        "Multi-scale Hessian\n");
    multipleHessianVariable.add_options()(
        "sigmaMin,m",
        boost::program_options::value<double>()->default_value(0.5),
        "The minimum sigma used for the multi-scale analysis. Usually the "
        "smallest resolution of the acquisition.")(
        "sigmaMax,M",
        boost::program_options::value<double>()->default_value(4.0),
        "The maximum sigma used for the multi-scale analysis. Usually the "
        "biggest expected size of a vessel.")(
        "numberOfScale,n",
        boost::program_options::value<int>()->default_value(5),
        "The number of scales created between the sigma min/max in the "
        "multi-scale analysis.");

    boost::program_options::options_description vesselnessVariable(
        "Frangi vesselness measure\n");
    vesselnessVariable.add_options()("darkBlood,d",
                                     "Flag to extract black blood vessel.")(
        "alpha,a", boost::program_options::value<double>()->default_value(0.5),
        "The alpha parameter used in Frangi vesselness equation to limit "
        "blob-like structure.")(
        "beta,b", boost::program_options::value<double>()->default_value(1.0),
        "The beta parameter used in Frangi vesselness "
        "equation to limit plate-like structure.")(
        "c,c", boost::program_options::value<double>()->default_value(0.00001),
        "The c parameter used in Frangi vesselness equation to limit the "
        "background comparison.");

    boost::program_options::options_description vedVariable(
        "Vessel Enhancing Diffusion\n");
    vedVariable.add_options()(
        "numberOfIteration,t",
        boost::program_options::value<int>()->default_value(1),
        "The number of diffusion iterations.")(
        "sensitivity,s",
        boost::program_options::value<double>()->default_value(5.0),
        "The sensitivity used in VED param.")(
        "wStrength,w",
        boost::program_options::value<double>()->default_value(15.0),
        "The weigthed strength used in VED param.")(
        "epsilon,e",
        boost::program_options::value<double>()->default_value(1.0),
        "The epsilon used in VED param.");

    boost::program_options::options_description flagVariable("Flags\n");
    flagVariable.add_options()("frangiOnly,f", "Flag to stop the pipeline "
                                               "after Frangi equation (No "
                                               "post-processing with VED).")(
        "scaleObject,O",
        "Flag to rescale vesselness measure based on eigen amplitude.")(
        "generateScale,S",
        "Flag to generate an output for the best sigma value per voxel.")(
        "generateHessian,H",
        "Flag to generate an output with the hessian matrix per voxel.")(
        "generateIterationFiles,I",
        "Flag to generate output iteration files and vesselness.");

    boost::program_options::options_description global;

    global.add(program)
        .add(requiredVariable)
        .add(multipleHessianVariable)
        .add(vesselnessVariable)
        .add(vedVariable)
        .add(flagVariable);

    boost::program_options::store(
        boost::program_options::parse_command_line(argc, argv, global), vm);

    if (vm.count("help"))
    {
      std::cout << global << std::endl;
      return false;
    }

    boost::program_options::notify(vm);
  }
  catch (std::exception& e)
  {
    std::cerr << "Error: " << e.what() << "\n";
    return false;
  }
  catch (...)
  {
    std::cerr << "Unknown error!\n";
    return false;
  }
  return true;
}

int main(int argc, char* argv[])
{
  const int Dimension = 3;
  typedef double InputPixelType;
  typedef double OutputPixelType;

  typedef itk::Image<InputPixelType, Dimension> InputImageType;
  typedef itk::Image<InputPixelType, Dimension> OutputImageType;
  typedef itk::ImageFileReader<InputImageType> ImageReaderType;

  typedef itk::SymmetricSecondRankTensor<double, Dimension> TensorPixelType;
  typedef itk::Image<TensorPixelType, Dimension> TensorImageType;

  typedef float ScalesPixelType;
  typedef itk::Image<ScalesPixelType, Dimension> ScalesImageType;

  boost::program_options::variables_map vm;
  if (!process_command_line(argc, argv, vm))
  {
    return 1;
  }

  ImageReaderType::Pointer reader = ImageReaderType::New();
  std::cout << "Reading input image : " << vm["input"].as<std::string>()
            << std::endl;

  reader->SetFileName(vm["input"].as<std::string>());

  try
  {
    reader->Update();
  }
  catch (itk::ExceptionObject& err)
  {
    std::cerr << "Exception thrown: " << err << std::endl;
    return EXIT_FAILURE;
  }

  typedef AnisotropicDiffusionVesselEnhancementImageFilter<
      InputImageType, OutputImageType> VesselnessFilterType;

  // Create a vesselness Filter.
  VesselnessFilterType::Pointer VesselnessFilter = VesselnessFilterType::New();
  VesselnessFilter->SetInput(reader->GetOutput());

  // Multi-scale Hessian parameters
  VesselnessFilter->SetSigmaMin(vm["sigmaMin"].as<double>());
  VesselnessFilter->SetSigmaMax(vm["sigmaMax"].as<double>());
  VesselnessFilter->SetNumberOfSigmaSteps(vm["numberOfScale"].as<int>());

  // Frangi vesselness equation parameters
  if (vm.count("darkBlood"))
  {
    VesselnessFilter->SetBrightBlood(false);
    std::cout << "Will extract dark blood.\n";
  }
  else
  {
    std::cout << "Will extract bright blood.\n";
  }

  VesselnessFilter->SetAlpha(vm["alpha"].as<double>());
  VesselnessFilter->SetBeta(vm["beta"].as<double>());
  VesselnessFilter->SetC(vm["c"].as<double>());

  // Vessel enhancing diffusion parameters

  // Add + 1 to iteration, to apply frangi to last VED iteration.
  VesselnessFilter->SetNumberOfIterations(vm["numberOfIteration"].as<int>() + 1);
  VesselnessFilter->SetSensitivity(vm["sensitivity"].as<double>());
  VesselnessFilter->SetWStrength(vm["wStrength"].as<double>());
  VesselnessFilter->SetEpsilon(vm["epsilon"].as<double>());

  // Flags
  if (vm.count("frangiOnly"))
  {
    VesselnessFilter->SetFrangiOnly(true);
    std::cout << "Will generate the Frangi vesselness measure image only.\n";
  }

  if (vm.count("scaleObject"))
  {
    VesselnessFilter->SetScaleObject(true);
    std::cout
        << "Will scale the vesselness based on the eigen value amplitude.\n";
  }

  if (vm.count("generateScale"))
  {
    VesselnessFilter->SetGenerateScale(true);
    std::cout << "Will generate the best scale image.\n";
  }

  if (vm.count("generateHessian"))
  {
    VesselnessFilter->SetGenerateHessian(true);
    std::cout << "Will generate the best hessian image.\n";
  }

  if (vm.count("generateIterationFiles"))
  {
    VesselnessFilter->SetGenerateIterationFiles(true);
    std::cout << "Will generate the iteration files\n";
  }

  try
  {
    VesselnessFilter->Update();
  }
  catch (itk::ExceptionObject& err)
  {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Writing out the enhanced image to "
            << vm["output"].as<std::string>() << std::endl;

  typedef itk::Image<float, Dimension> floatImageType;
  typedef itk::ImageFileWriter<floatImageType> ImageWriterType;

  typedef itk::CastImageFilter< OutputImageType, floatImageType > CastFilterType;
  typename CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(VesselnessFilter->GetOutput());

  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName(vm["output"].as<std::string>());
  writer->SetInput(castFilter->GetOutput());

  try
  {
    writer->Update();
  }
  catch (itk::ExceptionObject& err)
  {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
  }

  //DO SEGMENTATIONS!
  bool smooth = true;
  if (smooth == true) 
  {
    typedef float OutputPixelType;
    typedef itk::Image< OutputPixelType, 3 > OutputSmoothImageType;
 
    typedef itk::CurvatureFlowImageFilter< OutputSmoothImageType, OutputSmoothImageType > 
      CurvatureFlowImageFilterType;
    typedef itk::GradientAnisotropicDiffusionImageFilter<OutputSmoothImageType, OutputSmoothImageType > 
      GradientAnisotropicDiffusionFilterType ;
    typedef itk::AnisotropicDiffusionImageFilter<OutputSmoothImageType, OutputSmoothImageType> 
      AnisotropicDiffusionFilterType ;
  
    typedef std::map<std::string, itk::AnisotropicDiffusionImageFilter<OutputSmoothImageType, OutputSmoothImageType>::Pointer> FilterContainerType;
      FilterContainerType filterContainer;
    //filterContainer["CurvatureFlow"] = CurvatureFlowImageFilterType::New();
    filterContainer["GradientAnisotropicDiffusion"] = GradientAnisotropicDiffusionFilterType::New();
    //filterContainer["AnisotropicDiffusion"] = AnisotropicDiffusionFilterType::New();
    
    FilterContainerType::iterator it = filterContainer.begin();
    int iterations = 5;
    for (it = filterContainer.begin(); it != filterContainer.end(); ++it) 
    {
      for (int iter = 10; iter <= 11.0 ; iter += 3)
      {
      std::cout << "Writing out segmentations (" << (*it).first << "). \n" ;
      (*it).second->SetConductanceParameter( 0.3 ); //0.5
      (*it).second->SetTimeStep( 0.02 );
      (*it).second->SetInput( castFilter->GetOutput() );
      (*it).second->SetNumberOfIterations( iter );
      (*it).second->Update();

      typedef itk::ImageFileWriter<OutputSmoothImageType> SegWriterType;
      typename SegWriterType::Pointer writer = SegWriterType::New();
      //std::string message = std::string("Ved_") + (*it).first + std::to_string(fConductance) + std::string(".nii.gz") ;  
      std::string message = std::string("Ved_") + (*it).first + std::to_string(iter) + std::string(".nii.gz") ;  
      writer->SetFileName(message);
      writer->SetInput((*it).second->GetOutput());
      writer->Update();
      }
    }

  }
 
  //DO SEGMENTATIONS!
  bool provideSeg = false;
  if (provideSeg == true) 
  {
    typedef  unsigned char  OutputPixelType;
    typedef itk::Image< OutputPixelType, 3 >   OutputSegType;
 
    typedef itk::LiThresholdImageFilter<floatImageType, OutputSegType >
      LiFilterType;
    typedef itk::HuangThresholdImageFilter<floatImageType, OutputSegType >
      HuangFilterType;
    typedef itk::IntermodesThresholdImageFilter<floatImageType, OutputSegType >
      IntermodesFilterType;
    typedef itk::IsoDataThresholdImageFilter<floatImageType, OutputSegType >
      IsoDataFilterType;
    typedef itk::KittlerIllingworthThresholdImageFilter<floatImageType, OutputSegType >
      KittlerIllingworthFilterType;
    typedef itk::LiThresholdImageFilter<floatImageType, OutputSegType >
      LiFilterType;
    typedef itk::MaximumEntropyThresholdImageFilter<floatImageType, OutputSegType >
      MaximumEntropyFilterType;
    typedef itk::MomentsThresholdImageFilter<floatImageType, OutputSegType >
      MomentsFilterType;
    typedef itk::OtsuThresholdImageFilter<floatImageType, OutputSegType >
      OtsuFilterType;
    typedef itk::RenyiEntropyThresholdImageFilter<floatImageType, OutputSegType >
      RenyiEntropyFilterType;
    typedef itk::ShanbhagThresholdImageFilter<floatImageType, OutputSegType >
      ShanbhagFilterType;
    typedef itk::TriangleThresholdImageFilter<floatImageType, OutputSegType >
      TriangleFilterType;
    typedef itk::YenThresholdImageFilter<floatImageType, OutputSegType >
      YenFilterType;
  
    typedef itk::StatisticsImageFilter<floatImageType> StatisticsImageFilterType;
    StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New();
    statisticsImageFilter->SetInput(castFilter->GetOutput());
    statisticsImageFilter->Update();

    typedef std::map<std::string, itk::HistogramThresholdImageFilter<floatImageType, OutputSegType>::Pointer> FilterContainerType;
    FilterContainerType filterContainer;
  
    filterContainer["Huang"] = HuangFilterType::New();
    //filterContainer["Intermodes"] = IntermodesFilterType::New();
    //filterContainer["IsoData"] = IsoDataFilterType::New();
    //filterContainer["KittlerIllingworth"] = KittlerIllingworthFilterType::New();
    //filterContainer["Li"] = LiFilterType::New();
    //filterContainer["MaximumEntropy"] = MaximumEntropyFilterType::New();
    filterContainer["Moments"] = MomentsFilterType::New();
    filterContainer["Otsu"] = OtsuFilterType::New();
    //filterContainer["RenyiEntropy"] = RenyiEntropyFilterType::New();
    //filterContainer["Shanbhag"] = ShanbhagFilterType::New();
    filterContainer["Triangle"] = TriangleFilterType::New();
    //filterContainer["Yen"] = YenFilterType::New();
  
    FilterContainerType::iterator it = filterContainer.begin();
    for (it = filterContainer.begin(); it != filterContainer.end(); ++it)
      {
      std::cout << "Writing out segmentations (" << (*it).first << "). \n" ;
      (*it).second->SetInsideValue( statisticsImageFilter->GetMaximum() );
      (*it).second->SetOutsideValue( statisticsImageFilter->GetMinimum() );
      (*it).second->SetInput( castFilter->GetOutput() );
      (*it).second->SetNumberOfHistogramBins( 100 );
      (*it).second->Update();

      typedef itk::ImageFileWriter<OutputSegType> SegWriterType;
      typename SegWriterType::Pointer writer = SegWriterType::New();
      std::string message = std::string("Ved_") + (*it).first + std::string(".nii.gz") ;  
      writer->SetFileName(message);
      writer->SetInput((*it).second->GetOutput());
      writer->Update();
      }

  }
 

  if (vm.count("generateScale"))
  {

    std::cout << "Writing out the best sigma scale image. \n";
    typedef itk::ImageFileWriter<ScalesImageType> ImageWriterType;
    ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetFileName("ved_generated_best_scale.nii.gz");
    writer->SetInput(VesselnessFilter->GetScalesOutput());

    try
    {
      writer->Update();
    }
    catch (itk::ExceptionObject& err)
    {
      std::cerr << "Exception caught: " << err << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (vm.count("generateHessian"))
  {

    std::cout << "Writing out the best hessian matrix image. \n";

    typedef itk::ImageFileWriter<TensorImageType> ImageWriterType;
    ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetFileName("ved_generated_best_Hessian.nii.gz");
    writer->SetInput(VesselnessFilter->GetHessianOutput());

    try
    {
      writer->Update();
    }
    catch (itk::ExceptionObject& err)
    {
      std::cerr << "Exception caught: " << err << std::endl;
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}
