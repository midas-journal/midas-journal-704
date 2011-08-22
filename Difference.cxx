#include <stdlib.h>
#include <stdio.h>

#include "itkImageFileReader.h"
#include "itkConstrainedValueDifferenceImageFilter.h"
#include "itkImageFileWriter.h"

int main( int argc, char *argv[] )
{
  char   * in1_file_name = argv[1];
  char   * in2_file_name = argv[2];
  char   * out_file_name = argv[3];

  typedef float FloatPixelType;

  const unsigned int Dimension = 3;

  typedef itk::Image< FloatPixelType, Dimension >  FloatImageType;

  // Reader
  typedef  itk::ImageFileReader< FloatImageType > ReaderType;

  ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( in1_file_name );

  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( in2_file_name );

  // Difference
  typedef itk::ConstrainedValueDifferenceImageFilter
    < FloatImageType, FloatImageType, FloatImageType > DiffType;
  DiffType::Pointer diff = DiffType::New();
  diff->SetInput1( reader1->GetOutput() );
  diff->SetInput2( reader2->GetOutput() );

  // Writer
  typedef itk::ImageFileWriter< FloatImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( out_file_name );
  writer->SetInput( diff->GetOutput() );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
