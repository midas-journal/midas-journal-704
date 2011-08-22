#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#include "itkImage.h"
#include "itkImageRegion.h"

#include "itkImageRegionNonCubeSplitter.h"

#include "itkImageFileReader.h"
#include "itkExtractImageFilter.h"
#include "itkFastMarchingImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkPasteImageFilter.h"
#include "itkImageFileWriter.h"

int main( int argc, char *argv[] )
{

  char   * in_file_name     =       argv[1];
  char   * seed_file_name   =       argv[2];
  char   * out_file_name    =       argv[3];
  int      max_iterations   = atoi( argv[4] );

  typedef float FloatPixelType;

  const unsigned int Dimension = 3;

  typedef itk::Image< FloatPixelType, Dimension >  FloatImageType;
  typedef FloatImageType::RegionType               RegionType;
  typedef FloatImageType::IndexType                IndexType;
  typedef FloatImageType::SizeType                 SizeType;
  typedef FloatImageType::PointType                PointType;

  int mpi_rank;
  int mpi_size;

  // Initialise MPI
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );

  // Reader
  typedef  itk::ImageFileReader< FloatImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( in_file_name );
  reader->UpdateOutputInformation(); // ITK Bug 9766 requires that the
                                       // reader be fully updated. 
  //reader->Update();
  FloatImageType::Pointer input_image = reader->GetOutput();
  PointType  input_origin = input_image->GetOrigin();
  RegionType input_region = input_image->GetLargestPossibleRegion();
  IndexType  input_index  = input_region.GetIndex();
  SizeType   input_size   = input_region.GetSize();

  // Writer Create the file to which we will write in parallel.
  typedef itk::ImageFileWriter< FloatImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( reader->GetOutput() );
  writer->SetFileName( out_file_name );
  writer->SetNumberOfStreamDivisions( mpi_size );

  if ( mpi_rank == 0 )
    {
    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cerr << "Exception caught !" << std::endl;
      std::cerr << excep << std::endl;
      }
    }

  // Region Splitter
  typedef itk::ImageRegionNonCubeSplitter< Dimension > SplitterType;
  SplitterType::Pointer splitter = SplitterType::New();

  // Get the splits
  std::vector < RegionType > split_regions( mpi_size );
  std::vector < RegionType > padded_regions( mpi_size );
  for ( int split = 0; split < mpi_size; ++split )
    {
    split_regions[ split ] = splitter->GetSplit( split, mpi_size, input_region );
    IndexType pad_index = split_regions[ split ].GetIndex();
    SizeType  pad_size  = split_regions[ split ].GetSize();
    // Pad the splits by one pixel in the positive direction.
    for ( unsigned int dim = 0; dim < Dimension; ++dim )
      {
      if (( pad_index[ dim ] + pad_size[ dim ] ) < ( input_index[ dim ] + input_size[ dim ] ))
        {
        ++pad_size[ dim ];
        }
      }
    padded_regions[ split ].SetIndex( pad_index );
    padded_regions[ split ].SetSize(  pad_size );
    }

  // We need to know who our neighbors are 
  // and the regions of overlap
  std::vector< RegionType > com_regions( mpi_size );
  IndexType my_index     = padded_regions[  mpi_rank ].GetIndex();
  SizeType  my_size      = padded_regions[  mpi_rank ].GetSize();
  for ( int split = 0; split < mpi_size; ++split )
    {
    IndexType neighbor_index = padded_regions[  split ].GetIndex();
    SizeType  neighbor_size  = padded_regions[  split ].GetSize();
    
    IndexType com_index;
    SizeType  com_size;
    
    for ( unsigned int dim = 0; dim < Dimension; ++dim )
      {
      
      // What is the INDEX and SIZE of the com region?
      if ( 
        ( my_index[ dim ]                       <= neighbor_index[ dim ] ) && 
        ( my_index[ dim ] + (int)my_size[ dim ] >  neighbor_index[ dim ] ))
        {
        com_index[ dim ] = neighbor_index[ dim ];
        com_size[  dim ] = my_index[ dim ] + my_size[ dim ] - neighbor_index[ dim ];
        }
      else if ( 
        ( neighbor_index[ dim ]                             <= my_index[ dim ] ) &&
        ( neighbor_index[ dim ] + (int)neighbor_size[ dim ] >  my_index[ dim ] ))
        {
        com_index[ dim ] = my_index[ dim ];
        com_size[  dim ] = neighbor_index[ dim ] + neighbor_size[ dim ] - my_index[ dim ];
        }
      else
        {
        com_index[ dim ] = 0;
        com_size[  dim ] = 0;
        }
            
      // We don't need to talk to ourself.
      if ( split == mpi_rank )
        {
        com_size[ dim ] = 0;
        }

      } // End dim loop
    com_regions[ split ].SetIndex( com_index );
    com_regions[ split ].SetSize(  com_size );
    } // End split loop

  // Extract our Piece
  typedef itk::ExtractImageFilter< FloatImageType, FloatImageType > ExtractorType;
  ExtractorType::Pointer extractor = ExtractorType::New();
  extractor->SetExtractionRegion( padded_regions[ mpi_rank ] );
  extractor->SetInput( reader->GetOutput() );

  // Fast Marching
  typedef   itk::FastMarchingImageFilter
    < FloatImageType, FloatImageType >  FastFilterType;
  FastFilterType::Pointer fast = FastFilterType::New();
  fast->SetInput( extractor->GetOutput() );

  // Create the node Container
  typedef FastFilterType::NodeContainer NodeContainer;
  typedef NodeContainer::ElementIdentifier ElementIdType;
  NodeContainer::Pointer seeds = NodeContainer::New();
  seeds->Initialize();
  ElementIdType num_seeds = 0;

  // Type for node
  typedef FastFilterType::NodeType NodeType;

  // Read the seeds from the file.
  FILE * seed_file = fopen( seed_file_name, "r" );
  unsigned int num_read = Dimension;
  while ( num_read == Dimension )
    {
    num_read = 0;
    PointType point;
    IndexType index;
    NodeType node;
    for ( unsigned int dim = 0; dim < Dimension; ++dim )
      {
      num_read += fscanf( seed_file, "%lf", &( point[ dim ] ));
      }
    input_image->TransformPhysicalPointToIndex( point, index );
    // Check if the seed is in our region.
    if ( split_regions[ mpi_rank ].IsInside( index ) == true )
      {
      node.SetValue( 0.0 );
      node.SetIndex( index );
      seeds->InsertElement( num_seeds, node );
      ++num_seeds;
      }
    }
  fclose( seed_file );

  // Add the seeds to the fast marching filter.
  fast->SetTrialPoints( seeds );

  // We need arrays of smart pointers to 
  // images for the send and receive buffers.
  std::vector< FloatImageType::Pointer > send_images( mpi_size );
  std::vector< FloatImageType::Pointer > recv_images( mpi_size );
  
  // We must create the receive images before we can get any data.
  for ( int split = 0; split < mpi_size; ++split )
    {
    if ( com_regions[ split ].GetNumberOfPixels() != 0 )
      {
      recv_images[ split ] = FloatImageType::New();
      recv_images[ split ]->SetOrigin( input_origin );
      recv_images[ split ]->SetRegions( com_regions[ split ] );
      recv_images[ split ]->Allocate();
      }
    }

  // We extract the send images from the fast marching filter.
  // Images created during iteration loop.
  
  // Iterator Types
  typedef itk::ImageRegionConstIterator< FloatImageType >  IteratorType;
  typedef itk::ConstNeighborhoodIterator< FloatImageType > NeighborhoodIteratorType;
  NeighborhoodIteratorType::RadiusType radius;
  radius.Fill(1);


  int iteration = 0;
  while ( iteration < max_iterations )
    {
    // Update fast marching now even though the send extractor will do it.
    try
      {
      fast->Update();
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cerr << "Exception caught while updating Fast Marching!" << std::endl;
      std::cerr << excep << std::endl;
      }

    for ( int split = 0; split < mpi_size; ++split )
      {
      if ( com_regions[ split ].GetNumberOfPixels() != 0 )
        {

        // Send Extractor
        ExtractorType::Pointer send_extractor = ExtractorType::New();
        send_extractor->SetInput( fast->GetOutput() );
        send_extractor->SetExtractionRegion( com_regions[ split ] );
        try
          {
          send_extractor->Update();
          }
        catch( itk::ExceptionObject & excep )
          {
          std::cerr << "Exception caught while updating Send Extractor!" << std::endl;
          std::cerr << excep << std::endl;
          }
        send_images[ split ] = send_extractor->GetOutput();
        send_images[ split ]->DisconnectPipeline();


        // MPI Communication
        // Low rank sends then receives 
        // while high rank receives then sends.
        if ( mpi_rank < split )
          {
          MPI_Send( send_images[ split ]->GetBufferPointer(), 
            com_regions[ split ].GetNumberOfPixels(),
            MPI_FLOAT, split, 0, MPI_COMM_WORLD );
          MPI_Recv( recv_images[ split ]->GetBufferPointer(),
            com_regions[ split ].GetNumberOfPixels(),
            MPI_FLOAT, split, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
          }
        else
          {
          MPI_Recv( recv_images[ split ]->GetBufferPointer(),
            com_regions[ split ].GetNumberOfPixels(),
            MPI_FLOAT, split, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
          MPI_Send( send_images[ split ]->GetBufferPointer(), 
            com_regions[ split ].GetNumberOfPixels(),
            MPI_FLOAT, split, 0, MPI_COMM_WORLD );
          }

        } // End neighbors
      } // End split loop

    // Trial Points
    NodeContainer::Pointer new_trial_points = NodeContainer::New();
    new_trial_points->Initialize();
    ElementIdType num_new_trial_points = 0;

    // This will hold the minimum value at our boundary.
    double min_boundary_time = itk::NumericTraits< FloatPixelType >::max();


    // Take another pass over all of the neighbor splits and set the trial points.
    for ( int split = 0; split < mpi_size; ++split )
      {
      if ( com_regions[ split ].GetNumberOfPixels() != 0 )
        {
        IteratorType recv_iterator( recv_images[ split ], com_regions[ split ] );
        IteratorType send_iterator( send_images[ split ], com_regions[ split ] );
        for ( recv_iterator.GoToBegin(), send_iterator.GoToBegin(); 
          !recv_iterator.IsAtEnd(); 
          ++recv_iterator, ++send_iterator )
          {
          FloatImageType::PixelType recv_value = recv_iterator.Get();
          FloatImageType::PixelType send_value = send_iterator.Get();
          if ( recv_value < send_value )
            {
            IndexType index = recv_iterator.GetIndex();
            NodeType node;
            node.SetValue( recv_value );
            node.SetIndex( index );
            new_trial_points->InsertElement( num_new_trial_points, node );
            ++num_new_trial_points;
            if ( recv_value < min_boundary_time )
              {
              min_boundary_time = recv_value;
              }
            }
          } // End Com iterator
        } // End neighbors
      } // End split loop


    // Save the number of trial points we got from our edge.
    ElementIdType edge_trial_points = num_new_trial_points;

    // Alive Points
    NodeContainer::Pointer alive_points = NodeContainer::New();
    alive_points->Initialize();
    ElementIdType num_alive_points = 0;

    // We'll iterate over our Fast Marching output and set alive and more trial points.
    NeighborhoodIteratorType trial_alive_iterator( 
      radius, fast->GetOutput(), split_regions[ mpi_rank ] );
    for ( 
      trial_alive_iterator.GoToBegin(); 
      !trial_alive_iterator.IsAtEnd(); 
      ++trial_alive_iterator )
      {
      FloatImageType::PixelType value = trial_alive_iterator.GetCenterPixel();
      // Is this pixel an alive or trial point?
      if ( value < min_boundary_time )
        {
        NodeType node;
        IndexType index = trial_alive_iterator.GetIndex();
        unsigned long neighborhood_size = trial_alive_iterator.Size();
        int neighbor_above = 0;
        node.SetValue( value );
        node.SetIndex( index );
        for ( unsigned long neighbor = 0; neighbor < neighborhood_size; ++neighbor )
          {
          // Are any of this pixels neighbors above threshold?
          if ( trial_alive_iterator.GetPixel( neighbor ) > min_boundary_time )
            {
            neighbor_above = 1;
            }
          } // End for neighbors
        if ( neighbor_above == 1 )
          {
          // This pixel is a trial point!
          new_trial_points->InsertElement( num_new_trial_points, node );
          ++num_new_trial_points;
          }
        else
          {
          // This pixel is an alive point.
          alive_points->InsertElement( num_alive_points, node );
          ++num_alive_points;
          }
        } // End if value below
      } // End for Trial versus Alive Iterator


    std::cout 
      << "Iteration " << iteration 
      << " rank " << mpi_rank 
      << " has " << num_alive_points << " alive points, and "
      << num_new_trial_points << " trial points, " 
      << edge_trial_points << " of wich came from an edge." << std::endl;

    // Set the new trial and alive points for the fast marching filter.
    fast->SetTrialPoints( new_trial_points );
    fast->SetAlivePoints( alive_points );

    ++iteration;
    } // End itteration loop

  // This extractor removes the padding around our region.
  ExtractorType::Pointer final_extractor = ExtractorType::New();
  final_extractor->SetInput( fast->GetOutput() );
  final_extractor->SetExtractionRegion( split_regions[ mpi_rank ] );

  IndexType write_index = split_regions[ mpi_rank ].GetIndex();
  SizeType  write_size  = split_regions[ mpi_rank ].GetSize();

  // We paste our output into the input.
  typedef itk::PasteImageFilter < FloatImageType, FloatImageType, FloatImageType > PasteType;
  PasteType::Pointer paste = PasteType::New();
  paste->SetInput( 0, reader->GetOutput() );
  paste->SetInput( 1, final_extractor->GetOutput() );
  paste->SetDestinationIndex( write_index );
  paste->SetSourceRegion( split_regions[ mpi_rank ] );

  // Writer writes only the region we're responsible for.
  // The other regions will be written by the other MPI processes.
  itk::ImageIORegion write_region( Dimension );
  for ( unsigned int dim = 0; dim < Dimension; ++dim )
    {
    write_region.SetIndex( dim, write_index[ dim ] );
    write_region.SetSize(  dim, write_size[  dim ] );
    }
   writer->SetIORegion( write_region );
   writer->SetInput( paste->GetOutput() );


  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }
  
  MPI_Finalize();

  return EXIT_SUCCESS;
}
