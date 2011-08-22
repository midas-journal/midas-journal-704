
#ifndef _itkImageRegionNonCubeSplitter_txx
#define _itkImageRegionNonCubeSplitter_txx
#include "itkImageRegionNonCubeSplitter.h"
#include <math.h>

namespace itk
{

/**
 *
 */
template <unsigned int VImageDimension>
unsigned int 
ImageRegionNonCubeSplitter<VImageDimension>
::GetNumberOfSplits(const RegionType &region, unsigned int requestedNumber)
{
  const SizeType &regionSize = region.GetSize();

  // if a given region has fewer pixels than requestedNumber, then
  // only split number of pixels times
  unsigned int i, numPixels;
  numPixels = 1;
  for (i=0; i < VImageDimension; ++i)
    numPixels *= regionSize[i];
  if (numPixels > requestedNumber)
    return requestedNumber;

  return numPixels;
}

  
/**
 *
 */
template <unsigned int VImageDimension>
ImageRegion<VImageDimension>
ImageRegionNonCubeSplitter<VImageDimension>
::GetSplit(unsigned int i, unsigned int numberOfPieces,
           const RegionType &region)
{
  RegionType splitRegion;
  IndexType splitIndex;
  SizeType splitSize, regionSize;
  
  // Initialize the splitRegion to the requested region
  splitRegion = region;
  splitIndex = splitRegion.GetIndex();
  splitSize = splitRegion.GetSize();

  regionSize = region.GetSize();
  
  // Keep the number of pieces from going over the 
  // number of pixels.
  // The user will probably be surprised by this.
  unsigned long numPixels;
  numPixels = 1;
  for ( unsigned long dim = 0;
    dim < VImageDimension; 
    ++dim )
    numPixels *= regionSize[dim];
  if ( numberOfPieces > numPixels )
    numberOfPieces = numPixels;

  // I'm not as smart as whoever wrote the VTK splitter.
  // I'm going to make a big array of sizes and indecies 
  // and search through it a lot. I'll justify it by 
  // saying that hopefully I'll need all of the splits 
  // at once.

  // For some reason I can't make an array of RegionType, IndexType, 
  // or even SizeType::SizeValueType.
  typedef itk::Array< unsigned long > ArrayOfSizeType;
  ArrayOfSizeType sizes( numberOfPieces * VImageDimension );

  typedef itk::Array< unsigned long > ArrayOfIndexType;
  ArrayOfIndexType indecies( numberOfPieces * VImageDimension );

  unsigned long num_pieces_so_far = 0;

  // Set the first piece.
  for (unsigned long dim = 0; dim < VImageDimension; ++dim)
    {
    sizes[ num_pieces_so_far * VImageDimension + dim ] 
      = region.GetSize()[ dim ];

    indecies[ num_pieces_so_far * VImageDimension + dim ] 
      = region.GetIndex()[ dim ];
    }
  ++num_pieces_so_far;

  // Split untill we've made the right number.
  while ( num_pieces_so_far < numberOfPieces )
    {
    // Find the longest dim.
    unsigned long long_dim = 0;
    unsigned long long_dim_piece = 0;
    unsigned long long_dim_dim = 0;
    // Loop over image dim
    for ( unsigned long dim = 0;
      dim < VImageDimension;
      ++dim )
      {
      // Loop over piece
      for ( unsigned long piece = 0;
        piece < num_pieces_so_far;
        ++piece )
        {
        // Check the dim size
        if ( sizes[piece * VImageDimension + dim] > long_dim )
          {
          long_dim = sizes[piece * VImageDimension + dim];
          long_dim_piece = piece;
          long_dim_dim = dim;
          } // End if long_dim
        } // End for piece
      } // End for dim
    // Copy the longest piece
    for ( unsigned long dim = 0;
      dim < VImageDimension;
      ++dim )
      {
      sizes[num_pieces_so_far * VImageDimension + dim] 
        = sizes[long_dim_piece * VImageDimension + dim];
      indecies[num_pieces_so_far * VImageDimension + dim] 
        = indecies[long_dim_piece * VImageDimension + dim];
      }
    // Split the longest dim of the copy
    sizes[num_pieces_so_far * VImageDimension + long_dim_dim] 
      -= sizes[long_dim_piece * VImageDimension + long_dim_dim] / 2;
    indecies[num_pieces_so_far * VImageDimension + long_dim_dim] 
      += sizes[long_dim_piece * VImageDimension + long_dim_dim] / 2;
    // Split the origional piece
    sizes[long_dim_piece * VImageDimension + long_dim_dim] /= 2;
    ++num_pieces_so_far;
    } // End while num_pieces_so_far

  // Copy out the piece
  for ( unsigned long dim = 0; dim < VImageDimension; ++dim )
    {
    splitIndex[dim] = indecies[i * VImageDimension + dim];
    splitSize[dim] = sizes[i * VImageDimension + dim];
    }
  // set the split region ivars
  splitRegion.SetIndex( splitIndex );
  splitRegion.SetSize( splitSize );

  itkDebugMacro("  Split Piece: " << std::endl << splitRegion );

  return splitRegion;
}
  
  

/**
 *
 */
template <unsigned int VImageDimension>
void 
ImageRegionNonCubeSplitter<VImageDimension>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}


} // end namespace itk

#endif
