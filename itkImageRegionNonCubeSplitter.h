
#ifndef __itkImageRegionNonCubeSplitter_h
#define __itkImageRegionNonCubeSplitter_h

#include "itkObject.h"
#include "itkRegion.h"
#include "itkObjectFactory.h"
#include "itkIndex.h"
#include "itkSize.h"
#include "itkArray.h"
#include "itkImageRegionSplitter.h"

namespace itk
{

/** \class ImageRegionNonCubeSplitter
 * \brief Divide a region into several pieces.
 *
 * ImageRegionNonCubeSplitter divides an ImageRegion into
 * smaller regions.  ImageRegionNonCubeSplitter is meant to split the 
 * region into exactly N smaller regions.  This object has two
 * basic methods: GetNumberOfSplits() and GetSplit().
 * 
 * GetNumberOfSplits() is used to determine how may subregions a given
 * region can be divided.  You call GetNumberOfSplits with an argument
 * that is the number of subregions you want.  If the image region can
 * support that number of subregions, that number is returned.
 * Otherwise, the maximum number of splits a region can support will
 * be returned.  Here the maximum number of splits is the number of 
 * pixels
 *
 * GetSplit() returns the ith of N subregions (as an ImageRegion object).
 *
 *
 * \ingroup ITKSystemObjects
 * \ingroup DataProcessing
 */

template <unsigned int VImageDimension>
class ITK_EXPORT ImageRegionNonCubeSplitter: public ImageRegionSplitter<VImageDimension>
{
public:
  /** Standard class typedefs. */
  typedef ImageRegionNonCubeSplitter              Self;
  typedef ImageRegionSplitter<VImageDimension>  Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageRegionNonCubeSplitter,ImageRegionSplitter);

  /** Dimension of the image available at compile time. */
  itkStaticConstMacro(ImageDimension, unsigned int, VImageDimension);
  
  /** Index typedef support. An index is used to access pixel values. */
  typedef Index<VImageDimension>  IndexType;

  /** Size typedef support. A size is used to define region bounds. */
  typedef Size<VImageDimension>  SizeType;

  /** Region typedef support.   */
  typedef ImageRegion<VImageDimension> RegionType;

  /** How many pieces can the specifed region be split? A given region
   * cannot always be divided into the requested number of pieces.  
   * This method returns a number less than or equal to the requested 
   * number of pieces.  */
  virtual unsigned int GetNumberOfSplits(const RegionType &region,
                                         unsigned int requestedNumber);

  /** Get a region definition that represents the ith piece a specified region.
   * The "numberOfPieces" specified should be less than or equal to what
   * GetNumberOfSplits() returns. */
  virtual RegionType GetSplit(unsigned int i, unsigned int numberOfPieces,
                              const RegionType &region);

protected:
  ImageRegionNonCubeSplitter() {}
  ~ImageRegionNonCubeSplitter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  ImageRegionNonCubeSplitter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_ImageRegionNonCubeSplitter(_, EXPORT, x, y) namespace itk { \
  _(1(class EXPORT ImageRegionNonCubeSplitter< ITK_TEMPLATE_1 x >)) \
  namespace Templates { typedef ImageRegionNonCubeSplitter< ITK_TEMPLATE_1 x > \
                                               ImageRegionNonCubeSplitter##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkImageRegionNonCubeSplitter+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkImageRegionNonCubeSplitter.txx"
#endif

#endif

