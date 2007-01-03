#ifndef __itkSeededRegionGrowingImageFilter_h
#define __itkSeededRegionGrowingImageFilter_h

#include "itkSeededRegionGrowingBaseImageFilter.h"
#include "itkSeededRegionGrowingSupport.h"

namespace itk {

/** \class SeededRegionGrowingImageFilter
 * \brief Greedy region growing from markers
 * This filter implements a greedy seeded region growing algorithm
 * which results in all pixels being assigned to a class based on
 * similarity to the seed regions and connectivity. There is no need
 * to define a stopping point, as with the "traditional" region
 * growing approaches. However there must be enough seed regions in
 * the right places to achieve what you want. In many ways this is
 * similar to the morphological watershed filter, and much of the
 * implementation is borrowed from that filter.
 *
 * The result is of type TLabelImage and will have the same labels as
 * the input. Disconnected seeds with the same label will be treated
 * as the same region.
 *
 * Seeded region growing can be applied directly to multichannel
 * images (as soon as I figure out how to code that sort of thing) and
 * can be easier to use with strongly anisotropic pixels.
 *
 * "Seeded region growing", R. Adams and L. Bischof, PAMI, 1994.
 *
 * \author Richard Beare, Richard.Beare@med.monash.edu.au
 *
 * \sa WatershedImageFilter, MorphologicalWatershedImageFilter
 * \ingroup ImageEnhancement  MathematicalMorphologyImageFilters
 */
template<class TInputImage, class TLabelImage>
class ITK_EXPORT SeededRegionGrowingImageFilter : 
    public SeededRegionGrowingBaseImageFilter<TInputImage, TLabelImage, 
					      SRGRegionStatsType<typename TLabelImage::PixelType,
								 typename TInputImage::PixelType> >
{
public:
  /** Standard class typedefs. */
  typedef SeededRegionGrowingImageFilter Self;
  typedef SeededRegionGrowingBaseImageFilter<TInputImage, TLabelImage, 
    SRGRegionStatsType<typename TLabelImage::PixelType,
    typename TInputImage::PixelType> >  Superclass;

  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

 /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(SeededRegionGrowingImageFilter,
               SeededRegionGrowingBaseImageFilter);

 
protected:
  SeededRegionGrowingImageFilter();
  ~SeededRegionGrowingImageFilter() {};

private:
  SeededRegionGrowingImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

} ; // end of class


} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSeededRegionGrowingImageFilter.txx"
#endif

#endif


