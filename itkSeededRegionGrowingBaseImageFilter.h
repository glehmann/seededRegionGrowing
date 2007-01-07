#ifndef __itkSeededRegionGrowingBaseImageFilter_h
#define __itkSeededRegionGrowingBaseImageFilter_h

#include "itkImageToImageFilter.h"

namespace itk {

/** \class SeededRegionGrowingBaseImageFilter
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
template<class TInputImage, class TLabelImage, class TRegionStats>
class ITK_EXPORT SeededRegionGrowingBaseImageFilter : 
    public ImageToImageFilter<TInputImage, TLabelImage>
{
public:
  /** Standard class typedefs. */
  typedef SeededRegionGrowingBaseImageFilter Self;
  typedef ImageToImageFilter<TInputImage, TLabelImage>
  Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage InputImageType;
  typedef TLabelImage LabelImageType;
  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename InputImageType::RegionType     InputImageRegionType;
  typedef typename InputImageType::PixelType      InputImagePixelType;
  typedef typename LabelImageType::Pointer        LabelImagePointer;
  typedef typename LabelImageType::ConstPointer   LabelImageConstPointer;
  typedef typename LabelImageType::RegionType     LabelImageRegionType;
  typedef typename LabelImageType::PixelType      LabelImagePixelType;
  
  typedef typename LabelImageType::IndexType      IndexType;
  typedef typename NumericTraits<InputImagePixelType>::RealType RealType;
  /** ImageDimension constants */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Standard New method. */
  itkNewMacro(Self);  

  /** Runtime information support. */
  itkTypeMacro(SeededRegionGrowingBaseImageFilter, 
               ImageToImageFilter);
  

   /** Set the marker image */
  void SetMarkerImage(TLabelImage *input)
     {
     // Process object is not const-correct so the const casting is required.
     this->SetNthInput( 1, const_cast<TLabelImage *>(input) );
     }

  /** Get the marker image */
  LabelImageType * GetMarkerImage()
    {
    return static_cast<LabelImageType*>(const_cast<DataObject *>(this->ProcessObject::GetInput(1)));
    }

   /** Set the input image */
  void SetInput1(TInputImage *input)
     {
     this->SetInput( input );
     }

   /** Set the marker image */
  void SetInput2(TLabelImage *input)
     {
     this->SetMarkerImage( input );
     }

  /**
   * Set/Get whether the connected components are defined strictly by
   * face connectivity or by face+edge+vertex connectivity.  Default is
   * FullyConnectedOff.  For objects that are 1 pixel wide, use
   * FullyConnectedOn.
   */
  itkSetMacro(FullyConnected, bool);
  itkGetConstReferenceMacro(FullyConnected, bool);
  itkBooleanMacro(FullyConnected);

  /**
   * Set/Get whether the border between regions will be marked or not. Default
   * is true.
   */
  itkSetMacro(MarkBoundaryLine, bool);
  itkGetConstReferenceMacro(MarkBoundaryLine, bool);
  itkBooleanMacro(MarkBoundaryLine);

  itkSetMacro(PadImageBoundary, bool);
  itkGetConstReferenceMacro(PadImageBoundary, bool);
  itkBooleanMacro(PadImageBoundary);

protected:
  SeededRegionGrowingBaseImageFilter();
  ~SeededRegionGrowingBaseImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** SeededRegionGrowingBaseImageFilter needs to request enough of the
   * marker image to account for the elementary structuring element.
   * The mask image does not need to be padded. Depending on whether
   * the filter is configured to run a single iteration or until
   * convergence, this method may request all of the marker and mask
   * image be provided. */
  void GenerateInputRequestedRegion();

  /** This filter will enlarge the output requested region to produce
   * all of the output if the filter is configured to run to
   * convergence.
   * \sa ProcessObject::EnlargeOutputRequestedRegion() */
  void EnlargeOutputRequestedRegion(DataObject *itkNotUsed(output));

  /** Single-threaded version of GenerateData. */
  void GenerateData();
  
  /** Data types to support the region statistics */
  typedef TRegionStats StatsMapType;

  
private:
  SeededRegionGrowingBaseImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  bool m_FullyConnected;

  bool m_MarkBoundaryLine;
  bool m_PadImageBoundary;


  void InitMap(StatsMapType &SM, LabelImageConstPointer LabelIm, 
	       InputImageConstPointer InputIm);
} ; // end of class

} // end namespace itk
  

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSeededRegionGrowingBaseImageFilter.txx"
#endif

#endif


