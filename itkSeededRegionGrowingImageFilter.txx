#ifndef __itkSeededRegionGrowingImageFilter_txx
#define __itkSeededRegionGrowingImageFilter_txx

#include "itkSeededRegionGrowingImageFilter.h"
#include "itkProgressReporter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkSize.h"
#include "itkConnectedComponentAlgorithm.h"
#include "itkPriorityQueue.h"
#include "itkImageDuplicator.h"
namespace itk {

template <class TInputImage, class TLabelImage>
SeededRegionGrowingImageFilter<TInputImage, TLabelImage>
::SeededRegionGrowingImageFilter()
{
  this->SetNumberOfRequiredInputs(2);
  m_FullyConnected = false;
  m_MarkBoundaryLine = true;
  m_PadImageBoundary = true;
}


template <class TInputImage, class TLabelImage>
void 
SeededRegionGrowingImageFilter<TInputImage, TLabelImage>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  
  // get pointers to the inputs
  LabelImagePointer  markerPtr = this->GetMarkerImage();

  InputImagePointer  inputPtr = 
    const_cast< InputImageType * >( this->GetInput() );
  
  if ( !markerPtr || !inputPtr )
    { return; }

  // We need to
  // configure the inputs such that all the data is available.
  //
  markerPtr->SetRequestedRegion(markerPtr->GetLargestPossibleRegion());
  inputPtr->SetRequestedRegion(inputPtr->GetLargestPossibleRegion());
}


template <class TInputImage, class TLabelImage>
void 
SeededRegionGrowingImageFilter<TInputImage, TLabelImage>
::EnlargeOutputRequestedRegion(DataObject *)
{
  this->GetOutput()->SetRequestedRegion( this->GetOutput()->GetLargestPossibleRegion() );
}


template<class TInputImage, class TLabelImage>
void
SeededRegionGrowingImageFilter<TInputImage, TLabelImage>
::GenerateData()
{
  // there are 2 possible cases: with or without boundary lines.

  // the label used to find background in the marker image
  static const LabelImagePixelType bgLabel = NumericTraits< LabelImagePixelType >::Zero;
  // the label used to mark the watershed line in the output image
  static const LabelImagePixelType wsLabel = NumericTraits< LabelImagePixelType >::Zero;

  this->AllocateOutputs();
  // Set up the progress reporter
  // we can't found the exact number of pixel to process in the 2nd pass, so we use the maximum number possible.
  ProgressReporter progress(this, 0, this->GetMarkerImage()->GetRequestedRegion().GetNumberOfPixels()*2);
  
  // mask and marker must have the same size
  if ( this->GetMarkerImage()->GetRequestedRegion().GetSize() != this->GetInput()->GetRequestedRegion().GetSize() )
    { itkExceptionMacro( << "Marker and input must have the same size." ); }
  

  // the radius which will be used for all the shaped iterators
  Size< ImageDimension > radius;
  radius.Fill(1);

  // iterator for the marker image
  typedef ConstShapedNeighborhoodIterator<LabelImageType> MarkerIteratorType;
  typename MarkerIteratorType::ConstIterator nmIt;
  MarkerIteratorType markerIt(radius, this->GetMarkerImage(), this->GetMarkerImage()->GetRequestedRegion());
  // add a boundary constant to avoid adding pixels on the border in the fah
  ConstantBoundaryCondition<LabelImageType> lcbc;
  lcbc.SetConstant( NumericTraits<LabelImagePixelType>::max() );
  markerIt.OverrideBoundaryCondition(&lcbc);
  setConnectivity( &markerIt, m_FullyConnected );

  // iterator for the input image
  typedef ConstShapedNeighborhoodIterator<InputImageType> InputIteratorType;
  InputIteratorType inputIt(radius, this->GetInput(), this->GetInput()->GetRequestedRegion());
  typename InputIteratorType::ConstIterator niIt;
  setConnectivity( &inputIt, m_FullyConnected );
  
  // iterator for the output image
  typedef ShapedNeighborhoodIterator<LabelImageType> OutputIteratorType;
  typedef typename OutputIteratorType::OffsetType OffsetType;
  typename OutputIteratorType::Iterator noIt;
  OutputIteratorType outputIt(radius, this->GetOutput(), this->GetOutput()->GetRequestedRegion());
  setConnectivity( &outputIt, m_FullyConnected );

  // initialization - iterate over seed and control images to compute
  // the mean of each region

  // FAH (in french: File d'Attente Hierarchique)
  typedef PriorityQueue< InputImagePixelType, IndexType > PriorityQueueType;
  PriorityQueueType fah;
  //---------------------------------------------------------------------------
  // Meyer's algorithm
  //---------------------------------------------------------------------------
  if( m_MarkBoundaryLine )
    {
    // first stage:
    //  - set markers pixels to already processed status
    //  - copy markers pixels to output image
    //  - init FAH with indexes of background pixels with marker pixel(s) in their neighborhood
    
    ConstantBoundaryCondition<LabelImageType> lcbc2;
    lcbc2.SetConstant( wsLabel ); // outside pixel are watershed so they won't be use to find real watershed  pixels
    outputIt.OverrideBoundaryCondition(&lcbc2);
    
    // create a temporary image to store the state of each pixel (processed or not)
    typedef Image< bool, ImageDimension > StatusImageType;
    typename StatusImageType::Pointer statusImage = StatusImageType::New();
    statusImage->SetRegions( this->GetMarkerImage()->GetLargestPossibleRegion() );
    statusImage->Allocate();
  
    // iterator for the status image
    typedef ShapedNeighborhoodIterator<StatusImageType> StatusIteratorType;
    typename StatusIteratorType::Iterator nsIt;
    StatusIteratorType statusIt(radius, statusImage, this->GetOutput()->GetRequestedRegion());
    ConstantBoundaryCondition< StatusImageType > bcbc;
    bcbc.SetConstant( true );  // outside pixel are already processed
    statusIt.OverrideBoundaryCondition(&bcbc);
    setConnectivity( &statusIt, m_FullyConnected );
    
    // the status image must be initialized before the first stage. In the first stage, the
    // set to true are the neighbors of the marker (and the marker) so it's difficult
    // (impossible ?)to init the status image at the same time
    // the overhead should be small
    statusImage->FillBuffer( false );
    
    for ( markerIt.GoToBegin(), statusIt.GoToBegin(), outputIt.GoToBegin(), inputIt.GoToBegin();
	  !markerIt.IsAtEnd(); ++markerIt, ++outputIt)
      {
      LabelImagePixelType markerPixel = markerIt.GetCenterPixel();
      if ( markerPixel != bgLabel )
	{
        
	IndexType idx = markerIt.GetIndex();
        
	// move the iterators to the right place
	OffsetType shift = idx - statusIt.GetIndex();
	statusIt += shift;
	inputIt += shift;
        
	// this pixel belongs to a marker
	// mark it as already processed
	statusIt.SetCenterPixel( true );
	// copy it to the output image
	outputIt.SetCenterPixel( markerPixel );
	// and increase progress because this pixel will not be used in the flooding stage.
	progress.CompletedPixel();
        
	// search the background pixels in the neighborhood
	for ( nmIt= markerIt.Begin(), nsIt= statusIt.Begin(), niIt= inputIt.Begin();
	      nmIt != markerIt.End();
	      nmIt++, nsIt++, niIt++ )
	  {
	  if ( !nsIt.Get() && nmIt.Get() == bgLabel )
	    {
	    // this neighbor is a background pixel and is not already processed; add its
	    // index to fah
	    fah.Push( niIt.Get(), markerIt.GetIndex() + nmIt.GetNeighborhoodOffset() );
	    // mark it as already in the fah to avoid adding it several times
	    nsIt.Set( true );
	    }
	  }
	}
      else
	{
	// Some pixels may be never processed so, by default, non marked pixels
	// must be marked as watershed
	outputIt.SetCenterPixel( wsLabel );
	}
      // one more pixel done in the init stage
      progress.CompletedPixel();
      }
    // end of init stage
    
    // flooding
    // init all the iterators
    outputIt.GoToBegin();
    statusIt.GoToBegin();
    inputIt.GoToBegin();
    
    // and start flooding
    while( !fah.Empty() )
      {
      // store the current vars
      const InputImagePixelType & currentValue = fah.FrontKey();
      const IndexType & idx = fah.FrontValue();
      
      // move the iterators to the right place
      OffsetType shift = idx - outputIt.GetIndex();
      outputIt += shift;
      statusIt += shift;
      inputIt += shift;
      
      // iterate over the neighbors. If there is only one marker value, give that value
      // to the pixel, else keep it as is (watershed line)
      LabelImagePixelType marker = wsLabel;
      bool collision = false;
      for (noIt = outputIt.Begin(); noIt != outputIt.End(); noIt++)
	{
	LabelImagePixelType o = noIt.Get();
	if( o != wsLabel )
	  {
	  if( marker != wsLabel && o != marker )
	    { 
	    collision = true; 
	    break;
	    }
	  else
	    { marker = o; }
	  }
	}
      if( !collision )
	{
	// set the marker value
	outputIt.SetCenterPixel( marker );
	// and propagate to the neighbors
	for (niIt = inputIt.Begin(), nsIt = statusIt.Begin();
	     niIt != inputIt.End();
	     niIt++, nsIt++)
	  {
	  if ( !nsIt.Get() )
	    {
	    // the pixel is not yet processed. add it to the fah
	    InputImagePixelType grayVal = niIt.Get();
	    if ( grayVal <= currentValue )
	      { fah.Push( currentValue, inputIt.GetIndex() + niIt.GetNeighborhoodOffset() ); }
	    else
	      { fah.Push( grayVal, inputIt.GetIndex() + niIt.GetNeighborhoodOffset() ); }
	    // mark it as already in the fah
	    nsIt.Set( true );
	    }
	  }
	}
      // one more pixel in the flooding stage
      progress.CompletedPixel();
      // remove the processed pixel of the queue
      fah.Pop();
      }
    }


  //---------------------------------------------------------------------------
  // Beucher's algorithm
  //---------------------------------------------------------------------------
  else
    {
    // first stage:
    //  - copy markers pixels to output image
    //  - init FAH with indexes of pixels with background pixel in their neighborhood
    
    ConstantBoundaryCondition<LabelImageType> lcbc2;
    lcbc2.SetConstant( NumericTraits< LabelImagePixelType >::max() ); // outside pixel are watershed so they won't be use to find real watershed  pixels
    outputIt.OverrideBoundaryCondition(&lcbc2);
    
    for ( markerIt.GoToBegin(), outputIt.GoToBegin(), inputIt.GoToBegin();
	  !markerIt.IsAtEnd();
	  ++markerIt, ++outputIt ) 
      {
      LabelImagePixelType markerPixel = markerIt.GetCenterPixel();
      if ( markerPixel != bgLabel )
	{
	IndexType idx = markerIt.GetIndex();
	OffsetType shift = idx - inputIt.GetIndex();
	inputIt += shift;
        
	// this pixels belongs to a marker
	// copy it to the output image
	outputIt.SetCenterPixel( markerPixel );
	// search if it has background pixel in its neighborhood
	bool haveBgNeighbor = false;
	for ( nmIt= markerIt.Begin(); nmIt != markerIt.End(); nmIt++ )
	  {
	  if ( nmIt.Get() == bgLabel )
	    { 
	    haveBgNeighbor = true; 
	    break;
	    }
	  }
	if ( haveBgNeighbor )
	  {
	  // there is a background pixel in the neighborhood; add to fah
	  fah.Push( inputIt.GetCenterPixel(), markerIt.GetIndex() );
	  }
	else
	  {
	  // increase progress because this pixel will not be used in the flooding stage.
	  progress.CompletedPixel();
	  }
	}
      else
	{
	outputIt.SetCenterPixel( wsLabel );
	}
      progress.CompletedPixel();
      }
    // end of init stage
      
    // flooding
    // init all the iterators
    outputIt.GoToBegin();
    inputIt.GoToBegin();
    
    // and start flooding
    while( !fah.Empty() )
      {
      // store the current vars
      const InputImagePixelType & currentValue = fah.FrontKey();
      const IndexType & idx = fah.FrontValue();
      
      // move the iterators to the right place
      OffsetType shift = idx - outputIt.GetIndex();
      outputIt += shift;
      inputIt += shift;
      
      LabelImagePixelType currentMarker = outputIt.GetCenterPixel();
      // get the current value of the pixel
      // iterate over neighbors to propagate the marker
      for (noIt = outputIt.Begin(), niIt = inputIt.Begin();
	   noIt != outputIt.End();
	   noIt++, niIt++)
	{
	if ( noIt.Get() == wsLabel )
	  {
	  // the pixel is not yet processed. It can be labeled with the current label
	  noIt.Set( currentMarker );
	  const InputImagePixelType & grayVal = niIt.Get();
	  if ( grayVal <= currentValue )
	    { fah.Push( currentValue, inputIt.GetIndex() + niIt.GetNeighborhoodOffset() ); }
            else
              { fah.Push( grayVal, inputIt.GetIndex() + niIt.GetNeighborhoodOffset() ); }
	  progress.CompletedPixel();
	  }
	}
      // remove the processed pixel of the queue
      fah.Pop();
      }
    }
  
}



template<class TInputImage, class TLabelImage>
void
SeededRegionGrowingImageFilter<TInputImage, TLabelImage>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  
  os << indent << "FullyConnected: "  << m_FullyConnected << std::endl;
  os << indent << "MarkBoundaryLine: "  << m_MarkBoundaryLine << std::endl;
}
  
}// end namespace itk
#endif
