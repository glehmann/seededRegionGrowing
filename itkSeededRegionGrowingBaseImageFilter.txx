#ifndef __itkSeededRegionGrowingBaseImageFilter_txx
#define __itkSeededRegionGrowingBaseImageFilter_txx

#include "itkSeededRegionGrowingBaseImageFilter.h"
#include "itkProgressReporter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkSize.h"
#include "itkConnectedComponentAlgorithm.h"
//#include "itkPriorityQueue.h"
#include "itkImageDuplicator.h"

#include <queue>
#include <map>

namespace itk {

template <class TInputImage,  class TLabelImage, class TRegionStats>
SeededRegionGrowingBaseImageFilter<TInputImage, TLabelImage, TRegionStats>
::SeededRegionGrowingBaseImageFilter()
{
  this->SetNumberOfRequiredInputs(2);
  m_FullyConnected = false;
  m_MarkBoundaryLine = true;
}


template <class TInputImage,  class TLabelImage, class TRegionStats>
void 
SeededRegionGrowingBaseImageFilter<TInputImage, TLabelImage, TRegionStats>
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


template <class TInputImage,  class TLabelImage, class TRegionStats>
void 
SeededRegionGrowingBaseImageFilter<TInputImage, TLabelImage, TRegionStats>
::EnlargeOutputRequestedRegion(DataObject *)
{
  this->GetOutput()->SetRequestedRegion( this->GetOutput()->GetLargestPossibleRegion() );
}


template<class TInputImage,  class TLabelImage, class TRegionStats>
void
SeededRegionGrowingBaseImageFilter<TInputImage, TLabelImage, TRegionStats>
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
  
	typedef typename InputImageType::OffsetValueType OffScalar;
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
  StatsMapType StatsMap;
  this->InitMap(StatsMap, this->GetMarkerImage(), this->GetInput());

  // FAH (in french: File d'Attente Hierarchique)
  //typedef PriorityQueue< RealType, IndexType > PriorityQueueType;
//	typedef std::queue< IndexType >                    QueueType;
	typedef std::queue< OffScalar >                    QueueType;
  // typedef std::queue< IndexType, std::list<IndexType> >                    QueueType;
  typedef std::map< RealType, QueueType > MapType;
  MapType fah;

  //PriorityQueueType fah;
  //---------------------------------------------------------------------------
  // based on Meyer's algorithm for watershed transform
  //---------------------------------------------------------------------------
  //if( m_MarkBoundaryLine )
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
	    // add neighbors to the queue if they are background, and
	    // haven't been added to the queue already
	    RealType priority=StatsMap.computePriority(markerPixel, niIt.Get());
	    priority = (RealType)(InputImagePixelType)round(priority);
	    //fah.Push(priority, markerIt.GetIndex() +
	    //nmIt.GetNeighborhoodOffset() );
	    fah[priority].push(this->GetMarkerImage()->ComputeOffset(markerIt.GetIndex() + nmIt.GetNeighborhoodOffset()));
	    // mark it as already in the fah to avoid adding it
	    // several times  
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
    while( !fah.empty() )
      {
      // store the current vars
      //const IndexType & idx = fah.FrontValue();
      // remove the processed pixel of the queue -- different to the
      // watershed because new pixels may get added with higher
      // priority than the current one
      //fah.Pop();

//       size_t totalsize  = 0;
//       for (typename MapType::iterator mit = fah.begin(); mit != fah.end(); mit++)
// 	{
// 	totalsize += mit->second.size();
// 	}
// //      std::cout << maxmapsize << std::endl;
//       maxmapsize = std::max(maxmapsize, totalsize);

      //RealType currentValue = fah.begin()->first;
      //QueueType currentQueue = fah.begin()->second;
      // nasty - if this queue is empty, delete it and go back to the
      // beginning of the loop
      if (fah.begin()->second.empty()) 
	{
	// will this force disposal?
	//QueueType currentQueue = fah.begin()->second;
	fah.erase(fah.begin());
	continue;
	}

//       IndexType idx = currentQueue.front();
//       currentQueue.pop();
      IndexType idx = this->GetMarkerImage()->ComputeIndex(fah.begin()->second.front());
      fah.begin()->second.pop();
      
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
	    if( !m_MarkBoundaryLine )
	      {
	      // need to do something a bit dodgy - check which of the
	      // markers in the neigborhood is better
	      InputImagePixelType grayVal = inputIt.GetCenterPixel();
	      RealType p1=StatsMap.computePriority(marker, grayVal);
	      RealType p2=StatsMap.computePriority(o, grayVal);
	      if (p1 > p2)
		{
		marker = o;
		}
	      // keep looking for collisions if we aren't marking boundaries
	      collision=false;
	      }
	    else
	      {  
	      break;
	      }
	    }
	  else
	    { 
	    marker = o; 
	    }
	  }
	}
      if( !collision )
	{
	// set the marker value
	outputIt.SetCenterPixel( marker );
	// update the region statistics
	StatsMap.updateRegion(marker, inputIt.GetCenterPixel());
	// and propagate to the neighbors
	for (niIt = inputIt.Begin(), nsIt = statusIt.Begin();
	     niIt != inputIt.End();
	     niIt++, nsIt++)
	  {
	  if ( !nsIt.Get() )
	    {
	    // the pixel is not yet processed. add it to the fah
	    InputImagePixelType grayVal = niIt.Get();
	    // compute priority using marker and grayVal
	    RealType priority=StatsMap.computePriority(marker, grayVal);
// 	    fah.Push(priority, 
// 		     inputIt.GetIndex() + niIt.GetNeighborhoodOffset()
// 		     ); 

	    // this won't work well with floating point input
	    priority = (RealType)(InputImagePixelType)round(priority);
	    fah[priority].push(this->GetMarkerImage()->ComputeOffset(inputIt.GetIndex() + niIt.GetNeighborhoodOffset()));
	    // mark it as already in the fah
	    nsIt.Set( true );
	    }
	  }
	}
      else
	{
	if( !m_MarkBoundaryLine )
	  {
	  // need to do something a bit dodgy - we checked which of the
	  // markers in the neigborhood is better, so assign that one
	  outputIt.SetCenterPixel( marker );
	  }
	}
      // one more pixel in the flooding stage
      progress.CompletedPixel();
      }
    }
}

template<class TInputImage,  class TLabelImage, class TRegionStats>
void 
SeededRegionGrowingBaseImageFilter<TInputImage, TLabelImage, TRegionStats>
::InitMap(StatsMapType &SM, LabelImageConstPointer LabelIm, InputImageConstPointer InputIm)
{
  // iterate over both images and build the map
  typedef typename itk::ImageRegionConstIterator<InputImageType> RawIterType;
  typedef typename itk::ImageRegionConstIterator<LabelImageType> LabIterType;

  RawIterType rit(InputIm, InputIm->GetLargestPossibleRegion());
  LabIterType lit(LabelIm, LabelIm->GetLargestPossibleRegion());

  rit.GoToBegin();
  lit.GoToBegin();
  while(!lit.IsAtEnd())
    {
    LabelImagePixelType V = lit.Get();
    if (V)
      {
      // found a region
      SM.initRegion(V, rit.Get());
      }
    ++lit;
    ++rit;
    }
}

template<class TInputImage,  class TLabelImage, class TRegionStats>
void
SeededRegionGrowingBaseImageFilter<TInputImage, TLabelImage, TRegionStats>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  
  os << indent << "FullyConnected: "  << m_FullyConnected << std::endl;
  os << indent << "MarkBoundaryLine: "  << m_MarkBoundaryLine << std::endl;
}
  
}// end namespace itk
#endif
