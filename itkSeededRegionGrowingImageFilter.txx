#ifndef __itkSeededRegionGrowingImageFilter_txx
#define __itkSeededRegionGrowingImageFilter_txx

#include "itkSeededRegionGrowingImageFilter.h"

namespace itk {

template <class TInputImage, class TLabelImage>
SeededRegionGrowingImageFilter<TInputImage, TLabelImage>
::SeededRegionGrowingImageFilter()
{
 //  this->SetNumberOfRequiredInputs(2);
//   this->m_FullyConnected = false;
//   this->m_MarkBoundaryLine = true;
}


// template<class TInputImage, class TLabelImage>
// void
// SeededRegionGrowingImageFilter<TInputImage, TLabelImage>
// ::PrintSelf(std::ostream &os, Indent indent) const
// {
//   Superclass::PrintSelf(os, indent);
  
//   os << indent << "FullyConnected: "  << this->m_FullyConnected << std::endl;
//   os << indent << "MarkBoundaryLine: "  << this->m_MarkBoundaryLine << std::endl;
// }
  
}// end namespace itk
#endif
