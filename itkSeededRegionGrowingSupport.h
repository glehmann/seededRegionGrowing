#ifndef __itk_SeededRegionGrowingSupport_h
#define __itk_SeededRegionGrowingSupport_h

#include <map>
#include "itkNumericTraits.h"

namespace itk {

template <class LabType, class PixType>
class SRGRegionStatsType
{
public:
  typedef typename NumericTraits<PixType>::RealType RealType;

  RealType computePriority(LabType region,
			   PixType  candidate)
  {
    typename StatsMapType::const_iterator reg = m_StatsMap.find(region);
    RealType regmean = reg->second.m_Sum/((RealType)reg->second.m_Count);
//     if (region == 1) 
//       {

//       std::cout << "Reg 1" << regmean << " " << (int)candidate<< std::endl;
//       }
    return (RealType)fabs(regmean - (RealType)candidate);
  }

  void updateRegion(LabType region,
		    PixType candidate)
  {
    typename StatsMapType::iterator reg = m_StatsMap.find(region);
    ++(reg->second.m_Count);
    reg->second.m_Sum += (RealType)candidate;
  }

  void initRegion(LabType region,
		  PixType candidate)
  {
    ++(m_StatsMap[region].m_Count);
    m_StatsMap[region].m_Sum += (RealType)candidate;
   // std::cout << "Added " << (int) region << " " << m_StatsMap[region].m_Count << " " << m_StatsMap[region].m_Sum << std::endl;
  }

private:
  typedef class RegionStatisticsType
  {
  public:
    RegionStatisticsType()
    {
      m_Count = 0;
      m_Sum = NumericTraits<RealType>::Zero;
    }
    
    RealType m_Sum;
    unsigned long m_Count;
  } RegionStatisticsType;

  typedef typename std::map<LabType, RegionStatisticsType> StatsMapType;

  StatsMapType m_StatsMap;
};

}
#endif
