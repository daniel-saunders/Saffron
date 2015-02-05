/*
 * SafGeometry.cpp
 *
 *  Created on: Dec 13, 2014
 *      Author: ds9133
 */

#include "SafGeometry.h"

//_____________________________________________________________________________

SafGeometry::SafGeometry() :
  m_nGlibs(5),
  m_nChannels(76),
  m_maskedChannels(500, false)
{
//	m_maskedChannels[164] = true;
//	m_maskedChannels[155] = true;
//	m_maskedChannels[50] = true;
//	m_maskedChannels[11] = true;
//	m_maskedChannels[12] = true;
	m_maskedChannels[31] = true;
	m_maskedChannels[10] = true;
	m_maskedChannels[89] = true;
}


//_____________________________________________________________________________

SafGeometry::~SafGeometry()
{
	// TODO Auto-generated destructor stub
}


//_____________________________________________________________________________
