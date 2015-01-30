/*
 * SafGeometry.h
 *
 *  Created on: Dec 13, 2014
 *      Author: Daniel Saunders
 *        Desc: Container for all geometry information.
 */


#ifndef SAFGEOMETRY_H_
#define SAFGEOMETRY_H_
#include <vector>


class SafGeometry
{
private:
	// Members __________________________________________________________________
	unsigned int m_nGlibs;
	unsigned int m_nChannels;
	std::vector<bool> m_maskedChannels;


public:
	SafGeometry();
	virtual ~SafGeometry();


	// Setters and getters ______________________________________________________
	unsigned int nGlibs() {return m_nGlibs;}
	unsigned int nChannels() {return m_nChannels;}
	bool masked(unsigned int channelIndex) {return m_maskedChannels[channelIndex];}
};

#endif /* SAFGEOMETRY_H_ */
