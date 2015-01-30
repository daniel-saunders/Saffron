/*
 * SafTriggerDataSet.cpp
 *
 *  Created on: Dec 13, 2014
 *      Author: Daniel Saunders
 */

#include "SafTriggerDataSet.h"

//_____________________________________________________________________________

SafTriggerDataSet::SafTriggerDataSet(SafRunner * tempRunner) :
	SafDataSet(tempRunner, "SafTriggerDataSet")
{
	runner()->setTriggerData(this);
}


//_____________________________________________________________________________

SafTriggerDataSet::~SafTriggerDataSet()
{
	// TODO Auto-generated destructor stub
}


//_____________________________________________________________________________

void SafTriggerDataSet::clear()
{
  times()->clear();
	values()->clear();
	dipValues()->clear();
	peakValues()->clear();
	baseLines()->clear();
	channels()->clear();
}


//_____________________________________________________________________________
