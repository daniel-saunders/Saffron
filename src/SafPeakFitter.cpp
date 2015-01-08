/*
 * SafPeakFitter.cpp
 *
 *  Created on: Dec 22, 2014
 *      Author: Daniel Saunders
 */

#include "SafPeakFitter.h"


//_____________________________________________________________________________

SafPeakFitter::SafPeakFitter(SafRunner * runner) :
  SafAlgorithm(runner, "SafPeakFitter"),
  m_triggerWindowSizeA(16),
  m_triggerWindowSizeB(8),
  m_triggerWindowSizeC(8),
  m_triggerValueCut(100),
  m_nTriggers(0),
  m_caching(true)
{
	m_triggerWindowSizeTotal = m_triggerWindowSizeA + m_triggerWindowSizeB
			+ m_triggerWindowSizeC;
}


//_____________________________________________________________________________

SafPeakFitter::~SafPeakFitter()
{
}


//_____________________________________________________________________________

void SafPeakFitter::initialize()
{
	m_threading = true;

	h_triggerValues = initPerChannelPlots("FirstEventTriggerValues", "FirstEventTriggerValues", 
		runner()->eventTimeWindow(), 0.0, runner()->eventTimeWindow());

	TDirectory * instance_direc = runner()->saveFile()->mkdir(name().c_str());
	instance_direc->mkdir("FirstEventTriggerValues");
	for (unsigned int i=0; i<runner()->geometry()->nGlibs(); i++) {
		std::stringstream ssGlib; ssGlib<<i;
		instance_direc->mkdir(("FirstEventTriggerValues/Glib" + ssGlib.str()).c_str());
	}
}


//_____________________________________________________________________________

void SafPeakFitter::threadExecute(unsigned int iGlib, unsigned int iChannelLow, 
	unsigned int iChannelUp)
{
	for (unsigned int i=iChannelLow; i<std::min(iChannelUp, 
		runner()->geometry()->nChannels()); i++) {
		scanChannel(runner()->rawData()->channel(iGlib, i));
	}
}



//_____________________________________________________________________________

void SafPeakFitter::execute()
{
	unsigned int nGlibs = runner()->geometry()->nGlibs();
	unsigned int nChannels = runner()->geometry()->nChannels();
	
	for (unsigned int i=0; i<nGlibs; i++) {
		for (unsigned int j=0; j<nChannels; j++) {
			scanChannel(runner()->rawData()->channel(i, j));
		}
	}
	
	for (unsigned int i=0; i<nGlibs; i++) {
	  for (unsigned int j=0; j<nChannels; j++) {
			m_nTriggers += runner()->rawData()->channel(i, j)->nTriggers();
	  }
	}
}


//_____________________________________________________________________________

void SafPeakFitter::scanChannel(SafRawDataChannel * channel)
{
	unsigned int nTriggers = 0;
	std::vector<int> * times = channel->times();
	std::vector<double> * signals = channel->signals();
	if (times->size() <= m_triggerWindowSizeTotal) return;

	
	bool firstTimeEval = true;
	double triggerValue;
	double triggerDipValue;
	double triggerPeakValue;
	double triggerBaseLine;
	double cacheA = 0.0;
	double cacheB = 0.0; 
	double cacheC = 0.0;
	bool triggered = false;
  unsigned int plotIndex = channel->glibID() * runner()->geometry()->nChannels() +
		  channel->channelID();

	for (unsigned int i=0; i<times->size() - m_triggerWindowSizeTotal; i++) {
		// Temporary variables used to retrieve info from eval method. Also used 
		// for caching.		
		evalTimeWindow(signals, times, i, &triggerValue, &triggerDipValue, 
				&triggerPeakValue, &triggerBaseLine, &firstTimeEval, &cacheA, &cacheB,
				&cacheC);
		
	  if (event() == 0) h_triggerValues->at(plotIndex)->SetBinContent(i, triggerValue);
		
		if (triggerValue > m_triggerValueCut && !triggered) {
			channel->triggerTimes()->push_back(times->at(i));
			channel->triggerValues()->push_back(triggerValue);
			channel->triggerDipValues()->push_back(triggerDipValue);
			channel->triggerPeakValues()->push_back(triggerPeakValue);
			channel->triggerBaseLines()->push_back(triggerBaseLine);
			nTriggers++;
			triggered = true;
			//i += m_triggerWindowSizeB;
		}
		else if (triggerValue < m_triggerValueCut) triggered = false;

	}
	//std::exit(0);

	channel->setNTriggers(channel->triggerTimes()->size());
}


//_____________________________________________________________________________

void SafPeakFitter::evalTimeWindow(std::vector<double> * signals,
		std::vector<int> * times, unsigned int iStart, double * triggerValue, 
		double * triggerDipValue, double * triggerPeakValue, double * triggerBaseLine,
		bool * firstTimeEval, double * cacheA, double * cacheB, double * cacheC)
{
	// Cannot handle zero surpression in this version.
	// First eval case for this channel, and/or no caching cases.	
	
	if (*firstTimeEval || !m_caching) {
		(*cacheA) = 0.0;
		(*cacheB) = 0.0; 
		(*cacheC) = 0.0;

		unsigned int i;
		for (i = iStart; i<iStart + m_triggerWindowSizeA; i++)
			(*cacheA) += signals->at(i);

		for (; i<iStart + m_triggerWindowSizeA + m_triggerWindowSizeB; i++)
			(*cacheB) += signals->at(i);
	
		for (; i<iStart + m_triggerWindowSizeTotal; i++)
			(*cacheC) += signals->at(i);

		(*firstTimeEval) = false;
	}
	
	
	// Otherwise, fancy stuff with caching.
	else {
		unsigned int triggerWindowSizeAB = m_triggerWindowSizeA + m_triggerWindowSizeB;	
		(*cacheA) -= signals->at(iStart - 1);
		(*cacheA) += signals->at(iStart + m_triggerWindowSizeA - 1);
		
		(*cacheB) -= signals->at(iStart + m_triggerWindowSizeA - 1);
		(*cacheB) += signals->at(iStart + triggerWindowSizeAB - 1);
		
		(*cacheC) -= signals->at(iStart + triggerWindowSizeAB - 1);
		(*cacheC) += signals->at(iStart + m_triggerWindowSizeTotal - 1);
	}


	(*triggerBaseLine) = (*cacheA)/(1.0*m_triggerWindowSizeA);
	(*triggerDipValue) = (*cacheB)/(1.0*m_triggerWindowSizeB) - (*triggerBaseLine);
	(*triggerPeakValue) = (*cacheC)/(1.0*m_triggerWindowSizeC) - (*triggerBaseLine);
	(*triggerValue) = (*triggerPeakValue) - (*triggerDipValue);
}



//_____________________________________________________________________________

void SafPeakFitter::finalize()
{}


//_____________________________________________________________________________
