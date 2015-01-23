/*
 * SafTrigger.cpp
 *
 *  Created on: Dec 22, 2014
 *      Author: Daniel Saunders
 */

#include "SafTrigger.h"

struct Local {
    Local(std::vector<int> * vec) {this->vec = vec;}
    bool operator () (int i, int j) {return (vec->at(i) < vec->at(j));}
    std::vector<int> * vec;
};

template <typename T>
std::vector<int> sort_permutation(std::vector<T> & vec)
{
    std::vector<int> p(vec.size());
    for (unsigned int i=0; i<vec.size(); i++) p[i] = i;
    std::sort(p.begin(), p.end(), Local(&vec));
    return p;
}

template <typename T>
std::vector<T> apply_permutation(
    std::vector<T> const& vec,
    std::vector<int> const& p)
{
    std::vector<T> sorted_vec(p.size());
    for (unsigned int i=0; i<p.size(); i++) {
    	sorted_vec[i] = vec[p[i]];
    }
    return sorted_vec;
}

//_____________________________________________________________________________

SafTrigger::SafTrigger(SafRunner * runner) :
  SafAlgorithm(runner, "SafTrigger"),
  m_triggerWindowSizeA(16),
  m_triggerWindowSizeB(8),
  m_triggerWindowSizeC(8),
  m_nTriggers(0),
  m_caching(true),
  m_triggerMethod(1),
  m_nSamplesWritten(0)
{
	m_triggerWindowSizeTotal = m_triggerWindowSizeA + m_triggerWindowSizeB
			+ m_triggerWindowSizeC;
}


//_____________________________________________________________________________

SafTrigger::~SafTrigger()
{
}


//_____________________________________________________________________________

void SafTrigger::initialize()
{
	if (!runner()->triggerData()) new SafTriggerDataSet(runner());
	m_threading = true;

	h_triggerValues = initPerChannelPlots("FirstEventTriggerValues", "FirstEventTriggerValues", 
		runner()->eventTimeWindow(), 0.0, runner()->eventTimeWindow());

	TDirectory * instance_direc = runner()->saveFile()->mkdir(name().c_str());
	instance_direc->mkdir("FirstEventTriggerValues");
	for (unsigned int i=0; i<runner()->geometry()->nGlibs(); i++) {
		std::stringstream ssGlib; ssGlib<<i;
		instance_direc->mkdir(("FirstEventTriggerValues/Glib" + ssGlib.str()).c_str());
	}
	m_triggerValueCuts.push_back(90);
	m_triggerValueCuts.push_back(runner()->triggerThreshold);
}


//_____________________________________________________________________________

void SafTrigger::threadExecute(unsigned int iGlib, unsigned int iChannelLow, 
	unsigned int iChannelUp, int iThread)
{
	for (unsigned int i=iChannelLow; i<std::min(iChannelUp, 
		runner()->geometry()->nChannels()); i++) {
		scanChannel(runner()->rawData()->channel(iGlib, i));
	}
}


//_____________________________________________________________________________

void SafTrigger::execute()
{
	for (unsigned int i=0; i<runner()->geometry()->nGlibs(); i++){
		threadExecute(i, 0, runner()->geometry()->nChannels(), -1);
	}
}


//_____________________________________________________________________________

void SafTrigger::postExecute()
{
	std::vector<int> * triggerTimes = runner()->triggerData()->times();
	std::vector<SafRawDataChannel*> * channels = runner()->triggerData()->channels();
	std::vector<double> * values = runner()->triggerData()->values();
	std::vector<double> * dipValues = runner()->triggerData()->dipValues();
	std::vector<double> * peakValues = runner()->triggerData()->peakValues();
	std::vector<double> * baseLines = runner()->triggerData()->baseLines();

	auto p = sort_permutation(*triggerTimes);

	(*triggerTimes) = apply_permutation((*triggerTimes), p);
	(*channels) = apply_permutation((*channels), p);
	(*values) = apply_permutation((*values), p);
	(*dipValues) = apply_permutation((*dipValues), p);
	(*peakValues) = apply_permutation((*peakValues), p);
	(*baseLines) = apply_permutation((*baseLines), p);
}


//_____________________________________________________________________________

void SafTrigger::scanChannel(SafRawDataChannel * channel)
{
	unsigned int nTriggers = 0;
	unsigned int nTriggeredSamplesWritten = 0;
	unsigned int nSinceLastTrigger = 10;
	bool writing = false;
	unsigned int nWritten = 0;
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
		
	  double tempTriggerValue;
	  if (m_triggerMethod == 0) tempTriggerValue = triggerValue;
	  else if (m_triggerMethod == 1)
	  	tempTriggerValue = signals->at(i) - channel->baseLineEst();
		if (tempTriggerValue > m_triggerValueCuts[m_triggerMethod] && !triggered && channel->baseLineEst() > 5000) {
			std::vector<double>::iterator iSigMax = std::max_element(
					signals->begin() + i, signals->begin()+i+m_triggerWindowSizeTotal);
			double val = (*iSigMax) + (*(iSigMax-1)) + (*(iSigMax+1)) + (*(iSigMax-2)) + (*(iSigMax+2)) - 5*triggerBaseLine;
			double time = times->at(i);

			m_mtx.lock();
			runner()->triggerData()->times()->push_back(time);
			runner()->triggerData()->channels()->push_back(channel);
			runner()->triggerData()->values()->push_back(val);
			runner()->triggerData()->dipValues()->push_back(triggerDipValue);
			runner()->triggerData()->peakValues()->push_back(triggerPeakValue);
			runner()->triggerData()->baseLines()->push_back(triggerBaseLine);
			m_nTriggers++;
			m_mtx.unlock();


			nTriggers++;
			triggered = true;
			nSinceLastTrigger = 0;
			unsigned int prevWindow = 8;

			// Data counting _______
			if (!writing) {
				writing = true;
				if (nSinceLastTrigger < prevWindow) nTriggeredSamplesWritten+=nSinceLastTrigger;
				else nTriggeredSamplesWritten+=prevWindow;
			}

			nWritten = prevWindow;

		}
		else if (triggerValue < m_triggerValueCuts[m_triggerMethod]) triggered = false;


		// Data counting _______
		if (writing) {
			nWritten += 1;
			nTriggeredSamplesWritten += 1;
			unsigned int stop = 32;
			if (nWritten > stop) {
				nWritten = 0;
				writing = false;
			}
		}
		else nSinceLastTrigger++;
	}

	m_mtx.lock();
	m_nSamplesWritten += nTriggeredSamplesWritten;
	m_mtx.unlock();

	channel->addNTriggers(nTriggers);
	channel->addNTriggerSamplesWritten(nTriggeredSamplesWritten);
}


//_____________________________________________________________________________

void SafTrigger::evalTimeWindow(std::vector<double> * signals,
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

void SafTrigger::finalize()
{
	std::cout<<"nTriggers:\t"<<m_nTriggers<<"\t"<<"TriggerThreshold:\t"<<m_triggerValueCuts[m_triggerMethod]<<std::endl;
	std::cout<<"TriggerRate(kHz)\t"<<m_nTriggers/(1000.*runner()->realTimeElapsed())<<std::endl;
	std::cout<<"DataRate(MBs)\t"<<(m_nSamplesWritten*2)/(1000000.*runner()->realTimeElapsed());

	for (unsigned int i=0; i<h_triggerValues->size(); i++) {
		int iGlib = i/runner()->geometry()->nChannels();
		std::stringstream ssGlib; ssGlib << iGlib;
		runner()->saveFile()->cd((name() + "/FirstEventTriggerValues/Glib" + ssGlib.str()).c_str());
		h_triggerValues->at(i)->Write();
	}
}


//_____________________________________________________________________________

