/*
 * SafTrigger.cpp
 *
 *  Created on: Dec 22, 2014
 *      Author: Daniel Saunders
 */

#include "SafTrigger.h"

struct Local {
    Local(std::vector<long long int> * vec) {this->vec = vec;}
    bool operator () (int i, int j) {return (vec->at(i) < vec->at(j));}
    std::vector<long long int> * vec;
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

	double size;
	if (runner()->runMode() == 1) size = runner()->eventTimeWindow();
	else size = 300;

	h_triggerValues = initPerChannelPlots("FirstEventTriggerValues", "FirstEventTriggerValues", 
		size, 0.0, size);

	h_triggerOnOff = initPerChannelPlots("FirstEventTriggerOnOff", "FirstEventTriggerOnOff",
		size, 0.0, size);

	h_writingOnOff = initPerChannelPlots("FirstEventWritingOnOff", "FirstEventWritingOnOff",
		size, 0.0, size);

	TDirectory * instance_direc = runner()->saveFile()->mkdir(name().c_str());
	instance_direc->mkdir("FirstEventTriggerValues");
	instance_direc->mkdir("FirstEventTriggerOnOff");
	instance_direc->mkdir("FirstEventWritingOnOff");
	for (unsigned int i=0; i<runner()->geometry()->nGlibs(); i++) {
		std::stringstream ssGlib; ssGlib<<i;
		instance_direc->mkdir(("FirstEventTriggerValues/Glib" + ssGlib.str()).c_str());
		instance_direc->mkdir(("FirstEventTriggerOnOff/Glib" + ssGlib.str()).c_str());
		instance_direc->mkdir(("FirstEventWritingOnOff/Glib" + ssGlib.str()).c_str());
	}
	m_triggerValueCuts.push_back(90);
	m_triggerValueCuts.push_back(runner()->triggerThreshold);

	std::cout<<"Trigger threshold: "<<m_triggerValueCuts[m_triggerMethod]<<std::endl;
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
	std::vector<long long int> * triggerTimes = runner()->triggerData()->times();
	std::vector<SafRawDataChannel*> * channels = runner()->triggerData()->channels();
	std::vector<double> * values = runner()->triggerData()->values();
	std::vector<double> * integrals = runner()->triggerData()->integrals();

	auto p = sort_permutation(*triggerTimes);

	(*triggerTimes) = apply_permutation((*triggerTimes), p);
	(*channels) = apply_permutation((*channels), p);
	(*values) = apply_permutation((*values), p);
	(*integrals) = apply_permutation((*integrals), p);

	if (m_event == 0) {
		TH1F * h_triggerSort = new TH1F("triggerSorting", "triggerSorting", triggerTimes->size(), 0, triggerTimes->size());
		for (unsigned int i=0; i<triggerTimes->size(); i++) {
			h_triggerSort->SetBinContent(i, triggerTimes->at(i));
		}
		h_triggerSort->Write();
		delete h_triggerSort;
	}

	for (unsigned int i=0; i<channels->size(); i++)
		channels->at(i)->triggerIDs()->push_back(i);
}


//_____________________________________________________________________________

void SafTrigger::scanChannel(SafRawDataChannel * channel)
{
	unsigned int nTriggers = 0;
	unsigned int nTriggeredSamplesWritten = 0;
	unsigned int nSinceLastTrigger = 10;
	bool writing = false;
	unsigned int nWritten = 0;
	std::vector<long long int> * times = channel->times();
	std::vector<double> * signals = channel->signals();
	if (times->size() <= m_triggerWindowSizeTotal && m_triggerMethod == 0) return;

	
	bool firstTimeEval = true;
	double triggerValue;
	double cacheA = 0.0;
	double cacheB = 0.0; 
	double cacheC = 0.0;
	double triggerDipValue, triggerPeakValue, triggerBaseLine;
	bool triggered = false;
  unsigned int plotIndex = channel->glibID() * runner()->geometry()->nChannels() +
		  channel->channelID();

	for (unsigned int i=0; i<times->size(); i++) {
		if (m_triggerMethod == 0 && i > times->size() - m_triggerWindowSizeTotal) break;
		if (m_triggerMethod == 0) evalTimeWindow(signals, times, i, &triggerValue, &triggerDipValue,
				&triggerPeakValue, &triggerBaseLine, &firstTimeEval, &cacheA, &cacheB,
				&cacheC);
		
	  double tempTriggerValue;
	  if (m_triggerMethod == 0) tempTriggerValue = triggerValue;
	  else tempTriggerValue = signals->at(i) - channel->baseLineEst();
		if (tempTriggerValue > m_triggerValueCuts[m_triggerMethod] && !triggered && channel->baseLineEst() > 5000) {
			long long int time = times->at(i);

			double integral = 0.0;
			if (i+5<signals->size() && i>1) {
				std::vector<double>::iterator peakPos = std::max_element(signals->begin() + i, signals->begin() + i + 5);
				integral += *peakPos;
				integral += *(peakPos-1);
				integral += *(peakPos+1);
				integral -= 3*channel->baseLineEst();
			}

			m_mtx.lock();
			runner()->triggerData()->times()->push_back(time);
			runner()->triggerData()->channels()->push_back(channel);
			runner()->triggerData()->values()->push_back(tempTriggerValue);
			runner()->triggerData()->integrals()->push_back(integral);
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

//		if (m_event == 0) {
//			h_triggerValues stuff
//			h_triggerOnOff->at(channel->plotIndex())->SetBinContent(i, int(triggered));
//			h_writingOnOff->at(channel->plotIndex())->SetBinContent(i, int(writing));
//		}
	}

	m_mtx.lock();
	m_nSamplesWritten += nTriggeredSamplesWritten;
	m_mtx.unlock();

	channel->addNTriggers(nTriggers);
	channel->addNTriggerSamplesWritten(nTriggeredSamplesWritten);
}


//_____________________________________________________________________________

void SafTrigger::evalTimeWindow(std::vector<double> * signals,
		std::vector<long long int> * times, unsigned int iStart, double * triggerValue,
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

	for (unsigned int i=0; i<h_triggerValues->size(); i++) {
		int iGlib = i/runner()->geometry()->nChannels();
		std::stringstream ssGlib; ssGlib << iGlib;
		runner()->saveFile()->cd((name() + "/FirstEventTriggerOnOff/Glib" + ssGlib.str()).c_str());
		h_triggerOnOff->at(i)->Write();
	}

	for (unsigned int i=0; i<h_triggerValues->size(); i++) {
		int iGlib = i/runner()->geometry()->nChannels();
		std::stringstream ssGlib; ssGlib << iGlib;
		runner()->saveFile()->cd((name() + "/FirstEventWritingOnOff/Glib" + ssGlib.str()).c_str());
		h_writingOnOff->at(i)->Write();
	}
}


//_____________________________________________________________________________

