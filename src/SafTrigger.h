/*
 * SafTrigger.h
 *
 *  Created on: Dec 22, 2014
 *      Author: Daniel Saunders
 */

#ifndef SAFTRIGGER_H_
#define SAFTRIGGER_H_

#include "SafAlgorithm.h"
#include "SafRunner.h"
#include "SafRawDataChannel.h"
#include <thread>

class SafRawDataChannel;

class SafTrigger: public SafAlgorithm
{
private:
	// Members __________________________________________________________________
  unsigned int m_triggerWindowSizeTotal;
  unsigned int m_triggerWindowSizeA;
  unsigned int m_triggerWindowSizeB;
  unsigned int m_triggerWindowSizeC;
  double m_triggerValueCut;
  unsigned int m_nTriggers;
  bool m_caching;
  std::vector<TH1F*> * h_triggerValues;


public:
  // Methods __________________________________________________________________
	SafTrigger(SafRunner * runner);
	virtual ~SafTrigger();

	void initialize();
	void execute();
	void finalize();
  void threadExecute(unsigned int iGlib, unsigned int iLow, unsigned int iUp);

	void scanChannel(SafRawDataChannel * channel);
	void scanGlib(unsigned int iGlib);
	void evalTimeWindow(std::vector<double> * signals,
		std::vector<int> * times, unsigned int i, double * triggerValue, 
		double * triggerDipValue, double * triggerPeakValue, double * triggerBaseLine, 
		bool * firstTimeEval, double * cacheA, double * cacheB, double * cacheC);
};

#endif /* SAFTRIGGER_H_ */
