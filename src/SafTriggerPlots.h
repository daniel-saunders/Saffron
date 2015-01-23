/*
 * SafTriggerPlots.h
 *
 *  Created on: Dec 30, 2014
 *      Author: Daniel Saunders
 */

#ifndef SAFTRIGGERPLOTS_H_
#define SAFTRIGGERPLOTS_H_

#include "SafAlgorithm.h"
#include "SafRunner.h"
#include "SafRawDataChannel.h"

class SafTriggerPlots: public SafAlgorithm
{
private:
	// Members __________________________________________________________________

	std::vector<TH1F*> * h_firstEventPeaks;
	std::vector<TH1F*> * h_valuesPerChannel;
	TH1F * h_dipValues;
	TH1F * h_peakValues;
	TH1F * h_values;
	TH2F * h_dipVsPeakValues;
	TH2F * h_dipVsValues;
	TH2F * h_peakVsValues;
	TH1F * h_nTriggers;
	TH2F * h_valuesVsChannel;
	TH1F * h_dataRates;
	TH1F * h_triggerRate;
	TH2F * h_nTriggersVsEvents;

public:
  // Methods __________________________________________________________________
	SafTriggerPlots(SafRunner * runner);
	virtual ~SafTriggerPlots();

	void initialize();
	void execute();
	void finalize();

	void fill();
};

#endif /* SAFTRIGGERPLOTS_H_ */
