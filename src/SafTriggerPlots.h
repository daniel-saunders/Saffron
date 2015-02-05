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
	TH1F * h_valuesZoomed;
	TH2F * h_integralsVsValues;
	TH1F * h_nTriggers;
	TH2F * h_valuesVsChannel;
	TH1F * h_dataRates;
	TH1F * h_triggerRate;
	TH2F * h_nTriggersVsEvents;
	TH1F * h_integrals;
	TH2F * h_integralsVsChannel;
	TH2F * h_integralsVsValues0to38;
	TH2F * h_integralsVsValues38to76;
	TH2F * h_integralsVsValuesOdd;
	TH2F * h_integralsVsValuesEven;

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
