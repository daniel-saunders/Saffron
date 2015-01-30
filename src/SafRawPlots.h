/*
 * SafRawPlots.h
 *
 *  Created on: Dec 15, 2014
 *      Author: ds9133
 */

#ifndef SAFRAWPLOTS_H_
#define SAFRAWPLOTS_H_

#include "SafAlgorithm.h"
#include "SafRunner.h"
#include <sstream>
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TDirectory.h"
#include "TMath.h"


class SafRawPlots: public SafAlgorithm
{
private:
	// Members __________________________________________________________________
	std::vector<TH1F*> * h_firstEventWaveforms;
	std::vector<TH1F*> * h_signals;
	std::vector<TH1F*> * h_signalsDiff;
	std::vector<TH1F*> * h_signalsDoubleDiff;
	std::vector<TH1F*> * h_signalsTripleDiff;
	TH2F * h_allSignals;
	TH1F * h_signalMeans;
	TH1F * h_signalWidths;
	bool m_filtered;
	TH2F * h_avSignalPerEventPerChannel;
	TH2F * h_rmsSignalPerEventPerChannel;
	TH1F * h_avSignalOfEvents;
	TH2F * h_avSignalPerChannelOfEvents;
	unsigned int m_diffBinRange;
	unsigned int m_nSeekedRoots;
	TH1F * h_gainsPerChannel;
	TH1F * h_gains;
	bool m_calculateGains;
	unsigned int m_nFinalizeThreads;
	TH1F * h_nBaseLineEstVsChannel;
	bool m_smoothing;


public:
	// Methods __________________________________________________________________
  SafRawPlots(SafRunner * runner, bool filtered);
	virtual ~SafRawPlots();

	void initialize();
	void execute();
	void finalize();
  void threadExecute(unsigned int iGlib, unsigned int iLow, unsigned int iUp,
  	int iThread);

  void calculateGains(unsigned int iLow, unsigned iUp);
  void smooth(TH1F * h, int plotIndex);
};

#endif /* SAFRAWPLOTS_H_ */
