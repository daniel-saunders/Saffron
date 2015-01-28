/*
 * SafEventBuilder.h
 *
 *  Created on: Dec 13, 2014
 *      Author: Daniel Saunders
 */

#ifndef SAFEVENTBUILDER_H_
#define SAFEVENTBUILDER_H_
#include "SafAlgorithm.h"
#include "SafRunner.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TFile.h"


class SafRawDataChannel;


class SafEventBuilder: public SafAlgorithm
{
private:
	// Members __________________________________________________________________
	double m_mean; //MC mode.
	double m_rms; //MC mode.
	std::vector<unsigned int> m_treePos;
	std::vector<TTree*> m_trees;
	std::vector<TFile*> m_files;
	std::vector<std::string> m_fileNames; // one per thread.
	std::vector<int> m_glibs;
	std::vector<int> m_glibchans;
	std::vector<int> m_triggers;
	std::vector<int> m_layers;
	std::vector<int> m_chanxs;
	std::vector<int> m_chanys;
	std::vector< std::vector<int>* > m_waveforms;
	std::vector<int> * m_spareWaveform;
	bool m_firstTime;
	unsigned int m_nFileThreads; // Doesn't work, should == 1 always.
	TH1F * m_allSignals;


	// MC stuff.
	double m_gain;
	double m_crossTorqueRatio;
	unsigned int m_peakRate;
	bool m_uniformPeakRate; //n peaks per event per channel.
	unsigned int m_halfPeakWidth;
	double m_singlePAHeight;
	TH1F * m_MCPAs;
	double m_periodicNoiseA;
	double m_periodicNoisePeriod;


public:
	// Methods __________________________________________________________________
	SafEventBuilder(SafRunner * runner);
	virtual ~SafEventBuilder();

	void initialize();
	void execute();
	void threadExecute(unsigned int iGlib, unsigned int iLow, unsigned int iUp, int);
	void finalize();

	void monteCarlo(unsigned int iGlib, unsigned int iLow, unsigned int iUp, int iThread);
	void addPedestal(SafRawDataChannel * channel, TRandom3 * randGen);
	void addUniformPeaks(SafRawDataChannel * channel, TRandom3 * randGen);
	void threadFill(SafRawDataChannel * channel, std::vector<int> * waveform);
	void realData(unsigned int channelIndexUpper);
	TTree * tree() {return m_trees.back();}
	unsigned int treePos() {return m_treePos.back();}
};

#endif /* SAFEVENTBUILDER_H_ */
