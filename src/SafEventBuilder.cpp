/*
 * SafEventBuilder.cpp
 *
 *  Created on: Dec 13, 2014
 *      Author: ds9133
 */

#include "SafEventBuilder.h"

//_____________________________________________________________________________

SafEventBuilder::SafEventBuilder(SafRunner * runner) :
  SafAlgorithm(runner, "SafEventBuilder"),
  m_nFileThreads(1)
{
  m_mean = 8000;
  m_rms = 7;
  m_firstTime = true;
}


//_____________________________________________________________________________

SafEventBuilder::~SafEventBuilder()
{
	// TODO Auto-generated destructor stub
}


//_____________________________________________________________________________

void SafEventBuilder::initialize()
{
	m_threading = true;

	if (runner()->runMode() == 1) m_threading = false;
	m_uniformPeakRate = true;
	m_peakRate = 80;
	m_crossTorqueRatio = 0.55;
	m_singlePAHeight = 30;
	m_gain = 30;
	m_halfPeakWidth = 3;
	m_periodicNoiseA = 12;
	m_periodicNoisePeriod = 8;


	m_fileNames.push_back(runner()->fileName());
	if (runner()->runMode() == 1) {
		for (unsigned int i=0; i<m_nFileThreads; i++) {
			std::cout<<"Reading file:\t"<<m_fileNames[i]<<std::endl;
			m_files.push_back(new TFile(m_fileNames[i].c_str(), "READ"));
			m_trees.push_back((TTree*)m_files.back()->Get("waveforms"));


			m_waveforms.push_back(new std::vector<int>);
			m_glibs.push_back(0);
			m_glibchans.push_back(0);
			m_triggers.push_back(0);
			m_layers.push_back(0);
			m_chanxs.push_back(0);
			m_chanys.push_back(0);
			m_treePos.push_back(i);

			m_trees.back()->SetBranchAddress("glib",&(m_glibs.back()));
			m_trees.back()->SetBranchAddress("glibchan",&(m_glibchans.back()));
			m_trees.back()->SetBranchAddress("trigger",&(m_triggers.back()));
			m_trees.back()->SetBranchAddress("layer",&(m_layers.back()));
			m_trees.back()->SetBranchAddress("chanx",&(m_chanxs.back()));
			m_trees.back()->SetBranchAddress("chany",&(m_chanys.back()));
			m_trees.back()->SetBranchAddress("waveform",&(m_waveforms.back()));
			m_trees.back()->GetEvent(0);
		}
		m_spareWaveform = new std::vector<int>;
	}

	else {
		m_MCPAs = new TH1F("MCPAs", "MCPAs", 200, -5.5, 19.5);
	}
	m_allSignals = new TH1F("allSignals", "allSignals", 4500, 7000, 16000);

	runner()->saveFile()->cd();
}


//_____________________________________________________________________________

void SafEventBuilder::execute()
{
	if (runner()->runMode() == 0) {
		for (unsigned int i=0; i<runner()->geometry()->nGlibs(); i++)
			monteCarlo(i, 0, runner()->geometry()->nChannels(), -1);
	}

	else {
		for (unsigned int i=0; i<m_nFileThreads; i++) {
			realData(i);
		}
	}
}


//_____________________________________________________________________________

void SafEventBuilder::threadExecute(unsigned int iGlib, unsigned int iLow, 
	unsigned int iUp, int iThread) {
	monteCarlo(iGlib, iLow, iUp, iThread);
}


//_____________________________________________________________________________

void SafEventBuilder::finalize()
{
	if (runner()->runMode() == 0) {
		runner()->saveFile()->cd();
		m_MCPAs->Write();
		m_allSignals->Write();
	}
}


//_____________________________________________________________________________

void SafEventBuilder::realData(unsigned int iThread)
{
	int skip = runner()->triggerSkip();
	if (m_firstTime) {
		for (unsigned int i=0; i<m_nFileThreads; i++) {
			m_trees[i]->GetEntry(m_treePos[i]);
		}
		m_firstTime = false;
	}

	// Read tree.
	while (m_triggers[iThread] < skip) {
		if (m_triggers[iThread] % 10 == 0 && m_glibchans[iThread] == 0 && m_glibs[iThread] == 111) {
			m_mtx.lock();
			std::cout<<"Trigger, skip trigger, thread: "<< m_triggers[iThread]<<"\t"<<skip<<"\t"<<iThread<<std::endl;
			m_mtx.unlock();
		}
		m_trees[iThread]->GetEntry(m_treePos[iThread]);
		m_treePos[iThread]+=m_nFileThreads;
	}


	int limit = (int)runner()->event() + skip;
	while (m_triggers[iThread] <= limit) {
		m_trees[iThread]->GetEntry(m_treePos[iThread]);
		m_treePos[iThread]+=m_nFileThreads;
		SafRawDataChannel * channel = runner()->rawData()->channel(
				m_glibs[iThread]-111, m_glibchans[iThread]);

		unsigned int size = m_waveforms[iThread]->size();
		for (unsigned int i=0; i<size; i++) {
			channel->signals()->push_back(m_waveforms[iThread]->at(i));
			m_allSignals->Fill(m_waveforms[iThread]->at(i));
			channel->times()->push_back(i);
		}

		channel->addNEntries(channel->times()->size());
		if (m_event == 0) channel->calcBaseLineEst();
		if (m_treePos[iThread] >= m_trees[iThread]->GetEntries()) {
			m_eof = true;
			break;
		}
	}
}



//_____________________________________________________________________________

void SafEventBuilder::monteCarlo(unsigned int iGlib, unsigned int iLow, 
	unsigned int iUp, int iThread)
{
	TRandom3 * randGen = new TRandom3();
	randGen->SetSeed();
	for (unsigned int ic = iLow; ic < iUp; ic++) {
		SafRawDataChannel * channel = runner()->rawData()->channel(iGlib, ic);
		addPedestal(channel, randGen);
		if (m_uniformPeakRate) addUniformPeaks(channel, randGen);
		channel->addNEntries(runner()->eventTimeWindow());
		if (m_event == 0) channel->calcBaseLineEst();
	}
	delete randGen;
}


//_____________________________________________________________________________

void SafEventBuilder::addPedestal(SafRawDataChannel * channel, TRandom3 * randGen) {
	for (unsigned int i=0; i<runner()->eventTimeWindow(); i++) {
		double signal = randGen->Gaus(m_mean, m_rms);
		signal += m_periodicNoiseA*sin(i*3.14/m_periodicNoisePeriod);
		channel->times()->push_back(i);
		channel->signals()->push_back(signal);
		// m_mtx.lock();
		// m_allSignals->Fill(signal);
		// m_mtx.unlock();
	}
}


//_____________________________________________________________________________

void SafEventBuilder::addUniformPeaks(SafRawDataChannel * channel, TRandom3 * randGen) {
	// Get peak spacing, avoid edge effects, so plus 2.
	int spacing = runner()->eventTimeWindow()/(m_peakRate+2);
	for (unsigned int i=0; i<m_peakRate; i++) {
		int peakTime = spacing*(i+1);
		int nPA = (int) randGen->Exp(m_crossTorqueRatio);
		//m_mtx.lock();
		m_MCPAs->Fill(nPA);
		//m_mtx.unlock();
		double amplitude = m_singlePAHeight + nPA*m_gain;
		double grad = amplitude/(1.*m_halfPeakWidth);
		for (unsigned int j=1; j<=m_halfPeakWidth; j++) {
			if (peakTime + j + m_halfPeakWidth > runner()->eventTimeWindow()) break;
			channel->signals()->at(peakTime + j) += j*grad;
			if (j!=m_halfPeakWidth) channel->signals()->at(peakTime + 2*m_halfPeakWidth - j) += j*grad;
		}
	}
}


//_____________________________________________________________________________
