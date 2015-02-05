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
  m_nFileThreads(1),
  m_currentFileID(0),
  m_chainPos(0),
  m_triggerEventWindow(50000000) //16ns.
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
	if (runner()->runMode() == 1 || runner()->runMode() == 2) {
		m_chain = new TChain("waveforms");
		m_threading = false;
	}

	// MC Settings.
	m_uniformPeakRate = true;
	m_peakRate = 80;
	m_crossTorqueRatio = 0.55;
	m_singlePAHeight = 40;
	m_gain = 40;
	m_halfPeakWidth = 3;
	m_periodicNoiseA = 12;
	m_periodicNoisePeriod = 8;

	// Tree stuff.
	m_waveforms.push_back(new std::vector<unsigned short int>);
	m_glibs.push_back(0);
	m_glibchans.push_back(0);
	m_triggers.push_back(0);
	m_layers.push_back(0);
	m_chanxs.push_back(0);
	m_chanys.push_back(0);
	m_chanzs.push_back(0);
	m_triggerTime = 0;

	if (runner()->rawDataFileNames().size() == 0) {
		std::cout<<"No files given!"<<std::endl;
		exit(0);
	}

	m_fileNames = runner()->rawDataFileNames();
	if (runner()->runMode() == 1 || runner()->runMode() == 2) setupChain();

	else 	m_MCPAs = new TH1F("MCPAs", "MCPAs", 200, -5.5, 19.5);
	m_allSignals = new TH1F("allSignals", "allSignals", 4500, 7000, 16000);
	h_nSamplesPerEvent = new TH1F("nSamplesPerEvent", "nSamplesPerEvent",
			runner()->nEvents(), 0, runner()->nEvents());
	h_nSamplesPerEventPerChannel = new TH2F("nSamplesPerEventPerChannel", "nSamplesPerEventPerChannel",
			runner()->nCnG(), 0, runner()->nCnG(), runner()->nEvents(), 0, runner()->nEvents());

	runner()->saveFile()->cd();
}


//_____________________________________________________________________________

void SafEventBuilder::setupChain() {
	for (unsigned int i=0; i<m_fileNames.size(); i++) {
		std::cout<<"Reading file:\t"<<m_fileNames[i]<<std::endl;
		m_chain->Add(m_fileNames[i].c_str());
	}
	m_chain->Print();
	m_chain->SetBranchAddress("glib",&(m_glibs.back()));
	m_chain->SetBranchAddress("glibchan",&(m_glibchans.back()));
	if (runner()->runMode() == 1) {
		m_chain->SetBranchAddress("trigger",&(m_triggers.back()));
		m_chain->SetBranchAddress("layer",&(m_layers.back()));
	}
	else {
		m_chain->SetBranchAddress("triggertime", &m_triggerTime);
		m_chain->SetBranchAddress("chanz",&(m_chanzs.back()));
	}
	m_chain->SetBranchAddress("chanx",&(m_chanxs.back()));
	m_chain->SetBranchAddress("chany",&(m_chanys.back()));
	m_chain->SetBranchAddress("waveform",&(m_waveforms.back()));
	m_chain->GetEntry(0);
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
			if (runner()->runMode() == 1) scopeData(i);
			else triggerData();
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
	runner()->saveFile()->cd(name().c_str());
	if (runner()->runMode() == 0) {

		m_MCPAs->Write();
		m_allSignals->Write();
	}
	h_nSamplesPerEvent->Write();
	h_nSamplesPerEventPerChannel->Write();
}


//_____________________________________________________________________________

void SafEventBuilder::triggerData()
{
	// NOT SKIPPING SUPPORTED.
	while (m_triggerTime < m_event*m_triggerEventWindow) {
		m_chain->GetEntry(m_chainPos);
		m_chainPos++;
		SafRawDataChannel * channel = runner()->rawData()->channel(
				m_glibs.back()-1, m_glibchans.back());
		if (runner()->geometry()->masked(channel->plotIndex())) continue;

		unsigned int size = m_waveforms.back()->size();
		for (unsigned int i=0; i<size; i++) {
			channel->signals()->push_back(m_waveforms.back()->at(i));
			m_allSignals->Fill(m_waveforms.back()->at(i));
			//std::cout<<m_waveforms.back()->at(i)<<std::endl;
			channel->times()->push_back(m_triggerTime + i);
		}

		channel->addNEntries(size);
		channel->calcBaseLineEst();
		if (m_chainPos >= m_chain->GetEntries()) {
			m_eof = true;
			break;
		}
		h_nSamplesPerEvent->Fill(m_event, size);
		h_nSamplesPerEventPerChannel->Fill(channel->plotIndex(), m_event, size);
	}
}


//_____________________________________________________________________________

void SafEventBuilder::scopeData(unsigned int iThread)
{
	int skip = runner()->triggerSkip();


	// Skipping events (per file in chain) - NOT TESTED.
	while (m_triggers[iThread] < skip) {
		if (m_triggers[iThread] % 10 == 0 && m_glibchans[iThread] == 0 && m_glibs[iThread] == 111) {
			m_mtx.lock();
			std::cout<<"Trigger, skip trigger, thread: "<< m_triggers[iThread]<<"\t"<<skip<<"\t"<<iThread<<std::endl;
			m_mtx.unlock();
		}
		m_chain->GetEntry(m_chainPos);
		m_chainPos+=1;
	}


	// Reading.
	unsigned int trigger = m_triggers[iThread];
	while (m_triggers[iThread] == trigger) {
		m_chain->GetEntry(m_chainPos);
		m_chainPos++;
		SafRawDataChannel * channel = runner()->rawData()->channel(
				m_glibs[iThread]-111, m_glibchans[iThread]);

		unsigned int size = m_waveforms[iThread]->size();
		if (runner()->geometry()->masked(channel->plotIndex())) {
			for (unsigned int i=0; i<size; i++) {
				channel->signals()->push_back(0.0);
				channel->times()->push_back(i);
			}
		}
		else {
			for (unsigned int i=0; i<size; i++) {
				channel->signals()->push_back(m_waveforms[iThread]->at(i));
				m_allSignals->Fill(m_waveforms[iThread]->at(i));
				channel->times()->push_back(i);
			}
		}

		channel->addNEntries(channel->times()->size());
		if (m_event == 0) channel->calcBaseLineEst();
		if (m_chainPos >= m_chain->GetEntries()) {
			m_eof = true;
			break;
		}
		h_nSamplesPerEvent->Fill(m_event, size);
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
