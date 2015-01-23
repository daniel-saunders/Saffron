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
	m_randGen = new TRandom3();
  m_randGen->SetSeed();
  m_mean = 8000;
  m_rms = 100;
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
	m_fileNames.push_back(runner()->fileName());


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
	runner()->saveFile()->cd();
}


//_____________________________________________________________________________

void SafEventBuilder::execute()
{
	if (runner()->runMode() == 0) {
		monteCarlo();
	}
	else {
		for (unsigned int i=0; i<m_nFileThreads; i++) {
			realData(i);
		}
	}
}

//_____________________________________________________________________________

void SafEventBuilder::finalize()
{
	//m_files[0]->Close();
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
			channel->times()->push_back(i);
		}
		channel->setNEntries(channel->times()->size());
		if (m_treePos[iThread] >= m_trees[iThread]->GetEntries()) m_eof = true;
	}
}



//_____________________________________________________________________________

void SafEventBuilder::monteCarlo()
{
	for (unsigned int ig = 0; ig < runner()->geometry()->nGlibs(); ig++) {
		for (unsigned int ic = 0; ic < runner()->geometry()->nChannels(); ic++) {
			SafRawDataChannel * channel = runner()->rawData()->channel(ig, ic);
			for (unsigned int i=0; i<runner()->eventTimeWindow(); i++) {
				int time = i;
				double signal = m_randGen->Gaus(m_mean, m_rms);
				channel->times()->push_back(time);
				channel->signals()->push_back(signal);
				plot(signal, "SafEventBuilder/allSignals", "allSignals", 100, 7000, 9000);
			}
			channel->addNEntries(runner()->eventTimeWindow());
		}
	}
}


//_____________________________________________________________________________

