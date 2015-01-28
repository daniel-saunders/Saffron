/*
 * SafCoincidenceFinder.cpp
 *
 *  Created on: Jan 13, 2014
 *      Author: Daniel Saunders
 */

#include "SafCoincidenceFinder.h"


//_____________________________________________________________________________

SafCoincidenceFinder::SafCoincidenceFinder(SafRunner * runner) :
  SafAlgorithm(runner, "SafCoincidenceFinder"),
  m_timeWindow(5)
{
}


//_____________________________________________________________________________

SafCoincidenceFinder::~SafCoincidenceFinder()
{
}


//_____________________________________________________________________________

void SafCoincidenceFinder::initialize()
{
	m_threading = false;
	h_size = new TH1F("size", "size", 40, -0.5, 39.5);
	h_duration = new TH1F("duration", "duration", m_timeWindow+4, -0.5, m_timeWindow+3.5);
	h_sizeVsDuration = new TH2F("sizeVsDuration", "sizeVsDuration", 40, -0.5, 39.5,
			m_timeWindow+4, -0.5, m_timeWindow+3.5);
	h_channelRate = new TH1F("channelRate", "channelRate", runner()->nCnG(), -0.5, runner()->nCnG()-0.5);
	TDirectory * instance_direc = runner()->saveFile()->mkdir(name().c_str());
}


//_____________________________________________________________________________

void SafCoincidenceFinder::threadExecute(unsigned int iGlib, unsigned int iChannelLow,
	unsigned int iChannelUp)
{
	std::vector<int> * triggerTimes = runner()->triggerData()->times();
	std::vector<SafRawDataChannel*> * channels = runner()->triggerData()->channels();
	std::vector<SafRawDataChannel*> coinChannels;

	for (unsigned int iTime = 0; iTime < triggerTimes->size()-1; iTime++) {
		unsigned int size = 1;
		int duration = 0;
		int coinDuration = 0;

		coinChannels.push_back(channels->at(iTime));
		for (unsigned int jTime = iTime + 1; jTime < triggerTimes->size(); jTime++) {
			duration = triggerTimes->at(jTime) - triggerTimes->at(iTime);
			if (channels->at(iTime) == channels->at(jTime)) continue;
			if (duration < m_timeWindow)	{
				coinDuration = duration;
				size++;
				coinChannels.push_back(channels->at(jTime));
			}
			else break;
		}

		if (size > 1) {
			h_size->Fill(size);
			h_duration->Fill(coinDuration);
			h_sizeVsDuration->Fill(size, coinDuration);

			for (unsigned int i=0; i<size; i++) {
				h_channelRate->Fill(coinChannels[i]->plotIndex());
			}
		}
		coinChannels.clear();
	}
}



//_____________________________________________________________________________

void SafCoincidenceFinder::execute()
{
	for (unsigned int i=0; i<runner()->geometry()->nGlibs(); i++){
		threadExecute(i, 0, runner()->geometry()->nChannels());
	}
}



//_____________________________________________________________________________

void SafCoincidenceFinder::finalize()
{
	runner()->saveFile()->cd(name().c_str());
	h_size->Write();
	h_sizeVsDuration->Write();
	h_duration->Write();
	h_channelRate->Write();
}


//_____________________________________________________________________________

