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
  m_timeWindow(5),
  m_upperThreshold(7500)
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
	h_values = new TH1F("values", "values", 1000, 0, 1000);
	h_sizeVsDuration = new TH2F("sizeVsDuration", "sizeVsDuration", 40, -0.5, 39.5,
			m_timeWindow+4, -0.5, m_timeWindow+3.5);
	h_channelRate = new TH1F("channelRate", "channelRate", runner()->nCnG(), -0.5, runner()->nCnG()-0.5);
	TDirectory * instance_direc = runner()->saveFile()->mkdir(name().c_str());
	h_valuesVsChannel = new TH2F("triggerValueVsChannel", "triggerValueVsChannel",
			runner()->nCnG(), -0.5, runner()->nCnG()-0.5, 1000, 0, 10000);
	h_integrals = new TH1F("integrals", "integrals", 1000, 0, 3000);
	h_integralsVsChannel = new TH2F("integralsVsChannel", "integralsVsChannel",
			runner()->nCnG(), 0, runner()->nCnG(), 1000, 0, 3000);
	h_integralsVsValues = new TH2F("integralsVsValues", "integralsVsValues", 1000, 0, 3000,
			1000, 0, 1000);
}


//_____________________________________________________________________________

void SafCoincidenceFinder::threadExecute(unsigned int iGlib, unsigned int iChannelLow,
	unsigned int iChannelUp)
{
	std::vector<long long int> * triggerTimes = runner()->triggerData()->times();
	std::vector<SafRawDataChannel*> * channels = runner()->triggerData()->channels();
	std::vector<SafRawDataChannel*> coinChannels;
	std::vector<double> coinTriggers;
	std::vector<double> coinIntegrals;

	if (triggerTimes->size() < 3) return;

	for (unsigned int iTime = 0; iTime < triggerTimes->size()-1; iTime++) {
		unsigned int size = 1;
		int duration = 0;
		int coinDuration = 0;
		unsigned int plane = channels->at(iTime)->plane();
		bool samePlane = false;
		bool oppositeSide = false;
		unsigned int nBelowUpperThreshold = 0;
		if (runner()->triggerData()->values()->at(iTime) < m_upperThreshold)
			nBelowUpperThreshold++;

		coinChannels.push_back(channels->at(iTime));
		coinTriggers.push_back(runner()->triggerData()->values()->at(iTime));
		coinIntegrals.push_back(runner()->triggerData()->integrals()->at(iTime));
		for (unsigned int jTime = iTime + 1; jTime < triggerTimes->size(); jTime++) {
			duration = triggerTimes->at(jTime) - triggerTimes->at(iTime);
			if (channels->at(iTime) == channels->at(jTime)) continue;
			if (duration < m_timeWindow)	{
				coinDuration = duration;
				size++;
				coinChannels.push_back(channels->at(jTime));
				coinTriggers.push_back(runner()->triggerData()->values()->at(jTime));
				coinIntegrals.push_back(runner()->triggerData()->integrals()->at(jTime));
				if (runner()->triggerData()->values()->at(jTime) < m_upperThreshold)
					nBelowUpperThreshold++;
				if (channels->at(jTime)->plane() == plane) {
					samePlane = true;
					if (channels->at(iTime)->side() != channels->at(jTime)->side())
						oppositeSide = true;
				}
			}
			else break;
		}

		if (size > 1 && samePlane && oppositeSide && nBelowUpperThreshold > 1) {
			h_size->Fill(size);
			h_duration->Fill(coinDuration);
			h_sizeVsDuration->Fill(size, coinDuration);

			for (unsigned int i=0; i<size; i++) {
				h_channelRate->Fill(coinChannels[i]->plotIndex());
				h_valuesVsChannel->Fill(coinChannels[i]->plotIndex(), coinTriggers[i]);
				h_integrals->Fill(coinIntegrals[i]);
				h_integralsVsChannel->Fill(coinChannels[i]->plotIndex(), coinIntegrals[i]);
				h_integralsVsValues->Fill(coinIntegrals[i], coinTriggers[i]);
				h_values->Fill(coinTriggers[i]);
				h_integrals->Fill(coinIntegrals[i]);
			}
			iTime += size;
		}
		coinChannels.clear();
		coinTriggers.clear();
		coinIntegrals.clear();
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
	h_values->Write();
	h_integrals->Write();
	h_sizeVsDuration->Write();
	h_duration->Write();
	h_channelRate->Write();
	h_valuesVsChannel->Write();
	h_integralsVsValues->Write();
	for (int i=0; i<h_channelRate->GetNbinsX(); i++)
		h_channelRate->SetBinContent(i, h_channelRate->GetBinContent(i)/runner()->realTimeElapsed());
}


//_____________________________________________________________________________

