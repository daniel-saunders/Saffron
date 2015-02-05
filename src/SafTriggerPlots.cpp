/*
 * SafTriggerPlots.cpp
 *
 *  Created on: Dec 30, 2014
 *      Author: Daniel Saunders
 */

#include "SafTriggerPlots.h"


//_____________________________________________________________________________

SafTriggerPlots::SafTriggerPlots(SafRunner * runner) :
  SafAlgorithm(runner, "SafTriggerPlots")
{
	// TODO Auto-generated constructor stub

}


//_____________________________________________________________________________

SafTriggerPlots::~SafTriggerPlots()
{
	// TODO Auto-generated destructor stub
}


//_____________________________________________________________________________

void SafTriggerPlots::initialize()
{
	unsigned int nG = runner()->geometry()->nGlibs();
	unsigned int nC = runner()->geometry()->nChannels();
	
	h_values = new TH1F("TriggerValues", "TriggerValues", 1000, 0, 10000);
	h_valuesZoomed = new TH1F("TriggerValuesZoomed", "TriggerValuesZoomed", 1000, 0, 1000);
	h_integrals = new TH1F("TriggerIntegral", "TriggerIntegral", 1000, 0, 3000);
	h_integralsVsValues = new TH2F("TriggerIntegralVsValue", "TriggerIntegralVsValue",
			1000, 0, 1000, 1000, 0, 3000);
	
	h_integralsVsValues0to38 = new TH2F("TriggerIntegralVsValue0to38", "TriggerIntegralVsValue0to38",
			1000, 0, 1000, 1000, 0, 3000);

	h_integralsVsValues38to76 = new TH2F("TriggerIntegralVsValue38to76", "TriggerIntegralVsValue38to76",
			1000, 0, 1000, 1000, 0, 3000);

	h_integralsVsValuesOdd = new TH2F("TriggerIntegralVsValueOdd", "TriggerIntegralVsValueOdd",
			1000, 0, 1000, 1000, 0, 3000);

	h_integralsVsValuesEven = new TH2F("TriggerIntegralVsValueEven", "TriggerIntegralVsValueEven",
			1000, 0, 1000, 1000, 0, 3000);

  h_firstEventPeaks = initPerChannelPlots("FirstEventPeaks", "FirstEventPeaks", 
		runner()->eventTimeWindow(), 0.0, runner()->eventTimeWindow());
	h_valuesPerChannel = initPerChannelPlots("TriggerValues", "TriggerValues", 500, 100, 10000);
	
	h_nTriggers = new TH1F("AverageTriggerRate", "AverageTriggerRate", nC*nG, 
			0, nC*nG);

	int nChannels = nC*nG;
	h_valuesVsChannel = new TH2F("AllTriggerDist", "AllTriggerDist", nChannels,
			-0.5, nChannels-0.5, 500, 100, 2500);

	h_integralsVsChannel = new TH2F("TriggerIntVsChannel", "TriggerIntVsChannel", nChannels,
		-0.5, nChannels-0.5, 1000, 0, 3000);


	// Root file directories.
	TDirectory * instance_direc = runner()->saveFile()->mkdir(name().c_str());
	instance_direc->mkdir("FirstPeaks");
	for (unsigned int i=0; i<nG; i++) {
		std::stringstream ssGlib; ssGlib<<i;
		instance_direc->mkdir(("FirstPeaks/Glib" + ssGlib.str()).c_str());
	}

	instance_direc->mkdir("TriggerValues");
	for (unsigned int i=0; i<nG; i++) {
		std::stringstream ssGlib; ssGlib<<i;
		instance_direc->mkdir(("TriggerValues/Glib" + ssGlib.str()).c_str());
	}

	h_dataRates = new TH1F("TriggerDataRate", "TriggerDataRate; ChannelID; SamplesWrittenPerSec", nChannels, -0.5, nChannels-0.5);
	h_triggerRate = new TH1F("TriggerRate", "TriggerRate; ChannelID; TriggerRate(Hz)", nChannels, -0.5, nChannels-0.5);
	h_nTriggersVsEvents = new TH2F("nTriggersVsEvents", "nTriggersVsEvents",
			nChannels, -0.5, nChannels-0.5, runner()->nEvents(), 0, runner()->nEvents());
}


//_____________________________________________________________________________

void SafTriggerPlots::execute()
{
	fill();
}


//_____________________________________________________________________________

void SafTriggerPlots::fill()
{
	SafTriggerDataSet * data = runner()->triggerData();
	for (unsigned int i=0; i<data->nTriggers(); i++) {
		unsigned int plotIndex = data->channels()->at(i)->plotIndex();
		h_values->Fill(data->values()->at(i));
		h_valuesZoomed->Fill(data->values()->at(i));
		h_integrals->Fill(data->integrals()->at(i));
		h_integralsVsValues->Fill(data->values()->at(i),
				data->integrals()->at(i));

		if (data->channels()->at(i)->channelID() < 38)
			h_integralsVsValues0to38->Fill(data->values()->at(i),
				data->integrals()->at(i));
		else
			h_integralsVsValues38to76->Fill(data->values()->at(i),
				data->integrals()->at(i));

		if (data->channels()->at(i)->channelID() % 2 == 0)
			h_integralsVsValuesEven->Fill(data->values()->at(i),
				data->integrals()->at(i));
		else
			h_integralsVsValuesOdd->Fill(data->values()->at(i),
				data->integrals()->at(i));

		h_valuesPerChannel->at(plotIndex)->Fill(data->values()->at(i));
		h_valuesVsChannel->Fill(plotIndex, data->values()->at(i));
		h_integralsVsChannel->Fill(plotIndex, data->integrals()->at(i));
		if (runner()->event() == 0)
			h_firstEventPeaks->at(plotIndex)->Fill(data->times()->at(i), data->values()->at(i));
	}


	for (unsigned int i=0; i<runner()->geometry()->nChannels(); i++) {
		for (unsigned int j=0; j<runner()->geometry()->nGlibs(); j++) {
			SafRawDataChannel * channel = runner()->rawData()->channel(j, i);
			h_nTriggersVsEvents->SetBinContent(channel->plotIndex(), m_event, channel->nTriggers());
		}
	}
}


//_____________________________________________________________________________

void SafTriggerPlots::finalize()
{
	for (unsigned int i=0; i<h_firstEventPeaks->size(); i++) {
		int iGlib = i/runner()->geometry()->nChannels();
		std::stringstream ssGlib; ssGlib << iGlib;
		runner()->saveFile()->cd((name() + "/FirstPeaks/Glib" + ssGlib.str()).c_str());

		h_firstEventPeaks->at(i)->Write();
	}

	for (unsigned int i=0; i<h_valuesPerChannel->size(); i++) {
		int iGlib = i/runner()->geometry()->nChannels();
		std::stringstream ssGlib; ssGlib << iGlib;
		runner()->saveFile()->cd((name() + "/TriggerValues/Glib" + ssGlib.str()).c_str());

		h_valuesPerChannel->at(i)->Write();
	}
	
	unsigned int nGlibs = runner()->geometry()->nGlibs();
	unsigned int nChannels = runner()->geometry()->nChannels();
	
	for (unsigned int i=0; i<nGlibs; i++) {
		for (unsigned int j=0; j<nChannels; j++) {
			SafRawDataChannel * channel = runner()->rawData()->channel(i, j);
			unsigned int plotIndex = channel->plotIndex();
			h_nTriggers->SetBinContent(plotIndex, channel->nTriggersTotal()/(1.*runner()->nEvents()));
			double triggerRate = channel->nTriggersTotal()/runner()->realTimeElapsed();
			h_triggerRate->SetBinContent(plotIndex, triggerRate);
			double dataRate = channel->nTriggerSamplesWritten()/runner()->realTimeElapsed();
			h_dataRates->SetBinContent(plotIndex, dataRate);
		}
	}

	for (unsigned int i=0; i<h_valuesPerChannel->size(); i++){
		for (unsigned int j=0; j<h_valuesPerChannel->at(i)->GetNbinsX(); j++) {
			h_valuesVsChannel->SetBinContent(i, j,
					h_valuesPerChannel->at(i)->GetBinContent(j)/runner()->realTimeElapsed());
		}
	}

	
	runner()->saveFile()->cd(name().c_str());
	h_integralsVsValues->Write();
	h_values->Write();
	h_integralsVsChannel->Write();
	h_nTriggers->Write();
	h_valuesVsChannel->Write();
	h_dataRates->Write();
	h_triggerRate->Write();
	h_nTriggersVsEvents->Write();
	h_integrals->Write();
	h_valuesZoomed->Write();
	h_integralsVsValues0to38->Write();
	h_integralsVsValues38to76->Write();
	h_integralsVsValuesOdd->Write();
	h_integralsVsValuesEven->Write();
}


//_____________________________________________________________________________


