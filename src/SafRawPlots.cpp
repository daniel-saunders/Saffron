/*
 * SafRawPlots.cpp
 *
 *  Created on: Dec 15, 2014
 *      Author: ds9133
 */

#include "SafRawPlots.h"


//_____________________________________________________________________________

SafRawPlots::SafRawPlots(SafRunner * runner, bool filtered) :
  SafAlgorithm(runner, "SafRawPlots"),
  m_diffBinRange(6),
  m_nSeekedRoots(3),
  m_calculateGains(true),
  m_nFinalizeThreads(10),
  m_smoothing(true)
{
	m_filtered = filtered;
}


//_____________________________________________________________________________

SafRawPlots::~SafRawPlots()
{
	// TODO Auto-generated destructor stub
}



//_____________________________________________________________________________

void SafRawPlots::initialize()
{
	m_threading = true;

	std::string direcName = name();
	if (m_filtered) direcName += "-Filtered";

	unsigned int nG = runner()->geometry()->nGlibs();
	unsigned int nC = runner()->geometry()->nChannels();

	h_firstEventWaveforms = initPerChannelPlots("FirstEventWaveForm", "FirstEventWaveForm", 
		runner()->eventTimeWindow(), 0.0, runner()->eventTimeWindow());

	std::string name = "Signal";
	if (m_filtered) name += "-Filtered";
	h_signals = initPerChannelPlots(name.c_str(), name.c_str(), 4500, 7000, 16000);
	name = "SignalDifferential";
	if (m_filtered) name += "-Filtered";
	h_signalsDiff = initPerChannelPlots(name.c_str(), name.c_str(), 4500, 7000, 16000);
	name = "SignalDoubleDifferential";
	if (m_filtered) name += "-Filtered";
	h_signalsDoubleDiff = initPerChannelPlots(name.c_str(), name.c_str(), 4500, 7000, 16000);
	name = "SignalTripleDifferential";
	if (m_filtered) name += "-Filtered";
	h_signalsTripleDiff = initPerChannelPlots(name.c_str(), name.c_str(), 4500, 7000, 16000);


	int nChannels = nC*nG;
	h_allSignals = new TH2F("AllSignalDist", "AllSignalDist", nChannels, -0.5,
			nChannels-0.5, 4500, 7000., 16000.);
	h_signalMeans = new TH1F("SignalMeans", "SignalMeans", nChannels, -0.5, nChannels-0.5);
	h_signalWidths = new TH1F("SignalWidths", "SignalWidths", nChannels, -0.5, nChannels-0.5);
	h_nBaseLineEstVsChannel = new TH1F("BaseLineEstimates", "BaseLineEstimates", nChannels, -0.5, nChannels-0.5);


	// Root file directories.
	TDirectory * instance_direc = runner()->saveFile()->mkdir(direcName.c_str());
	instance_direc->mkdir("FirstWaveforms");
	for (unsigned int i=0; i<nG; i++) {
		std::stringstream ssGlib; ssGlib<<i;
		instance_direc->mkdir(("FirstWaveforms/Glib" + ssGlib.str()).c_str());
	}

	instance_direc->mkdir("Signals");
	for (unsigned int i=0; i<nG; i++) {
		std::stringstream ssGlib; ssGlib<<i;
		instance_direc->mkdir(("Signals/Glib" + ssGlib.str()).c_str());
	}

	instance_direc->mkdir("SignalsDifferentiated");
	for (unsigned int i=0; i<nG; i++) {
		std::stringstream ssGlib; ssGlib<<i;
		instance_direc->mkdir(("SignalsDifferentiated/Glib" + ssGlib.str()).c_str());
	}

	instance_direc->mkdir("SignalsDoubleDifferentiated");
	for (unsigned int i=0; i<nG; i++) {
		std::stringstream ssGlib; ssGlib<<i;
		instance_direc->mkdir(("SignalsDoubleDifferentiated/Glib" + ssGlib.str()).c_str());
	}

	instance_direc->mkdir("SignalsTripleDifferentiated");
	for (unsigned int i=0; i<nG; i++) {
		std::stringstream ssGlib; ssGlib<<i;
		instance_direc->mkdir(("SignalsTripleDifferentiated/Glib" + ssGlib.str()).c_str());
	}


	h_avSignalPerEventPerChannel = new TH2F("AvSignalPerEventPerChannel",
		"AvSignalPerEventPerChannel", nC*nG, 0, nC*nG, runner()->nEvents(), 0, runner()->nEvents());
  h_rmsSignalPerEventPerChannel = new TH2F("RMSSignalPerEventPerChannel",
		"RMSSignalPerEventPerChannel", nC*nG, 0, nC*nG, runner()->nEvents(), 0, runner()->nEvents());

  h_avSignalOfEvents = new TH1F("AvSignalPerEvent", "AvSignalPerEvent", 500, 0.998, 1.002);
  h_avSignalPerChannelOfEvents = new TH2F("AvSignalPerEvent", "AvSignalPerEvent",
  		nC*nG, 0., 1.*nC*nG, 500, 0.998, 1.002);

  std::string nameX = "Gains";
  if (m_filtered) nameX += "-Filtered";
  h_gains = new TH1F(nameX.c_str(), nameX.c_str(), 80, 0., 100.);

  nameX = "GainsPerChannel";
  if (m_filtered) nameX += "-Filtered";
  h_gainsPerChannel = new TH1F(nameX.c_str(), nameX.c_str(), runner()->nCnG(), 0., 1.*runner()->nCnG());
}


//_____________________________________________________________________________

void SafRawPlots::execute()
{
	for (unsigned int i=0; i<runner()->geometry()->nGlibs(); i++){
		threadExecute(i, 0, runner()->geometry()->nChannels(), -1);
	}
}


//_____________________________________________________________________________

void SafRawPlots::finalize()
{
//	if (m_calculateGains && m_threading && !m_forceSingleThread) {
//		unsigned int step = h_signals->size()/m_nFinalizeThreads;
//		unsigned int iUp = step;
//		std::vector<std::thread> threads;
//		while (iUp < h_signals->size()) {
//			std::cout<<"Starting thread, with values: "<<iUp-step<<"\t"<<iUp<<std::endl;
//			threads.push_back(std::thread(&SafRawPlots::calculateGains, this, iUp-step, iUp));
//			iUp += step;
//		}
//		for (std::vector<std::thread>::iterator i = threads.begin();
//			i != threads.end(); i++)
//			(*i).join();
//	}
//	else if (m_calculateGains)

	if (m_calculateGains) calculateGains(0, 20);


	std::string direcName = name();
	if (m_filtered) direcName += "-Filtered";
	for (unsigned int i=0; i<h_firstEventWaveforms->size(); i++) {
		int iGlib = i/runner()->geometry()->nChannels();
		std::stringstream ssGlib; ssGlib << iGlib;
		runner()->saveFile()->cd((direcName + "/FirstWaveforms/Glib" + ssGlib.str()).c_str());

		h_firstEventWaveforms->at(i)->Write();
	}

	for (unsigned int i=0; i<h_signals->size(); i++) {
		for (int j=0; j<h_signals->at(i)->GetNbinsX(); j++) {
			h_allSignals->SetBinContent(i, j, h_signals->at(i)->GetBinContent(j));
		}
		int iGlib = i/runner()->geometry()->nChannels();
		std::stringstream ssGlib; ssGlib << iGlib;
		runner()->saveFile()->cd((direcName + "/Signals/Glib" + ssGlib.str()).c_str());
		h_signals->at(i)->Write();

		int fitStatus = h_signals->at(i)->Fit("gaus", "Q");
		if (fitStatus != 0) continue;
		h_signalMeans->SetBinContent(i, h_signals->at(i)->GetFunction("gaus")->GetParameter(1));
		h_signalMeans->SetBinError(i, h_signals->at(i)->GetFunction("gaus")->GetParError(1));

		h_signalWidths->SetBinContent(i, h_signals->at(i)->GetFunction("gaus")->GetParameter(2));
		h_signalWidths->SetBinError(i, h_signals->at(i)->GetFunction("gaus")->GetParError(2));

		runner()->saveFile()->cd((direcName + "/SignalsDifferentiated/Glib" + ssGlib.str()).c_str());
		h_signalsDiff->at(i)->Write();

		runner()->saveFile()->cd((direcName + "/SignalsDoubleDifferentiated/Glib" + ssGlib.str()).c_str());
		h_signalsDoubleDiff->at(i)->Write();

		runner()->saveFile()->cd((direcName + "/SignalsTripleDifferentiated/Glib" + ssGlib.str()).c_str());
		h_signalsTripleDiff->at(i)->Write();
	}

	for (int i=0; i<h_avSignalPerEventPerChannel->GetNbinsX(); i++) {
		for (int j=0; j<h_avSignalPerEventPerChannel->GetNbinsY(); j++) {
			double av = h_avSignalPerEventPerChannel->GetBinContent(i, j);
			double rms = h_rmsSignalPerEventPerChannel->GetBinContent(i, j);
			av /= h_signals->at(i)->GetMean();
			rms /= h_signals->at(i)->GetMean();
			h_avSignalPerEventPerChannel->SetBinContent(i, j, av);
			h_rmsSignalPerEventPerChannel->SetBinContent(i, j, rms);
			h_avSignalOfEvents->Fill(av);
			h_avSignalPerChannelOfEvents->Fill(i, av);
		}
	}

	for (unsigned int i=0; i<runner()->geometry()->nChannels(); i++) {
		for (unsigned int j=0; j<runner()->geometry()->nGlibs(); j++) {
			SafRawDataChannel * channel = runner()->rawData()->channel(j, i);
  		h_nBaseLineEstVsChannel->SetBinContent(channel->plotIndex(), channel->baseLineEst());
		}
	}

	runner()->saveFile()->cd(direcName.c_str());
	h_allSignals->Write();
	h_signalMeans->Write();
	h_signalWidths->Write();
	h_avSignalPerEventPerChannel->Write();
	h_rmsSignalPerEventPerChannel->Write();
	h_avSignalOfEvents->Write();
	h_avSignalPerChannelOfEvents->Write();
	h_gains->Write();
	h_gainsPerChannel->Write();
	h_nBaseLineEstVsChannel->Write();

	runner()->saveFile()->cd();
	std::cout<<"REF: "<<h_signals->at(0)->GetMean()<<std::endl;
}


//_____________________________________________________________________________

void SafRawPlots::calculateGains(unsigned int iLow, unsigned int iUp) {
	std::vector<TH1F*>::iterator ih;
	int iPlot = iLow;
	std::cout<<"Calculating gains - this can take some minutes..."<<std::endl;

	for (ih = (h_signals->begin()+iLow); ih!=(h_signals->begin()+iUp); ih++) {
		if (iPlot % 10 == 0) std::cout<<"Progress: "<<iPlot<<"\t/"<<iUp*4<<std::endl;

		if (m_smoothing) smooth(*ih, iPlot);
		int istart = h_signals->at(iPlot)->GetMaximumBin();
	  int iend = istart + 200;
		for (int i=istart; i<iend; i++) {
			double rangeLow = (*ih)->GetBinCenter(i);
			double rangeHigh = (*ih)->GetBinCenter(i+m_diffBinRange);
			std::stringstream ss;
			ss<<i;
			std::string name = "Signal" + ss.str();
			m_mtx.lock();
			TF1 * fit = new TF1(name.c_str(), "pol2", rangeLow, rangeHigh);
			int fitStatus = (*ih)->Fit(name.c_str(), "RQ");
			m_mtx.unlock();
			double x = (*ih)->GetBinCenter(i+m_diffBinRange/2);
			if (fitStatus == 0) {
				h_signalsDiff->at(iPlot)->SetBinContent(i+m_diffBinRange/2, fit->Derivative(x));
			}
			delete fit;
		}
		iPlot++;
	}

	iPlot = iLow;
	for (ih = (h_signalsDiff->begin()+iLow); ih!=(h_signalsDiff->begin()+iUp); ih++) {
		if (iPlot % 10 == 0) std::cout<<"Progress: "<<iPlot+iUp<<"\t/"<<iUp*4<<std::endl;
		int istart = h_signals->at(iPlot)->GetMaximumBin();
	  int iend = istart + 200;
		for (int i=istart; i<iend; i++) {
			double rangeLow = (*ih)->GetBinCenter(i);
			double rangeHigh = (*ih)->GetBinCenter(i+m_diffBinRange);
			std::stringstream ss;
			ss<<i;
			std::string name = "SignalDiff" + ss.str();
			m_mtx.lock();
			TF1 * fit = new TF1(name.c_str(), "pol2", rangeLow, rangeHigh);
			int fitStatus = (*ih)->Fit(name.c_str(), "RQ");
			m_mtx.unlock();
			double x = (*ih)->GetBinCenter(i+m_diffBinRange/2);
			if (fitStatus == 0) {
				h_signalsDoubleDiff->at(iPlot)->SetBinContent(i+m_diffBinRange/2, fit->Derivative(x));
			}
			delete fit;
		}
		iPlot++;
	}

	iPlot = iLow;
	for (ih = (h_signalsDoubleDiff->begin()+iLow); ih!=(h_signalsDoubleDiff->begin()+iUp); ih++) {
		if (iPlot % 10 == 0) std::cout<<"Progress: "<<iPlot + 2*iUp<<"\t/"<<iUp*4<<std::endl;
		int istart = h_signals->at(iPlot)->GetMaximumBin();
	  int iend = istart + 200;
		for (int i=istart; i<iend; i++) {
			double rangeLow = (*ih)->GetBinCenter(i);
			double rangeHigh = (*ih)->GetBinCenter(i+m_diffBinRange);
			std::stringstream ss;
			ss<<i;
			std::string name = "SignalDiff" + ss.str();
			m_mtx.lock();
			TF1 * fit = new TF1(name.c_str(), "pol2", rangeLow, rangeHigh);
			int fitStatus = (*ih)->Fit(name.c_str(), "RQ");
			m_mtx.unlock();
			double x = (*ih)->GetBinCenter(i+m_diffBinRange/2);
			if (fitStatus == 0) {
				h_signalsTripleDiff->at(iPlot)->SetBinContent(i+m_diffBinRange/2, fit->Derivative(x));
			}
			delete fit;
		}
		iPlot++;
	}

	iPlot = iLow;
	for (ih = (h_signalsTripleDiff->begin()+iLow); ih!=(h_signalsTripleDiff->begin()+iUp); ih++) {
		bool gainSet = false;
		if (iPlot % 10 == 0) std::cout<<"Progress 4: "<<iPlot + 3*iUp<<"\t/"<<iUp*4<<std::endl;
		std::vector<double> roots;
		int istart = h_signals->at(iPlot)->GetMaximumBin();
	  int iend = istart + 500;
		for (int i=istart; i<iend; i++) {
			if ((*ih)->GetBinContent(i) < 0 &&
					(*ih)->GetBinContent(i+1) < 0 &&
					(*ih)->GetBinContent(i+2) > 0 &&
					(*ih)->GetBinContent(i+3) > 0) {
//				std::stringstream ss;
//				ss<<i;
//				std::string name = "Fit" + ss.str();
//				m_mtx.lock();
//				TF1 * fit = new TF1(name.c_str(), "pol1", (*ih)->GetBinCenter(i), (*ih)->GetBinCenter(i+4));
//				int fitStatus = (*ih)->Fit(name.c_str(), "RQ");
//				m_mtx.unlock();
				double root = (*ih)->GetBinLowEdge(i+2);
//				if (fitStatus == 0) root = -fit->GetParameter(0)/fit->GetParameter(1);
//				else root = (*ih)->GetBinLowEdge(i+2);
				roots.push_back(root);
				//delete fit;
				if (roots.size() >= m_nSeekedRoots)  {
					if (roots.size() == m_nSeekedRoots) {
						std::vector<double> seps;
						for (unsigned int j=0; j<roots.size()-1; j++)
							seps.push_back(roots[j+1] - roots[j]);

						double gain = 0;
						for (unsigned int j=0; j<seps.size(); j++)
							gain += seps[j];

						gain/=(1.*seps.size());
						m_mtx.lock();
						h_gainsPerChannel->SetBinContent(iPlot, gain);
						h_gains->Fill(gain);
						std::cout<<"Found a gain: "<<gain<<"\t"<<iPlot<<std::endl;
						gainSet = true;
						m_mtx.unlock();
						i+=10;
					}
				}
			}
			if (gainSet) break;
		}
		iPlot++;
	}

	std::cout<<"Gain calculating complete. About 1 min left."<<std::endl;
}


//_____________________________________________________________________________

void SafRawPlots::threadExecute(unsigned int iGlib, unsigned int iLow, 
	unsigned int iUp, int iThread)
{
	unsigned int nChannels = runner()->geometry()->nChannels();
	for (unsigned int i=iLow; i<iUp; i++) {
		double sum = 0.0;
		double sumSq = 0.0;
		SafRawDataChannel * channel = runner()->rawData()->channel(iGlib, i);
		unsigned int plotIndex = iGlib*nChannels + i;

		for (unsigned int k=0; k<channel->nEntries(); k++) {
			if (runner()->event() == 0) {
				h_firstEventWaveforms->at(plotIndex)->SetBinContent(
						channel->times()->at(k), channel->signals()->at(k));
			}

			double signal = channel->signals()->at(k);
			h_signals->at(plotIndex)->Fill(signal);
			sum += signal;
			sumSq += signal*signal;
		}
		double av = sum/(1.*channel->nEntries());
		double rms = pow(sumSq/(1.*channel->nEntries()) - av*av, 0.5);


		unsigned int iChannel = i+iGlib*nChannels;
		m_mtx.lock();
		h_avSignalPerEventPerChannel->SetBinContent(iChannel, m_event, av);
		h_rmsSignalPerEventPerChannel->SetBinContent(iChannel, m_event, rms);
		m_mtx.unlock();
	}
}


//_____________________________________________________________________________

void SafRawPlots::smooth(TH1F * h, int plotIndex) {
	if (plotIndex == 379) return;
	int istart = 10;
	if (h->GetMaximumBin() > istart) istart = h->GetMaximumBin();
	int iend = istart + 500;
	for (unsigned int i=istart; i<iend; i++) {
		double content = h->GetBinContent(i);
		content += h->GetBinContent(i+1);
		content += h->GetBinContent(i-1);
		content += h->GetBinContent(i+2);
		content += h->GetBinContent(i-2);
		content += h->GetBinContent(i+3);
		content += h->GetBinContent(i-3);
		content /= 7.0;
		h->SetBinContent(i, content);

		if (iend > 2000) break;
	}
}


//_____________________________________________________________________________

