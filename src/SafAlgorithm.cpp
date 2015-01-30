/*
 * SafAlgorithm.cpp
 *
 *  Created on: Dec 13, 2014
 *      Author: ds9133
 */

#include "SafAlgorithm.h"

//____________________________________________________________________________

SafAlgorithm::SafAlgorithm(SafRunner * runner, std::string name) :
  m_eof(false),
  m_nThreadsPerGlib(2),
  m_threading(false),
  m_forceSingleThread(false)
{
	std::mutex m_mtx;
  m_runner = runner;
  m_childClassName = name;
  m_event = 0;
  m_totalTime = 0;
}


//____________________________________________________________________________

void SafAlgorithm::parentInitialize(unsigned int nGlibs, unsigned int nChannels)
{
  m_nGlibs = nGlibs;
  m_nChannels = nChannels;
	std::string output = "Initializer of " + m_childClassName;
  std::cout<<output<<std::endl;
  m_nChannelsPerThread = m_nChannels/m_nThreadsPerGlib;
  initialize();
}


//____________________________________________________________________________

void SafAlgorithm::parentExecute()
{
	struct timeval  tv1, tv2;
	gettimeofday(&tv1, NULL);
	
	preExecute();
  if (m_threading && !m_forceSingleThread) parentThreadExecute();
	else  execute();
  postExecute();
	m_event++;
	
	gettimeofday(&tv2, NULL);
	double timeElapsed = (unsigned long long) (tv2.tv_usec - tv1.tv_usec) 
			+ (unsigned long long) ((tv2.tv_sec - tv1.tv_sec)*1000000ULL);
	if (m_event > 1) m_totalTime += timeElapsed;
}


//____________________________________________________________________________

void SafAlgorithm::parentThreadExecute()
{
  std::vector<std::thread> threads;
  int ithread = 0;
  for (unsigned int i=0; i<m_nGlibs; i++) {
    for (unsigned int j=0; j<m_nChannels; j+=m_nChannelsPerThread) {
      threads.push_back(std::thread(&SafAlgorithm::threadExecute, this, i, j, 
        j+m_nChannelsPerThread, ithread));
      ithread++;
    }
  }

  for (std::vector<std::thread>::iterator i = threads.begin(); 
    i != threads.end(); i++)
    (*i).join();
}


//____________________________________________________________________________

void SafAlgorithm::parentFinalize()
{
	finalize();
	m_avTime = m_totalTime/(1000.*(m_event-1));

  // Save plots.
	for (std::vector<TH1F*>::iterator iPlot = h_th1sSLOW.begin();
  		iPlot != h_th1sSLOW.end(); iPlot++) (*iPlot)->Write();

	for (std::vector<TH2F*>::iterator iPlot = h_th2sSLOW.begin();
  		iPlot != h_th2sSLOW.end(); iPlot++) (*iPlot)->Write();
}


//____________________________________________________________________________

SafAlgorithm::~SafAlgorithm()
{
	// TODO Auto-generated destructor stub
}


//____________________________________________________________________________

void SafAlgorithm::plot(double x, std::string name, std::string title,
			unsigned int nBinsX, double xLow, double xHigh)
{
	TH1F * plot = findPlot(name);
	if (plot != NULL) {
		// Plot already exists.
		plot->Fill(x);
	}
	else {
		plot = new TH1F(name.c_str(), title.c_str(), nBinsX, xLow, xHigh);
		plot->Fill(x);
		h_th1sSLOW.push_back(plot);
	}
}


//____________________________________________________________________________

TH1F * SafAlgorithm::findPlot(std::string name)
{
  for (std::vector<TH1F*>::iterator iPlot = h_th1sSLOW.begin();
  		iPlot != h_th1sSLOW.end(); iPlot++) {
  	if (std::string((*iPlot)->GetName()) == name) {
  		return (*iPlot);
  	}
  }

  // No plot found.
  return NULL;
}



//____________________________________________________________________________

std::vector<TH1F*> * SafAlgorithm::initPerChannelPlots(std::string name, 
  std::string title, unsigned int nBins, double xLow, double xUp) 
{
  std::vector<TH1F*> * plots = new std::vector<TH1F*>();
  for (unsigned int i=0; i<m_nGlibs; i++) {
    std::stringstream ssGlib; ssGlib<<i;
    for (unsigned int j=0; j<m_nChannels; j++) {
      std::stringstream ssChan; ssChan<<j;

      // First set peaks.
      std::string plotName = name + ssChan.str() + "-" + ssGlib.str();
      std::string plotTitle = title + ssChan.str() + "-" + ssGlib.str();
      TH1F * plot = new TH1F(plotName.c_str(), plotTitle.c_str(),
          nBins, xLow, xUp);
      plots->push_back(plot);
    }
  }
  return plots;
}
