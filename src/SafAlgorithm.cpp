/*
 * SafAlgorithm.cpp
 *
 *  Created on: Dec 13, 2014
 *      Author: ds9133
 */

#include "SafAlgorithm.h"

//____________________________________________________________________________

SafAlgorithm::SafAlgorithm(SafRunner * runner, std::string name) :
  m_eof(false)
{
	std::mutex m_mtx;
  m_runner = runner;
  m_childClassName = name;
  m_event = 0;
  m_totalTime = 0;
}


//____________________________________________________________________________

void SafAlgorithm::detailedInitialize()
{
	std::string output = "Initializer of " + m_childClassName;
  std::cout<<output<<std::endl;
  initialize();
}


//____________________________________________________________________________

void SafAlgorithm::detailedExecute()
{
	struct timeval  tv1, tv2;
	gettimeofday(&tv1, NULL);
	
	execute();
	m_event++;
	
	gettimeofday(&tv2, NULL);
	double timeElapsed = (unsigned long long) (tv2.tv_usec - tv1.tv_usec) 
			+ (unsigned long long) ((tv2.tv_sec - tv1.tv_sec)*1000000ULL);
	if (m_event > 1) m_totalTime += timeElapsed;
}


//____________________________________________________________________________

void SafAlgorithm::detailedFinalize()
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