/*
 * SafCoincidenceFinder.h
 *
 *  Created on: Dec 22, 2014
 *      Author: Daniel Saunders
 */

#ifndef SAFCOINCIDENCEFINDER_H_
#define SAFCOINCIDENCEFINDER_H_

#include "SafAlgorithm.h"
#include "SafRunner.h"
#include "SafRawDataChannel.h"


class SafRawDataChannel;

class SafCoincidenceFinder: public SafAlgorithm
{
private:
	// Members __________________________________________________________________
	TH1F * h_size;
	TH1F * h_duration;
	int m_timeWindow;
	TH1F * h_values;
	TH2F * h_sizeVsDuration;
	TH1F * h_channelRate;
	TH2F * h_valuesVsChannel;
	TH1F * h_integrals;
	TH2F * h_integralsVsChannel;
	double m_upperThreshold;
	TH2F * h_integralsVsValues;

public:
  // Methods __________________________________________________________________
  SafCoincidenceFinder(SafRunner * runner);
  virtual ~SafCoincidenceFinder();

  void initialize();
  void execute();
  void finalize();
  void threadExecute(unsigned int iGlib, unsigned int iLow, unsigned int iUp);
};

#endif /* SAFCOINCIDENCEFINDER_H_ */
