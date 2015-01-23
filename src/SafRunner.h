/*
 * SafRunner.h
 *
 *  Created on: Dec 13, 2014
 *      Author: Daniel Saunders
 *        Desc: Runner of Saffron, including various options.
 */


#ifndef SAFRUNNER_H_
#define SAFRUNNER_H_
#include "SafAlgorithm.h"
#include "SafEventBuilder.h"
#include "SafTrigger.h"
#include "SafRawDataSet.h"
#include "SafTriggerDataSet.h"
#include "SafGeometry.h"
#include <vector>
#include <iostream>
#include <string>
#include "TFile.h"
#include "SafRawPlots.h"
#include "SafTriggerPlots.h"
#include "SafPeakFitter.h"
#include "SafFilter.h"
#include "SafCoincidenceFinder.h"


// Forward declarations.
class SafAlgorithm;
class SafRawDataSet;
class SafTriggerDataSet;
class SafEventBuilder;


class SafRunner
{
private:
  // Members __________________________________________________________________
	unsigned int m_printThreshold; // Sets the level of print outs.
	std::vector< SafAlgorithm* > m_algorithms;
	SafRawDataSet * m_rawData; // Set by the initializer of the event builder.
	SafTriggerDataSet * m_triggerData;

	// SafPeakSet * m_peaks;
	unsigned int m_eventTimeWindow; // Units of 16 ns.
	unsigned int m_event;
	unsigned int m_nEvents; // Number of events to loop over.
	SafGeometry * m_geometry;
	unsigned int m_timeZero;
	unsigned int m_runMode;
	TFile * m_saveFile;
	unsigned int m_printRate;
	unsigned int m_threadsPerGlib;
	unsigned int m_triggerSkip;
	std::string m_fileName;
	std::string m_saveFileName;


public:
	double triggerThreshold;
	// Methods __________________________________________________________________
	SafRunner();
	virtual ~SafRunner();
	void run();
	void safPrint(std::string output, int printLevel);
	void eventLoop();
	double realTimeElapsed() {return 1.*(m_event-m_triggerSkip)*2048*16e-9;}


	// Setters and getters ______________________________________________________
	void setFileName(std::string s) {m_fileName = s;}
	void setSaveFileName(std::string s) {m_saveFileName = s;}
	std::string fileName() {return m_fileName;}
	unsigned int printThreshold() {return m_printThreshold;}
	void setPrintThreshold(unsigned int printThreshold) {
		m_printThreshold = printThreshold;}

	SafRawDataSet * rawData() {return m_rawData;}
	SafTriggerDataSet * triggerData() {return m_triggerData;}
	void setRawData(SafRawDataSet * rawData) {m_rawData = rawData;}
	void setTriggerData(SafTriggerDataSet * triggerData) {m_triggerData = triggerData;}

//	SafPeakSet * peaks() {return m_peaks;}
	SafGeometry * geometry() {return m_geometry;}
	unsigned int eventTimeWindow() {return m_eventTimeWindow;}
	unsigned int runMode() {return m_runMode;}
	unsigned int event() {return m_event;}
	TFile * saveFile() {return m_saveFile;}
	unsigned int nEvents() {return m_nEvents;}
	unsigned int triggerSkip() {return m_triggerSkip;}
	unsigned int nCnG() {return geometry()->nChannels() * geometry()->nGlibs();}
};

#endif /* SAFRUNNER_H_ */
