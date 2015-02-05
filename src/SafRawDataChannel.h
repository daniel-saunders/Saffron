/*
 * SafRawDataChannel.h
 *
 *  Created on: Dec 13, 2014
 *      Author: Daniel Saunders
 */


#ifndef SAFRAWDATACHANNEL_H_
#define SAFRAWDATACHANNEL_H_
#include "SafDataSet.h"
#include "SafRawDataSet.h"
#include <vector>


// Forward declarations.
class SafRunner;
class SafRawDataSet;
class SafDataSet;


class SafRawDataChannel : SafDataSet
{
private:
	// Members __________________________________________________________________
	unsigned int m_nEntries;
	unsigned int m_nTriggers;
	unsigned int m_glibID;
	unsigned int m_channelID;
	std::vector<long long int> m_times;
	std::vector<double> m_signals;
	std::vector<unsigned int> m_triggerIDs;
	unsigned int m_nEntriesTotal;
	unsigned int m_nTriggersTotal;
	double m_baseLineEst;
	bool m_baseLineEstSet;
	unsigned int m_nTriggerSamplesWritten;



public:
	// Methods __________________________________________________________________
	SafRawDataChannel(unsigned int glibID, unsigned int channelID,
			SafRawDataSet * rawData, SafRunner * runner);
	virtual ~SafRawDataChannel();
	void clear();
	void calcBaseLineEst();


	// Setters and getters ______________________________________________________
	void addNTriggerSamplesWritten(unsigned int n) {m_nTriggerSamplesWritten += n;}
	unsigned int glibID() {return m_glibID;}
  unsigned int nEntries() {return m_nEntries;}
  unsigned int nTriggers() {return m_nTriggers;}
  bool baseLineEstSet() {return m_baseLineEstSet;}
  double baseLineEst() {
  	if (!m_baseLineEstSet) return -1.0;
  	else return m_baseLineEst;
  }
  void setBaseLineEst(double b) {m_baseLineEst = b;}
  void addNEntries(unsigned int n) {
  	m_nEntries = n;
    m_nEntriesTotal += n;
  }

  void setNEntries(unsigned int n) {
  	m_nEntries = n;
  }

  void addNTriggers(unsigned int n) {
  	m_nTriggers = n;
    m_nTriggersTotal += n;
  }

	unsigned int channelID() {return m_channelID;}
	std::vector<long long int> * times() {return &m_times;}
	std::vector<double> * signals() {return &m_signals;}
	std::vector<unsigned int> * triggerIDs() {return &m_triggerIDs;}
	unsigned int plotIndex();
	unsigned int nEntriesTotal() {return m_nEntriesTotal;}
	unsigned int nTriggersTotal() {return m_nTriggersTotal;}
	unsigned int nTriggerSamplesWritten() {return m_nTriggerSamplesWritten;}
	unsigned int plane() {return m_glibID*2 + (m_channelID/36);}
	unsigned int side() {return m_glibID*4 + (m_channelID/16);}
};

#endif /* SAFRAWDATACHANNEL_H_ */
