/*
 * SafTriggerDataSet.h
 *
 *  Created on: Jan 13, 2014
 *      Author: Daniel Saunders
 *        Desc: Trigger data container. Implemented in a data oriented way,
 *              as opposed to object oriented, for performance reasons.
 *              Created by the initlizer of the event builder.
 */


#ifndef SAFTRIGGERDATASET_H_
#define SAFTRIGGERDATASET_H_
#include "SafRunner.h"
#include "SafDataSet.h"
#include "SafRawDataChannel.h"


// Forward declarations.
class SafRunner;
class SafRawDataChannel;


class SafTriggerDataSet : SafDataSet
{
private:
	// NOT SORTED CORRECTLY.
	std::vector<long long int> m_times;
	std::vector<SafRawDataChannel*> m_channels;
	std::vector<double> m_values;
	std::vector<double> m_integrals;


public:
	// Methods __________________________________________________________________
	SafTriggerDataSet(SafRunner * runner);
	virtual ~SafTriggerDataSet();

	void clear();
	unsigned int nTriggers() {return m_times.size();}

	// Setters and getters ______________________________________________________
	std::vector<long long int> * times() {return &m_times;}
	std::vector<SafRawDataChannel*> * channels() {return &m_channels;}
	std::vector<double> * values() {return &m_values;}
	std::vector<double> * integrals() {return &m_integrals;}
};

#endif /* SAFTRIGGERDATASET_H_ */
