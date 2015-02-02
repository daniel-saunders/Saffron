/*
 * SafRunner.cpp
 *
 *  Created on: Dec 13, 2014
 *      Author: Daniel Saunders
 */

#include "SafRunner.h"

//_____________________________________________________________________________

SafRunner::SafRunner() :
  m_printThreshold(0),
  m_algorithms(),
  m_rawData(NULL),
  m_triggerData(NULL),
  m_eventTimeWindow(2048),
  m_event(0),
  m_timeZero(0),
  m_printRate(10),
  m_triggerSkip(0),
  triggerThreshold(200),
  m_saveFileName("Saffron-histos.root"),
  m_fileNamesPassed(false)
{
	// Default file name (removed if arguments passed in).
	m_rawDataFileNames.push_back("/storage/SOLID/SM1_06Jan2015_1023_run0_scoperun_slowcontrol-small.root");
	// Default algorithm list.
	m_algorithms.push_back(new SafEventBuilder(this));
	m_algorithms.push_back(new SafRawPlots(this, false));
//  m_algorithms.push_back(new SafFilter(this));
//  m_algorithms.push_back(new SafRawPlots(this, true));
	m_algorithms.push_back(new SafTrigger(this));
	m_algorithms.push_back(new SafTriggerPlots(this));
//	m_algorithms.push_back(new SafPeakFitter(this));
	m_algorithms.push_back(new SafCoincidenceFinder(this));

	// Geometry.
	m_geometry = new SafGeometry();

	// Options.
	m_nEvents = 8000; // Total over all input files. EOF is safe.
	m_runMode = 1; // 0 for MC, 1 for real data.
}


//_____________________________________________________________________________

void SafRunner::evalArg(std::string arg) {
	if (arg.substr(0, 14) == "--rawDataFile=") {
		if (!m_fileNamesPassed) m_rawDataFileNames.clear();
		m_fileNamesPassed = true;
		std::string fileName = arg.substr(14, arg.size());
		std::cout<<"Adding file to analysis: "<<fileName<<std::endl;
		m_rawDataFileNames.push_back(fileName);
	}

	else if (arg.substr(0, 19) == "--rawDataDirectory=") {
		if (!m_fileNamesPassed) m_rawDataFileNames.clear();
		m_fileNamesPassed = true;
		std::string dirName = arg.substr(19, arg.size());
		std::cout<<"Opening directory for analysis: "<<dirName<<std::endl;
		DIR *dir;
		struct dirent *ent;
		if ((dir = opendir (dirName.c_str())) != NULL) {
			while ((ent = readdir (dir)) != NULL) {
				std::cout<<"Found file:\t"<<ent->d_name<<std::endl;
				std::string fileName(ent->d_name);
				std::cout<<"Adding file to analysis: "<<fileName<<std::endl;
				m_rawDataFileNames.push_back(dirName + fileName);
			}
			closedir (dir);
		}
		else std::cout<<"Could not open directory: "<<dirName<<std::endl;
	}

	else if (arg.substr(0, 22) == "--outputHistoFileName=") {
		std::string saveFileName = arg.substr(22, arg.size());
		std::cout<<"Setting save file name as: "<<saveFileName<<std::endl;
		setSaveFileName(saveFileName);
	}

	else {
		std::cout<<"Unknown argument: "<<arg<<std::endl;
		exit(0);
	}
}


//_____________________________________________________________________________

SafRunner::~SafRunner()
{
	m_saveFile->Close();
}


//_____________________________________________________________________________

void SafRunner::safPrint(std::string text, int level)
{
	std::cout<<text<<std::endl;
}


//_____________________________________________________________________________

void SafRunner::run()
{
	// Save file.
	m_saveFile = new TFile(m_saveFileName.c_str(), "RECREATE");

	// Default runner of Saffron. In turn, calls all initialise, execute and
	// finalize methods of all SafAlgorithms used in this run.

	new SafRawDataSet(this);
	new SafTriggerDataSet(this);

	std::cout<<"Initializing all algorithms..."<<std::endl;
	std::vector< SafAlgorithm* >::iterator ialgo;
	for (ialgo = m_algorithms.begin(); ialgo != m_algorithms.end(); ialgo++)
		(*ialgo)->parentInitialize(geometry()->nGlibs(), geometry()->nChannels());

	struct timeval  tv1, tv2;
	gettimeofday(&tv1, NULL);

	std::cout<<"Executing event loop."<<std::endl;
	eventLoop();

	gettimeofday(&tv2, NULL);
	double totExTime = (unsigned long long) (tv2.tv_usec - tv1.tv_usec) 
			+ (unsigned long long) ((tv2.tv_sec - tv1.tv_sec)*1000000ULL);

	std::cout<<"Finalizing all algorithms..."<<std::endl;
	for (ialgo = m_algorithms.begin(); ialgo != m_algorithms.end(); ialgo++)
		(*ialgo)->parentFinalize();
	
	double totAvTime = 0;
	std::cout<<"\n\n--------- Algorithm Average Execute Time (us) ---------"<<std::endl;
	for (ialgo = m_algorithms.begin(); ialgo != m_algorithms.end(); ialgo++) {
		std::cout<<(*ialgo)->name()<<"\t\t"<<(*ialgo)->avTime()<<std::endl;
		totAvTime += (*ialgo)->avTime();
	}

  std::cout<<"\n"<<std::endl;
	double scopeEquivRead = m_nEvents*2048*geometry()->nGlibs()*geometry()->nChannels()*8/1000000.;
	std::cout<<"\nReal time scanned (scope mode + skips): \t"<<2048*m_nEvents*16/1000000.<<" (ms)"<<std::endl;
	if (runMode()==1) std::cout<<"Fraction tree read: \t\t\t\t"<<((SafEventBuilder*)m_algorithms[0])->treePos()/(1.*((SafEventBuilder*)m_algorithms[0])->chain()->GetEntries())<<std::endl;
	std::cout<<"Total algorithm average (per event): \t\t"<<totAvTime<<" (ms)"<<std::endl;
	std::cout<<"Total execution time (per event): \t\t"<<totExTime/(1000.*m_nEvents)<<" (ms)"<<std::endl;
	std::cout<<"Total execution time: \t\t\t\t"<<totExTime/1000000<<" (s)"<<std::endl;
}


//_____________________________________________________________________________

void SafRunner::eventLoop() {
  // Event loop.
	for (unsigned int i=0; i<m_nEvents; i++) {
		bool eof = false;
  	if (m_event % m_printRate == 0) {
      if (runMode()==1) {
    		double frac = ((SafEventBuilder*)m_algorithms[0])->treePos()/(1.*((SafEventBuilder*)m_algorithms[0])->chain()->GetEntries());
    		std::cout<<"Event: "<<m_event<<"\tFrac read: "<<frac<<"\t"<<"(/"<<nEvents()<<" events or EOF) \t Current File:\t"<<((SafEventBuilder*)m_algorithms[0])->chain()->GetFile()->GetName()<<std::endl;
    	}
      else std::cout<<"Event: "<<m_event<<"\t/\t"<<nEvents()<<std::endl;
    }

		for (std::vector< SafAlgorithm* >::iterator ialgo = m_algorithms.begin(); 
			ialgo != m_algorithms.end(); ialgo++) {
		  (*ialgo)->parentExecute();
			if ((*ialgo)->eof()) {
        std::cout<<"EOF"<<std::endl;
        eof = true;
        break;
      }
	  }

		rawData()->clear();
		triggerData()->clear();

		if (eof) break;
	  m_event++;
	}
}


//_____________________________________________________________________________
