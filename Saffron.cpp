//============================================================================
// Name        : Saffron.cpp
// Author      : Daniel Saunders
// Description : Top level runner of Saffron
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include "src/SafRunner.h"
#include <TROOT.h>

int main(int argc, char *argv[])
{
  gROOT->ProcessLine("gErrorIgnoreLevel = kFatal;");
  SafRunner runner;


  if (argc == 2) {
  	std::string fileName(argv[1]);
  	runner.setFileName(fileName);
  }

  else if (argc == 3) {
  	std::string fileName(argv[1]);
  	runner.setFileName(fileName);
  	runner.triggerThreshold = atof(argv[2]);
  	std::string ext(argv[2]);
  	std::string saveFileName = "Saffron-histos" + ext + ".root";
  	runner.setSaveFileName(saveFileName);
  }

	runner.safPrint("Running Saffron", 0);
	runner.run();

	return 0;
}
