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
  	std::cout<<fileName<<std::endl;
  	runner.setFileName(fileName);
  }

	runner.safPrint("Running Saffron", 0);
	runner.run();
	return 0;
}
