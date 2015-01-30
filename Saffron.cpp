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
  for (unsigned int i=1; i<argc; i++) {
  	std::string arg(argv[i]);
  	runner.evalArg(arg);
  }

	std::cout<<"Running Saffron."<<std::endl;
	runner.run();
	return 0;
}




