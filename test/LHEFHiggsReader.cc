// -*- mode: C++ -*-
// to compile:
// c++ -o LHEFHiggsReader `root-config --glibs --cflags` LHEFHiggsReader.cc

#include "LHEF.h"
#include "TFile.h"
#include "TTree.h"
#include <iomanip>
#include <cmath>

int main(int argc, char ** argv) {

  if (argc < 3) {
    std::cout << "requires 2 arguments <input lhe filename> <output root filename>." << std::endl;
    return 1;
  }

  // Create Reader and Writer object
  LHEF::Reader reader(argv[1]);

  TFile output(argv[2], "recreate");
  TTree outTree("outTree", "outTree");

  double Higgs_mass(-999.), NonRunning_wgt(1.), Running_wgt(1.);

  outTree.Branch("Higgs_mass", &Higgs_mass, "Higgs_mass/D");
  outTree.Branch("NonRunning_wgt", &NonRunning_wgt, "NonRunning_wgt/D");
  outTree.Branch("Running_wgt", &Running_wgt, "Running_wgt/D");

  long neve = 0;
  int Higgs_i(-1);

  // Read each event and write them out again.
  while ( reader.readEvent() ) {
    ++neve;
    if (neve%10000 == 0) std::cout << "event: " << neve << std::endl;
    Higgs_mass = -999.;
    Higgs_i = -1;

    for (int parti = 0; parti < reader.hepeup.NUP; ++parti) {
      if (abs(reader.hepeup.IDUP[parti]) == 25) {
	Higgs_i = parti;
	Higgs_mass = reader.hepeup.PUP[parti][4];
	// std::cout << neve << " parti: " << parti
	// 	  << " id: " << reader.hepeup.IDUP[parti]
	// 	  << " mass: " << Higgs_mass
	// 	  << " status: " << reader.hepeup.ISTUP[parti]
	// 	  << "\n";
      }
    }
    outTree.Fill();

    // if (neve >= 100) break;
  }

  output.Write();
  output.Close();

  std::cerr << neve << " events were found." << std::endl;

}
