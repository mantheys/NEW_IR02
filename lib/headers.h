#include <fstream>
#include <iostream>
#include <numeric>

#include <TTimer.h>
#include <TH1D.h>
#include <TFile.h>
#include <TF1.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TKey.h>
#include <TGaxis.h>
#include <TProfile.h>
#include <TH2D.h>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <TText.h>
#include <TMultiGraph.h>
#include <TText.h>
#include <TNtuple.h>
#include <TSpectrum.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TFitResult.h>

using namespace std;
void lets_pause()
{
  TTimer *timer = new TTimer("gSystem->ProcessEvents();", 50, kFALSE);
  timer->TurnOn();
  timer->Reset();
  std::cout << "q/Q to quit, other to continuee: ";
  char kkey;
  std::cin.get(kkey);
  if (kkey == 'q' || kkey == 'Q')
    throw std::exception(); // std::exit(0); //gSystem->Exit(0); //
  timer->TurnOff();
  delete timer;
}
int lets_pause(unsigned int &counter)
{
  TTimer *timer = new TTimer("gSystem->ProcessEvents();", 50, kFALSE);
  timer->TurnOn();
  timer->Reset();
  std::cout << "q/Q to quit, r to go back to last item, p to print and keep going, other to not print and keep going.";
  char kkey;
  std::cin.get(kkey);
  if (kkey == 'q' || kkey == 'Q')
    throw std::exception(); // std::exit(0); //gSystem->Exit(0); //
  if (kkey == 'r' || kkey == 'R')
    counter = counter - 2;
  timer->TurnOff();
  delete timer;
  if (kkey == 'P' || kkey == 'p')
    return 1;
  if( kkey == 'C' || kkey == 'c') 
    return 2;
  else
    return 0;
}
#include "mylib.h"
#include "lib.h"
#include "Waveform.h"
#include "Event.h"
#include "Cuts.h"
#include "NoiseTools.h"
#include "Run.h"
#include "RunCollection.h"
#include "Decon.h"
#include "InputConfig.h"
