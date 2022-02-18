#ifndef Run_h
#define Run_h
#include <stdio.h>
#include "HistogramCollection.h"
#include "Decon.h"
#include <iostream>
#include <string>
#include <cstring>

namespace ana
{
   float Square(float value){
    // Multiply value two times
    return value*value;
   }
   
   class Run_t
   {
   public:
      int Number = 0;
      string FolderName;
      double range1, range2;
      int PedRange;
      int TriggerRange;
      string adcmode;

      int firsttime;
      int lasttime;

      unsigned int NEvents;
      unsigned int NEventsALL;
      int MaximumEventsToProcess;
      EventReader_t event;
      int NSamples;
      int NChannels;

      double Sampling = 4e-9;
      int ADCDynamicRange;

      std::vector<string> PMT_SN = {"SiPM1", "SiPM2", "PMT", "Supercells"};
      std::vector<double> PMT_Voltages;
      std::vector<double> PMT_Gains;
      std::vector<double> PMT_SPEAmp;

      std::vector<int> adcchannels;

      std::vector<Event_t> EventList;
      double RunDuration;

      std::vector<double> PedestalPeak, PedestalSTDCUT, PedestalSigma, PedestalSTD;
      Cuts_t MyCuts;

      string calibration = "Gains/gains_20190912.csv";

      waveana::WaveAnaParameters *ParSet = new waveana::WaveAnaParameters();

      Run_t() {}
      void getADC(int n)
      {
         if (n > 14 && n <= 310)
            adcmode = "v1720";
         else
            adcmode = "v1740";
         // https://pddpelog.web.cern.ch/elisa/display/119
      }

      void SanityCheck()
      {
         if (range1 > range2)
         {
            std::cout << "Error in ranges set up by the user. Range2 must be larger than Range 1, please check this!" << std::endl;
            throw std::exception();
         }
      }
      Run_t(int r, std::vector<std::pair<string, string>> ChannelsToRead, string adcm = "", double r1 = 1.e-6, double r2 = 2.e-6, int pr = 30, int mevtp = 200, string cc = "Gains/gains_20190912.csv", int tr = 50) : Number(r), adcmode(adcm), range1(r1), range2(r2), PedRange(pr), TriggerRange(tr), MaximumEventsToProcess(mevtp), calibration(cc)
      {

         cout << "Processing run " << Number << endl;
         SanityCheck();

         ParSet->setPedRange(PedRange);
         ParSet->setTriggerRange(TriggerRange);
         ParSet->setRange1(range1);
         ParSet->setRange2(range2);

         if (adcmode == "DT5725")
         {
            //ParSet->setConversionFactor(-(16384.0 / 2.0) * 50.0);
            ParSet->setConversionFactor(-(16384.0/2.0)*1030);
            ADCDynamicRange = 16384;
            ParSet->setBaselineAmpLimit(24);
            Sampling = 4e-9;
         }
         else if (adcmode == "v1720")
         {
            ParSet->setConversionFactor(-(4096.0 / 2.0) * 50.0);
            ADCDynamicRange = 4096;
            ParSet->setBaselineAmpLimit(6);
         }
         else if (adcmode == "Osc")
         {
            ParSet->setConversionFactor(-1.0);
            ADCDynamicRange = 1000;
            ParSet->setBaselineAmpLimit(6);
         }
         //      EventReader_t ev;
         std::vector<string> ch;
         std::vector<int> adcch;
         int i = 0;
         for (auto c : ChannelsToRead)
         {
            if (c.first == "")
            {
               i++;
               continue;
            }
            if (adcmode == "Osc")
               event.SetFile(c.first.c_str(), c.second.c_str(), true);
            else
               event.SetFile(c.first.c_str(), c.second.c_str());
            adcch.push_back(i);
            ch.push_back(c.second);
            i++;
         }
         if (adcmode == "Osc")
         {
            GetEntry(0);
            Sampling = event.GetSampling();
         }
         SelectChannels(adcch);
         SetPMTs(ch); // Manually set the PTMs of interest
         NChannels = ch.size();

         NEvents = event.GetNEvents();
         NEventsALL = event.GetNEvents();
         if (MaximumEventsToProcess == -1)
            MaximumEventsToProcess = NEventsALL;
         if (MaximumEventsToProcess < (int)NEvents)
            NEvents = MaximumEventsToProcess;

         event.GetEntry(0);
         firsttime = event.GetPCTimeStamp();
         cout << "First Event time" << firsttime << endl;
         event.GetEntry(NEventsALL - 1);
         lasttime = event.GetPCTimeStamp();
         cout << "Last Event time" << lasttime << endl;
         cout << "Duration " << lasttime - firsttime << endl;
         NSamples = event.GetNSamples();
         //      Process();
         std::cout << "Run Loaded: " << NEventsALL << " events. " << std::endl;
         std::cout << MaximumEventsToProcess << " events will be processed. " << std::endl;
         std::cout << NSamples << " samples, " << Sampling * 1.e9 << "ns sampling. " << std::endl;
      }
      void Process()
      {
         //      int scale=1;
         EventList.clear();
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            GetEntry(ev); // if(ev==1000) scale=1;
            if (ev % 1000 == 0)
               std::cout << "Processing event " << ev << " out of " << NEvents << std::endl;
            EventList.push_back(ProcessEvent(ev));
         }
         TH1F *haux, *thped;
         PedestalSTDCUT.resize(NChannels);
         PedestalPeak.resize(NChannels);
         PedestalSigma.resize(NChannels);
         PedestalSTD.resize(NChannels);
         for (auto pm : adcchannels)
         {
            haux = TH1Pedestal(pm);
            //         double rmin= haux->GetBinCenter(haux->GetMaximumBin())-2*haux->GetStdDev();
            //         double rmax= haux->GetBinCenter(haux->GetMaximumBin())+2*haux->GetStdDev();

            double rmin = haux->GetMean() - 2 * haux->GetStdDev();
            double rmax = haux->GetMean() + 2 * haux->GetStdDev();
            TF1 *f = new TF1(Form("f%i", pm), "gaus");
            f->SetRange(rmin, rmax);
            f->ReleaseParameter(0);
            f->ReleaseParameter(1);
            f->ReleaseParameter(2);
            haux->Draw();
            haux->Fit(Form("f%i", pm), "LMERQ", "SAME", rmin, rmax); // if(pm==44){cout << rmin << " " << rmax << endl; lets_pause();}

            PedestalPeak[pm] = f->GetParameter(1);
            PedestalSigma[pm] = f->GetParameter(2);
            //         if(PedestalPeak[pm]>rmax || PedestalPeak[pm]<rmin) {cout << "PedesetalPeak out of range!!! What is happening??? " << endl; lets_pause();}
            //         cout << " ch " << pm << " - PedPEAK on " << PedestalPeak[pm] << endl;
            // if(pm==20){haux->Draw(); cout << " range " << rmin << " " << rmax<< endl; TH2PedestalVsPedSTD(pm)->Draw("COLZ");lets_pause();}
            thped = new TH1F(Form("Run%i_OpChannel%i_PedSTD", Number, pm), Form("Run%i_OpChannel%i_PedSTD;Pedestal STD (ADC);Events", Number, pm), 1000, 0, 20);
            for (unsigned int i = 0; i < NEvents; i++)
            {
               thped->Fill(getWaveform(i, pm)->getPedestalSTD());
            }

            rmin = thped->GetMean() - 2 * thped->GetStdDev();
            rmax = thped->GetMean() + 2 * thped->GetStdDev();
            TF1 *f2 = new TF1(Form("f2%i", pm), "gaus");
            f2->SetRange(rmin, rmax);
            f2->ReleaseParameter(0);
            f2->ReleaseParameter(1);
            f2->ReleaseParameter(2);
            thped->Draw();
            thped->Fit(Form("f2%i", pm), "LMERQ", "SAME", rmin, rmax);
            /*
                     Double_t xq[1];  // position where to compute the quantiles in [0,1]
                     Double_t yq[1];  // array to contain the quantiles
                     xq[0] = 0.88;
                     thped->GetQuantiles(1,yq,xq);
                     PedestalSTDCUT[pm]=yq[0];
                     PedestalSTD[pm]=thped->GetMean();
            */
            PedestalSTDCUT[pm] = thped->GetMean() + 2 * f2->GetParameter(2);
            PedestalSTD[pm] = thped->GetMean();
            cout << " ch " << pm << " Baseline on " << PedestalPeak[pm] << " - Sigma: " << PedestalSigma[pm] << " - PedSTD: " << PedestalSTD[pm] << " - PedSTDCUT: " << PedestalSTDCUT[pm] << endl;
            delete thped;
            delete haux;
            delete f;
         }
         cout << "Process done. " << endl;
      }

      Event_t ProcessEvent(int ev)
      {
         GetEntry(ev);
         Event_t myevt;
         myevt.setTimeStamp(event.GetTimeStamp());
         for (auto pm : adcchannels)
         {
            waveana::Waveform_t mywvf(event.GetTime(), event.GetAmp(pm), ParSet);
            myevt.AddWaveform(mywvf, pm); // std::cout << EventList[i].getChannel(pm)->getPedestalMean() << endl;
         }
         return myevt;
      }
      waveana::Waveform_t ProcessWaveform(int ev, int pm)
      {
         GetEntry(ev);
         waveana::Waveform_t mywvf(event.GetTime(), event.GetAmp(pm), ParSet);
         return mywvf;
      }

      waveana::Waveform_t ReProcessWaveform(TH1D *h)
      {
         std::vector<double> x, y;
         for (int i = 1; i <= h->GetSize(); i++)
         {
            x.push_back(h->GetBinCenter(i));
            y.push_back(h->GetBinContent(i));
         }
         waveana::Waveform_t mywvf(&x, &y, ParSet);
         return mywvf;
      }

      waveana::Waveform_t *getWaveform(int ev, int pm)
      {
         if (EventList.size() == 0)
            Process();
         if (ev < MaximumEventsToProcess)
         {
            return EventList[ev].getChannel(pm);
         }
         else
         {
            std::cout << "Error, using waveana::Waveform_t* getWaveform(int ev, int pm) with ev>MaximumEventsToProcess!!! ev = " << ev << " " << MaximumEventsToProcess << std::endl;
            throw std::exception();
         }
      }

      Event_t *getEvent(int ev)
      {
         if (EventList.size() == 0)
            Process();
         if (ev <= MaximumEventsToProcess)
            return &(EventList[ev]);
         else
         {
            std::cout << "Error, you are trying to get an event in EventList with ev>MaximumEventsToProcess!!!" << ev << " " << MaximumEventsToProcess << std::endl;
            throw std::exception();
         }
      }

      TH1D *TH1getWaveform(int ev, int pm)
      {
         GetEntry(ev);
         return event.GetWaveform(pm);
      }
      double getPedestal(int pm)
      {
         if (EventList.size() == 0)
            Process();
         return PedestalPeak[pm];
      }
      double getPedestalSTD(int pm)
      {
         if (EventList.size() == 0)
            Process();
         return PedestalSTD[pm];
      }
      double getPedestalSigma(int pm)
      {
         if (EventList.size() == 0)
            Process();
         return PedestalSigma[pm];
      }

      std::vector<double> getPedestals()
      {
         if (EventList.size() == 0)
            Process();
         return PedestalPeak;
      }
      std::vector<double> getPedestalSTDs()
      {
         if (EventList.size() == 0)
            Process();
         return PedestalSTD;
      }
      std::vector<double> getPedestalSigmas()
      {
         if (EventList.size() == 0)
            Process();
         return PedestalSigma;
      }

      string getPMTChannel(int pm) { return PMT_SN[pm]; }
      double getPMT_Voltage(int pm) { return PMT_Voltages[pm]; }
      double getPMT_Gain(int pm) { return PMT_Gains[pm]; }
      double getSampling() { return Sampling; }

      std::vector<string> getPMTChannels() { return PMT_SN; }
      std::vector<double> getPMT_Voltages() { return PMT_Voltages; }
      std::vector<double> getPMT_Gains() { return PMT_Gains; }

      std::vector<int> getADCChannels() { return adcchannels; }
      int getRunNumber() { return Number; }

      void GetEntry(int ev) { event.GetEntry(ev); }

      /* -- Functions to set cuts -- */
      void ResetCuts()
      {
         cout << "RESET CUTS." << endl;
         if (EventList.size() == 0)
            Process();
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : adcchannels)
               getWaveform(ev, pm)->setCut(false);
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            EventList[ev].setCut(false);
         MyCuts.ResetCuts();
      }
      void SetCutPedestalStatus()
      {
         cout << "CUT: Removing any WAVEFORM with PedestalStatus FALSE." << endl;
         if (EventList.size() == 0)
            Process();
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : adcchannels)
            {
               if (!getWaveform(ev, pm)->getPedestalStatus())
               {
                  getWaveform(ev, pm)->setCut(true);
               }
            }
         MyCuts.SetCutPedestalStatus();
         CutStatus();
      }
      void SetCutPedSTD(double std = 0)
      {
         if (EventList.size() == 0)
            Process();
         if (std == 0)
         {
            cout << "CUT: Removing any WAVEFORM with PedestalRMS larger than PedestalSTDCUT (2xsigma from the PedestalSTD Mean)." << endl;
            for (unsigned int ev = 0; ev < EventList.size(); ev++)
               for (auto pm : adcchannels)
               {
                  if (getWaveform(ev, pm)->getPedestalSTD() > PedestalSTDCUT[pm])
                  {
                     getWaveform(ev, pm)->setCut(true);
                  }
               }
            MyCuts.SetCutPedSTD(PedestalSTDCUT);
         }
         else
         {
            cout << "CUT: Removing any WAVEFORM with PedestalRMS larger than " << std << "ADC." << endl;
            for (unsigned int ev = 0; ev < EventList.size(); ev++)
               for (auto pm : adcchannels)
               {
                  if (getWaveform(ev, pm)->getPedestalSTD() > std)
                  {
                     getWaveform(ev, pm)->setCut(true);
                  }
               }
            MyCuts.SetCutPedSTD(std::vector<double>(PedestalSTDCUT.size(), std));
         }
         CutStatus();
      }
      void SetCutBaselineStability(double sigma = 1.0)
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Removing any WAVEFORM with baseline outside range of BaselineMean +- " << sigma << "sigma. " << endl;
         if (sigma < 0)
         {
            cout << "Error, sigma must be positive: " << sigma << endl;
            throw std::exception();
         }
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : adcchannels)
            {
               if (getWaveform(ev, pm)->getPedestalMean() > PedestalPeak[pm] + sigma * PedestalSigma[pm] || getWaveform(ev, pm)->getPedestalMean() < PedestalPeak[pm] - sigma * PedestalSigma[pm])
               {
                  getWaveform(ev, pm)->setCut(true);
               }
            }
         MyCuts.SetCutBaselineStability(sigma, PedestalPeak, PedestalSigma);
         CutStatus();
      }
      void SetCutSaturated()
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Removing any WAVEFORM with saturated signals (at least one sample of value 1ADC or ADC) " << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : adcchannels)
               if (getWaveform(ev, pm)->getSaturation())
                  getWaveform(ev, pm)->setCut(true);
         MyCuts.SetSaturationCut();
         CutStatus();
      }
      void SetCutMaximumMaximumSample(double Sample)
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Selecting WAVEFORM with the Maximum sample equal or below " << Sample << ". " << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : adcchannels)
               if (getWaveform(ev, pm)->getMaximumSample() > Sample)
                  getWaveform(ev, pm)->setCut(true);
         //       MyCuts.SetCutMaximumMaximumSample(Sample,TriggerChannels);
         CutStatus();
      }
      void SetCutOvershootingRange(int Sample)
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Removing WAVEFORM with the Maximum sample within range is " << Sample << " ADC above baseline (overshooting). " << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : adcchannels)
               if (getWaveform(ev, pm)->getMaximumSample_Range() - getWaveform(ev, pm)->getPedestalMean() > Sample)
                  getWaveform(ev, pm)->setCut(true);
         MyCuts.SetCutOvershootingRange(Sample);
         CutStatus();
      }
      void SetCutNsamplesBelowThresholdRebin(int nsamples, int ADCThreshold, int Rebin)
      {
         if (EventList.size() == 0)
            Process();
         MyCuts.SetCutNsamplesBelowThresholdRebin(nsamples, ADCThreshold, Rebin);
         cout << "CUT: Selecting WAVEFORM with " << nsamples << " below " << ADCThreshold << "ADC after doing a rebin of " << Rebin << " bins. " << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : adcchannels)
               if (!MyCuts.CheckCutNsamplesBelowThresholdRebin(*(getWaveform(ev, pm)), TH1getWaveform(ev, pm)))
                  getWaveform(ev, pm)->setCut(true);
         CutStatus();
      }
      void SetCutMinimumMaximumSample(double Sample)
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Selecting WAVEFORM with the Maximum sample above " << Sample << ". " << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : adcchannels)
               if (getWaveform(ev, pm)->getMaximumSample() <= Sample)
                  getWaveform(ev, pm)->setCut(true);
         MyCuts.SetCutMinimumMaximumSample(Sample);
         CutStatus();
      }
      void SetCutMaxAmplitude(double minAmp, double maxAmp)
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Selecting any WAVEFORM with maximum amplitude between " << minAmp << "ADC and " << maxAmp << "ADC." << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : adcchannels)
               if (getWaveform(ev, pm)->getMaxAmplitude() < minAmp || getWaveform(ev, pm)->getMaxAmplitude() > maxAmp)
                  getWaveform(ev, pm)->setCut(true);
         MyCuts.SetCutMaxAmplitude(minAmp, maxAmp);
         CutStatus();
      }

      void SetCutMaxAmplitudeRange(double minAmp, double maxAmp)
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Selecting any WAVEFORM with maximum amplitude in the range (" << range1 << "s," << range2 << "s) - between " << minAmp << "ADC and " << maxAmp << "ADC." << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : adcchannels)
               if (getWaveform(ev, pm)->getMaxAmplitudeRange() < minAmp || getWaveform(ev, pm)->getMaxAmplitudeRange() > maxAmp)
                  getWaveform(ev, pm)->setCut(true);
         MyCuts.SetCutMaxAmplitudeRange(minAmp, maxAmp);
         CutStatus();
      }

      void SetCutRemoveOutRangePeaksAbove(double maxAmp)
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Removing WAVEFORMS with a peak larger than " << maxAmp << "ADC, outside the range " << range1 << "s and " << range2 << "s." << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : adcchannels)
               if (getWaveform(ev, pm)->getMaxAmplitudeOutOfRange() > maxAmp)
                  getWaveform(ev, pm)->setCut(true);
         MyCuts.SetCutRemoveOutRangePeaksAbove(maxAmp);
         CutStatus();
      }

      void SetCutRemovePeaksOutOfMaxPeakAbove(double maxAmp) // cut events with secondary signals above maxAmp
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Removing WAVEFORMS with a peak larger than " << maxAmp << "ADC, 20ns before or 100ns after the first peak found within the range " << range1 << "s and " << range2 << "s." << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : adcchannels)
               if (getWaveform(ev, pm)->getMaxAmplitudeOutOfMaxPeak() > maxAmp)
                  getWaveform(ev, pm)->setCut(true);
         MyCuts.SetCutRemovePeaksOutOfMaxPeakAbove(maxAmp);
         CutStatus();
      }
      /*

                          ,getWaveform(ev,pm)->getMaxAmplitudeRange()
                         ,getWaveform(ev,pm)->getMaxAmplitudeRange()*Sampling*2/ADCDynamicRange/50.0/1.602e-19/PMT_SPEAmp[pm]
                         ,getWaveform(ev,pm)->getQ2MaxPeakRange()*1e12
                         ,getWaveform(ev,pm)->getQ3MaxPeakRange()*1e12
                         ,getWaveform(ev,pm)->getQ2MaxPeakRange()*1.0/1.602e-19/PMT_Gains[pm]
                         ,getWaveform(ev,pm)->getQ3MaxPeakRange()*1.0/1.602e-19/PMT_Gains[pm] );
      */

      void SetCutVariable(string var, double minAmp, double maxAmp)
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Selecting any WAVEFORM with variable " << var << " - between " << minAmp << " and " << maxAmp << " ." << endl;

         if (var == "MaxAmplitudeRange")
         {
            for (unsigned int ev = 0; ev < EventList.size(); ev++)
               for (auto pm : adcchannels)
                  if (getWaveform(ev, pm)->getMaxAmplitudeRange() < minAmp || getWaveform(ev, pm)->getMaxAmplitudeRange() > maxAmp)
                     getWaveform(ev, pm)->setCut(true);
         }
         else if (var == "Q2MaxPeakRange")
         {
            for (unsigned int ev = 0; ev < EventList.size(); ev++)
               for (auto pm : adcchannels)
                  if (getWaveform(ev, pm)->getQ2MaxPeakRange() < minAmp || getWaveform(ev, pm)->getQ2MaxPeakRange() > maxAmp)
                     getWaveform(ev, pm)->setCut(true);
         }
         else if (var == "Q3MaxPeakRange")
         {
            for (unsigned int ev = 0; ev < EventList.size(); ev++)
               for (auto pm : adcchannels)
                  if (getWaveform(ev, pm)->getQ3MaxPeakRange() < minAmp || getWaveform(ev, pm)->getQ3MaxPeakRange() > maxAmp)
                     getWaveform(ev, pm)->setCut(true);
         }
         else if (var == "PreTriggerSTD")
         {
            for (unsigned int ev = 0; ev < EventList.size(); ev++)
               for (auto pm : adcchannels)
                  if (getWaveform(ev, pm)->getPreTriggerSTD() < minAmp || getWaveform(ev, pm)->getPreTriggerSTD() > maxAmp)
                     getWaveform(ev, pm)->setCut(true);
         }
         else if (var == "ChargeMaxPeakRange")
         {
            for (unsigned int ev = 0; ev < EventList.size(); ev++)
               for (auto pm : adcchannels)
                  if (getWaveform(ev, pm)->getChargeMaxPeakRange() < minAmp || getWaveform(ev, pm)->getChargeMaxPeakRange() > maxAmp)
                     getWaveform(ev, pm)->setCut(true);
         }
         else if (var == "TEndMaxPeakRange")
         {
            for (unsigned int ev = 0; ev < EventList.size(); ev++)
               for (auto pm : adcchannels)
                  if (getWaveform(ev, pm)->getTEndMaxPeakRange() < minAmp || getWaveform(ev, pm)->getTEndMaxPeakRange() > maxAmp)
                     getWaveform(ev, pm)->setCut(true);
         }
         else
         {
            cout << "Error in cut, var not defined" << endl;
            throw std::exception();
         }
         MyCuts.SetCutVariable(var, minAmp, maxAmp);
         CutStatus();
      }
      void SetCutVariableVector(string var, std::map<int, std::pair<double, double>> map)
      {
         if (EventList.size() == 0)
            Process();
         for (auto pm : adcchannels)
            cout << "CUT: Selecting any WAVEFORM with variable " << var << " - between " << map[pm].first << " and " << map[pm].second << " ." << endl;
         if (var == "PreTriggerSTD")
         {
            for (unsigned int ev = 0; ev < EventList.size(); ev++)
               for (auto pm : adcchannels)
                  if (getWaveform(ev, pm)->getPreTriggerSTD() < map[pm].first || getWaveform(ev, pm)->getPreTriggerSTD() > map[pm].second)
                     getWaveform(ev, pm)->setCut(true);
         }
         else if (var == "MaxAmplitudeRange")
         {
            for (unsigned int ev = 0; ev < EventList.size(); ev++)
               for (auto pm : adcchannels)
                  if (getWaveform(ev, pm)->getMaxAmplitudeRange() < map[pm].first || getWaveform(ev, pm)->getMaxAmplitudeRange() > map[pm].second)
                     getWaveform(ev, pm)->setCut(true);
         }
         else if (var == "ChargeMaxPeakRange")
         {
            for (unsigned int ev = 0; ev < EventList.size(); ev++)
               for (auto pm : adcchannels)
                  if (getWaveform(ev, pm)->getChargeMaxPeakRange() < map[pm].first || getWaveform(ev, pm)->getChargeMaxPeakRange() > map[pm].second)
                     getWaveform(ev, pm)->setCut(true);
         }
         else if (var == "TEndMaxPeakRange")
         {
            for (unsigned int ev = 0; ev < EventList.size(); ev++)
               for (auto pm : adcchannels)
                  if (getWaveform(ev, pm)->getTEndMaxPeakRange() < map[pm].first || getWaveform(ev, pm)->getTEndMaxPeakRange() > map[pm].second)
                     getWaveform(ev, pm)->setCut(true);
         }
         else
         {
            cout << "Error in cut, var not defined" << endl;
            throw std::exception();
         }
         MyCuts.SetCutVariableVector(var, map);
         CutStatus();
      }
      void SetCutPeakTimeRange(double minTime = 0, double maxTime = 0)
      {
         if (minTime == 0)
         {
            minTime = range1;
            maxTime = range2;
         }
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Selecting WAVEFORMS with peaktime between " << minTime << "s and " << maxTime << "s on range (" << range1 << "s," << range2 << "s)" << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : adcchannels)
               if (getWaveform(ev, pm)->getPeakTime() < minTime || getWaveform(ev, pm)->getPeakTime() > maxTime)
                  getWaveform(ev, pm)->setCut(true);
         MyCuts.SetCutPeakTimeRange(minTime, maxTime);
         CutStatus();
      }
      void SetCutTriggerAmplitudeRange(double minAmp, double maxAmp, std::vector<int> TriggerChannels)
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Selecting EVENTS with a signal amplitude between " << minAmp << " and " << maxAmp << "ADC on range (" << range1 << "s," << range2 << "s) on ALL channels: ";
         for (auto pm : TriggerChannels)
            cout << pm << ", ";
         cout << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : TriggerChannels)
               if (getWaveform(ev, pm)->getMaxAmplitudeRange() < minAmp || getWaveform(ev, pm)->getMaxAmplitudeRange() > maxAmp)
               {
                  EventList[ev].setCut(true);
                  break;
               }
         MyCuts.SetCutTriggerAmplitudeRange(minAmp, maxAmp, TriggerChannels);
         CutStatus();
      }
      void SetCutTriggerFirstSampleBelowADCTriggerThresholdRange(double ADCThreshold, std::vector<int> TriggerChannels)
      {

         if (ParSet->ADCTriggerThreshold != ADCThreshold)
         {
            cout << ParSet->ADCTriggerThreshold << " != " << ADCThreshold << " -> Warning: ADCTriggerThreshold different in WaveAnaParameters. Reprocessing waveforms with the new configuration..." << endl;
            ParSet->setADCTriggerThreshold(ADCThreshold);
            Process();
         }
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Selecting EVENTS with an ADC sample below " << ADCThreshold << " on range (" << range1 << "s," << range2 << "s) on channels: ";
         for (auto pm : TriggerChannels)
            cout << pm << ", ";
         cout << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
         {
            unsigned int t = 0;
            for (auto pm : TriggerChannels)
               if (getWaveform(ev, pm)->getMinimumSampleRange() <= ADCThreshold)
                  t++;
            if (t == TriggerChannels.size())
               EventList[ev].setTriggerTime(getWaveform(ev, TriggerChannels[0])->getFirstSampleBelowADCTriggerThresholdRange());
            else
               EventList[ev].setCut(false);
         }
         MyCuts.SetCutTriggerFirstSampleBelowADCTriggerThresholdRange(ADCThreshold, TriggerChannels);
         CutStatus();
      }
      void SetCutTriggerAmplitude(double minAmp, double maxAmp, std::vector<int> TriggerChannels)
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Selecting EVENTS with a signal amplitude between " << minAmp << " and " << maxAmp << "ADC on ALL channels: ";
         for (auto pm : TriggerChannels)
            cout << pm << ", ";
         cout << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : TriggerChannels)
               if (getWaveform(ev, pm)->getMaxAmplitude() < minAmp || getWaveform(ev, pm)->getMaxAmplitude() > maxAmp)
               {
                  EventList[ev].setCut(true);
                  break;
               }
         //       MyCuts.SetCutTriggerAmplitude(minAmp,maxAmp,TriggerChannels);
         CutStatus();
      }
      void SetCutTriggerMinimumThreshold(double maxSample, std::vector<int> TriggerChannels)
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Selecting EVENTS with at least a sample equal or below " << maxSample << "ADC on channels: ";
         for (auto pm : TriggerChannels)
            cout << pm << ", ";
         cout << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : TriggerChannels)
               if (getWaveform(ev, pm)->getMinimumSample() > maxSample)
               {
                  EventList[ev].setCut(true);
                  break;
               }
         //       MyCuts.SetCutTriggerMinimumThreshold(minAmp,TriggerChannels);
         CutStatus();
      }
      void SetCutTriggerMaximumThreshold(double minSample, std::vector<int> TriggerChannels)
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Selecting EVENTS with at least a sample equal or below " << minSample << "ADC on channels: ";
         for (auto pm : TriggerChannels)
            cout << pm << ", ";
         cout << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : TriggerChannels)
               if (getWaveform(ev, pm)->getMinimumSample() > minSample)
               {
                  EventList[ev].setCut(true);
                  break;
               }
         //       MyCuts.SetCutTriggerMaximumThreshold(minSample,TriggerChannels);
         CutStatus();
      }
      void SetCutTriggerMinimumThresholdRange(double maxSample, std::vector<int> TriggerChannels)
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Selecting EVENTS with at least a sample equal or below " << maxSample << "ADC on range (" << range1 << "s," << range2 << "s) on channels: ";
         for (auto pm : TriggerChannels)
            cout << pm << ", ";
         cout << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : TriggerChannels)
               if (getWaveform(ev, pm)->getMinimumSampleRange() > maxSample)
               {
                  EventList[ev].setCut(true);
                  break;
               }
         MyCuts.SetCutTriggerMinimumThresholdRange(maxSample, TriggerChannels);
         CutStatus();
      }
      void SetCutTriggerMaximumThresholdRange(double minSample, std::vector<int> TriggerChannels)
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Removing EVENTS with at least a sample equal or below " << minSample << "ADC on range (" << range1 << "s," << range2 << "s) on channels: ";
         for (auto pm : TriggerChannels)
            cout << pm << ", ";
         cout << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : TriggerChannels)
               if (getWaveform(ev, pm)->getMinimumSampleRange() < minSample)
               {
                  EventList[ev].setCut(true);
                  break;
               }
         //       MyCuts.SetCutTriggerMinimumThresholdRange(minSample,TriggerChannels);
         CutStatus();
      }
      void SetCutTriggerSaturated(std::vector<int> TriggerChannels)
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Removing all EVENTS with saturated signals (at least one sample of value 1ADC or ADC) on channels: ";
         for (auto pm : TriggerChannels)
            cout << pm << ", ";
         cout << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : TriggerChannels)
               if (getWaveform(ev, pm)->getSaturation())
               {
                  EventList[ev].setCut(true);
                  break;
               }
         //       MyCuts.SetCutTriggerSaturated(minAmp,TriggerChannels);
         CutStatus();
      }

      void SetCutTriggerWaveformCuts(std::vector<int> TriggerChannels)
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Removing all EVENTS that have been cut from any of the channels: ";
         for (auto pm : TriggerChannels)
            cout << pm << ", ";
         cout << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : TriggerChannels)
               if (getWaveform(ev, pm)->getCut())
               {
                  EventList[ev].setCut(true);
                  break;
               }
         MyCuts.SetCutTriggerWaveformCuts(TriggerChannels);
         CutStatus();
      }

      void SetCutTriggerQuantile(string option, double percent, bool top, std::vector<int> TriggerChannels)
      {
         if (EventList.size() == 0)
            Process();
         if (top)
         {
            cout << "CUT: Removing " << percent * 100 << "/100 of the events on the end of the " << option << " distribution on channels: ";
            for (auto pm : TriggerChannels)
               cout << pm << ", ";
            cout << endl;
         }
         else
         {
            cout << "CUT: Removing " << percent * 100 << "/100 of the events on the start of the " << option << " distribution on channels: ";
            for (auto pm : TriggerChannels)
               cout << pm << ", ";
            cout << endl;
         }
         std::map<int, double> Quant;
         for (auto pm : TriggerChannels)
         {
            TH1D *thaux = new TH1D("aux", "aux", 10000, 0, 0);
            for (unsigned int ev = 0; ev < EventList.size(); ev++)
            {
               if (getEvent(ev)->getCut())
                  continue;
               if (getWaveform(ev, pm)->getCut())
                  continue;
               if (option == "AmpRange")
                  thaux->Fill(getWaveform(ev, pm)->getMaxAmplitudeRange());
               if (option == "ChargeRange")
                  thaux->Fill(getWaveform(ev, pm)->getChargeRange());
               if (option == "ChargeMaxPeak")
                  thaux->Fill(getWaveform(ev, pm)->getChargeMaxPeak());
               if (option == "MaxPeakRange")
                  thaux->Fill(getWaveform(ev, pm)->getChargeMaxPeakRange());
               if (option == "ChargeQ1MaxPeakRange")
                  thaux->Fill(getWaveform(ev, pm)->getQ1MaxPeakRange());
               if (option == "ChargeQ2MaxPeakRange")
                  thaux->Fill(getWaveform(ev, pm)->getQ2MaxPeakRange());
               if (option == "ChargeQ3MaxPeakRange")
                  thaux->Fill(getWaveform(ev, pm)->getQ3MaxPeakRange());
            }

            Double_t xq2[1]; // position where to compute the quantiles in [0,1]
            Double_t yq2[1]; // array to contain the quantiles
            if (!top)
               xq2[0] = percent;
            else
               xq2[0] = 1 - percent;
            thaux->GetQuantiles(2, yq2, xq2);
            Quant[pm] = yq2[0];
            cout << "Cut threshold stablished on " << Quant[pm] << "for Trigger channel " << pm << endl;
            thaux->Draw("HIST");
            lets_pause();
            delete thaux;
         }
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : TriggerChannels)
            {
               if (top)
               {
                  if (option == "AmpRange")
                     if (getWaveform(ev, pm)->getMaxAmplitudeRange() > Quant[pm])
                     {
                        EventList[ev].setCut(true);
                        break;
                     }
                  if (option == "ChargeRange")
                     if (getWaveform(ev, pm)->getChargeRange() > Quant[pm])
                     {
                        EventList[ev].setCut(true);
                        break;
                     }
                  if (option == "ChargeMaxPeak")
                     if (getWaveform(ev, pm)->getChargeMaxPeak() > Quant[pm])
                     {
                        EventList[ev].setCut(true);
                        break;
                     }
                  if (option == "ChargeMaxPeakRange")
                     if (getWaveform(ev, pm)->getChargeMaxPeakRange() > Quant[pm])
                     {
                        EventList[ev].setCut(true);
                        break;
                     }
                  if (option == "ChargeQ1MaxPeakRange")
                     if (getWaveform(ev, pm)->getQ1MaxPeakRange() > Quant[pm])
                     {
                        EventList[ev].setCut(true);
                        break;
                     }
                  if (option == "ChargeQ2MaxPeakRange")
                     if (getWaveform(ev, pm)->getQ2MaxPeakRange() > Quant[pm])
                     {
                        EventList[ev].setCut(true);
                        break;
                     }
                  if (option == "ChargeQ3MaxPeakRange")
                     if (getWaveform(ev, pm)->getQ3MaxPeakRange() > Quant[pm])
                     {
                        EventList[ev].setCut(true);
                        break;
                     }
               }
               else
               {
                  if (option == "AmpRange")
                     if (getWaveform(ev, pm)->getMaxAmplitudeRange() < Quant[pm])
                     {
                        EventList[ev].setCut(true);
                        break;
                     }
                  if (option == "ChargeRange")
                     if (getWaveform(ev, pm)->getChargeRange() < Quant[pm])
                     {
                        EventList[ev].setCut(true);
                        break;
                     }
                  if (option == "ChargeMaxPeak")
                     if (getWaveform(ev, pm)->getChargeMaxPeak() < Quant[pm])
                     {
                        EventList[ev].setCut(true);
                        break;
                     }
                  if (option == "ChargeMaxPeakRange")
                     if (getWaveform(ev, pm)->getChargeMaxPeakRange() < Quant[pm])
                     {
                        EventList[ev].setCut(true);
                        break;
                     }
                  if (option == "ChargeQ1MaxPeakRange")
                     if (getWaveform(ev, pm)->getQ1MaxPeakRange() < Quant[pm])
                     {
                        EventList[ev].setCut(true);
                        break;
                     }
                  if (option == "ChargeQ2MaxPeakRange")
                     if (getWaveform(ev, pm)->getQ2MaxPeakRange() < Quant[pm])
                     {
                        EventList[ev].setCut(true);
                        break;
                     }
                  if (option == "ChargeQ3MaxPeakRange")
                     if (getWaveform(ev, pm)->getQ3MaxPeakRange() < Quant[pm])
                     {
                        EventList[ev].setCut(true);
                        break;
                     }
               }
            }
         //       MyCuts.SetCutTriggerSaturated(minAmp,TriggerChannels);
         CutStatus();
      }

      void SetCutThresholdRange(double maxSample)
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Selecting WAVEFORMS with at least a sample equal or below " << maxSample << "ADC on range (" << range1 << "s," << range2 << "s)." << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : adcchannels)
               if (getWaveform(ev, pm)->getMinimumSampleRange() > maxSample)
               {
                  getWaveform(ev, pm)->setCut(true);
               }
         MyCuts.SetCutThresholdRange(maxSample);
         CutStatus();
      }

      void SetCutThreshold(double maxSample)
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Selecting any WAVEFORM with at least a sample equal or below " << maxSample << "ADC." << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : adcchannels)
               if (getWaveform(ev, pm)->getMinimumSample() > maxSample)
               {
                  getWaveform(ev, pm)->setCut(true);
               }
         MyCuts.SetCutThreshold(maxSample);
         CutStatus();
      }
      void SetCutS2Amplitude(double minAmp)
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Selecting any WAVEFORM with S2 Amplitude larger than " << minAmp << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : adcchannels)
               if (getWaveform(ev, pm)->getS2AverageAmplitude() < minAmp)
                  getWaveform(ev, pm)->setCut(true);
         MyCuts.SetCutS2Amplitude(minAmp);
         CutStatus();
      }
      void SetCutCharge(double minCharge, double maxCharge)
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Selecting any WAVEFORM with total charge between " << minCharge << "Coulombs and " << maxCharge << "Coulombs." << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : adcchannels)
               if (getWaveform(ev, pm)->getCharge() < minCharge || getWaveform(ev, pm)->getCharge() > maxCharge)
                  getWaveform(ev, pm)->setCut(true);
         MyCuts.SetCutCharge(minCharge, maxCharge);
         CutStatus();
      }
      void SetCutChargeRange(double minCharge, double maxCharge)
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Selecting any WAVEFORM with integrated charge in the range (" << range1 << "s," << range2 << "s) - between " << minCharge << "Coulombs and " << maxCharge << "Coulombs." << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : adcchannels)
               if (getWaveform(ev, pm)->getChargeRange() < minCharge || getWaveform(ev, pm)->getChargeRange() > maxCharge)
                  getWaveform(ev, pm)->setCut(true);
         MyCuts.SetCutChargeRange(minCharge, maxCharge);
         CutStatus();
      }
      void SetCutCustomChargeRange(double minCharge, double maxCharge)
      {
         if (EventList.size() == 0)
            Process();
         cout << "CUT: Selecting any WAVEFORM with total charge between " << minCharge << "Coulombs and " << maxCharge << "Coulombs. Charge integrated over a constant baseline." << endl;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
            for (auto pm : adcchannels)
            {
               double totalcharge = getWaveform(ev, pm)->getSampling() * (getWaveform(ev, pm)->getNSamples() * PedestalPeak[pm] - getWaveform(ev, pm)->getSumOfSamples()) * 2 / 4096 / 50;
               if (totalcharge < minCharge || totalcharge > maxCharge)
                  getWaveform(ev, pm)->setCut(true);
            }
         // MyCuts.SetCutCustomChargeRange(minCharge,maxCharge);
         CutStatus();
      }
      void CutStatus()
      {
         if (EventList.size() == 0)
            Process();
         cout << "Cut Status: " << endl;
         std::vector<int> counter(64, 0);
         int totcounter = 0;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
         {
            if (EventList[ev].getCut())
               continue;
            for (auto pm : adcchannels)
               if (!getWaveform(ev, pm)->getCut())
                  counter[pm]++;
            totcounter++;
         }
         for (auto pm : adcchannels)
            cout << "\tchannel " << pm << " - " << PMT_SN[pm] << ": " << counter[pm] << " events after cuts out of " << totcounter << ". " << endl;
         cout << "\tMultiChannelCut Status: " << totcounter << " events after cuts out of " << EventList.size() << ". " << endl;
      }

      int getWaveformsAfterCuts(int pm)
      {
         if (EventList.size() == 0)
            Process();
         int counter = 0;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
         {
            if (EventList[ev].getCut())
               continue;
            if (!getWaveform(ev, pm)->getCut())
               counter++;
         }
         return counter;
      }

      int getEventsAfterCuts()
      {
         if (EventList.size() == 0)
            Process();
         int totcounter = 0;
         for (unsigned int ev = 0; ev < EventList.size(); ev++)
         {
            if (EventList[ev].getCut())
               continue;
            totcounter++;
         }
         return totcounter;
      }

      /* -- Other Analysis Functions -- */

      TProfile *TProfile_ChargeVsTime(int pm)
      {
         TProfile *th = new TProfile(Form("Run%i_OpChannel%i_%s_%.0fV_TProfile_ChargeVsTime", Number, pm, PMT_SN[pm].c_str(), PMT_Voltages[pm]), ";Time;Charge (PE)", lasttime - firsttime, firsttime, lasttime);
         th->SetErrorOption("S");
         for (unsigned int ev = 0; ev < NEventsALL; ev++)
         {
            GetEntry(ev);
            th->Fill(event.GetTimeStamp(), waveana::Waveform_t(event.GetTime(), event.GetAmp(pm), ParSet).getCharge() / 1.602e-19 / PMT_Gains[pm]);
         }
         return th;
      }

      TProfile *TPVariableVsTime(string var, int pm)
      {
         TProfile *th = new TProfile(Form("Run%i_OpChannel%i_%s_%.0fV_TProf_%sVsTime", Number, pm, PMT_SN[pm].c_str(), PMT_Voltages[pm], var.c_str()), Form(";Time;%s", var.c_str()), lasttime - firsttime, firsttime, lasttime);
         th->SetErrorOption("S");
         if ("MaxAmplitudeRange")
            for (unsigned int ev = 0; ev < NEventsALL; ev++)
            {
               GetEntry(ev);
               th->Fill(event.GetPCTimeStamp(), waveana::Waveform_t(event.GetTime(), event.GetAmp(pm), ParSet).getMaxAmplitudeRange());
            }
         return th;
      }

      TH2D *TH2VariableVsTime(string var, int pm)
      {
         TH2D *th = new TH2D(Form("Run%i_OpChannel%i_%s_%.0fV_TH2_%sVsTime", Number, pm, PMT_SN[pm].c_str(), PMT_Voltages[pm], var.c_str()), Form(";Time;%s", var.c_str()), lasttime - firsttime, firsttime, lasttime, ADCDynamicRange, 0, ADCDynamicRange);
         if ("MaxAmplitudeRange")
            for (unsigned int ev = 0; ev < NEventsALL; ev++)
            {
               waveana::Waveform_t wv(event.GetTime(), event.GetAmp(pm), ParSet);
               GetEntry(ev);
               th->Fill(event.GetPCTimeStamp(), wv.getMaxAmplitudeRange());
               cout << "HEY, loot at event " << event.GetPCTimeStamp() << " " << event.GetTimeStamp() << " " << wv.getMaxAmplitudeRange() << endl;
               if (wv.getMaxAmplitudeRange() < 700)
               {
                  cout << "HEY, loot at event " << event.GetPCTimeStamp() << " " << event.GetTimeStamp() << " " << wv.getMaxAmplitudeRange() << endl;
               }
               if (ev % 1000 == 0)
               {
                  th->Draw("COLZ");
                  lets_pause();
               }
            }
         return th;
      }

      TH2D *TH2VariableVsEvent(string var, int pm)
      {
         TH2D *th = new TH2D(Form("Run%i_OpChannel%i_%s_%.0fV_TH2_%sVsTime", Number, pm, PMT_SN[pm].c_str(), PMT_Voltages[pm], var.c_str()), Form(";Time;%s", var.c_str()), 0.05 * NEventsALL, 0, NEventsALL, ADCDynamicRange, 0, ADCDynamicRange);
         if ("MaxAmplitudeRange")
            for (unsigned int ev = 0; ev < NEventsALL; ev++)
            {
               waveana::Waveform_t wv(event.GetTime(), event.GetAmp(pm), ParSet);
               GetEntry(ev);
               th->Fill(ev, wv.getMaxAmplitudeRange());
            }
         return th;
      }

      TH1D *TH1RateOfEventsPerPMT(int pm) // rate of events for this samples.
      {
         TH1D *th = new TH1D(Form("Run%i_OpChannel%i_%s_%.0fV_RateOfEventsPerPMT", Number, pm, PMT_SN[pm].c_str(), PMT_Voltages[pm]), ";Rate of Events(Hz);Counts", 10000, 0, 1000);
         th->Fill(NEvents / (lasttime - firsttime));
         return th;
      }

      int GetNEvents() { return NEventsALL; }
      int GetEventsToProcess() { return MaximumEventsToProcess; }

      double getAverage(int pm, string var = "Q1MaxPeakRange")
      {
         double Av = 0;
         int counter = 0;
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            if (EventList[ev].getCut())
               continue;
            if (getWaveform(ev, pm)->getCut())
               continue;
            if (var == "Q1MaxPeakRange")
               Av += getWaveform(ev, pm)->getQ1MaxPeakRange();
            if (var == "Q2MaxPeakRange")
               Av += getWaveform(ev, pm)->getQ2MaxPeakRange();
            if (var == "Q3MaxPeakRange")
               Av += getWaveform(ev, pm)->getQ3MaxPeakRange();
            if (var == "MaxAmplitudeRange")
               Av += getWaveform(ev, pm)->getMaxAmplitudeRange();
            if (var == "MinimumSampleRange")
               Av += getWaveform(ev, pm)->getMinimumSampleRange();
            counter++; // cout << getWaveform(ev,pm)->getMaxAmplitudeRange() << endl;
         }
         return Av / counter;
      }

      double getMedian(int pm, string var = "Q1MaxPeakRange")
      {
         TH1D *haux = new TH1D("haux", "haux", 10000, 0, 0);
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            if (EventList[ev].getCut())
               continue;
            if (getWaveform(ev, pm)->getCut())
               continue;
            if (var == "Q1MaxPeakRange")
               haux->Fill(getWaveform(ev, pm)->getQ1MaxPeakRange());
            if (var == "Q2MaxPeakRange")
               haux->Fill(getWaveform(ev, pm)->getQ2MaxPeakRange());
            if (var == "Q3MaxPeakRange")
               haux->Fill(getWaveform(ev, pm)->getQ3MaxPeakRange());
            if (var == "MaxAmplitudeRange")
               haux->Fill(getWaveform(ev, pm)->getMaxAmplitudeRange());
            if (var == "MinimumSampleRange")
               haux->Fill(getWaveform(ev, pm)->getMinimumSampleRange());
         }
         Double_t xq[1]; // position where to compute the quantiles in [0,1]
         Double_t yq[1]; // array to contain the quantiles
         xq[0] = 0.5;
         haux->GetQuantiles(1, yq, xq);
         haux->Delete();
         return yq[0];
      }

      int getRunDuration()
      {
         return lasttime - firsttime;
      }
      int getRunDurationNEvents(int nevt)
      {
         event.GetEntry(nevt - 1);
         return event.GetPCTimeStamp() - firsttime;
      }
      int getFirstTime()
      {
         return firsttime;
      }
      int getLastTime()
      {
         return lasttime;
      }
      TH1D *TH1RateOfEventsCut(int pm)
      {
         int bins = ((lasttime - firsttime) / 60);
         TH1D *th = new TH1D(Form("Run%i_OpChannel%i_%s_%.0fV_TH1RateOfEventsCut", Number, pm, PMT_SN[pm].c_str(), PMT_Voltages[pm]), "Rate after cuts;Time; Rate(Hz)", bins, firsttime, lasttime);
         // th->SetErrorOption("S");
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            //         waveana::Waveform_t mywf(event.GetTime(),event.GetAmp(pm),PedRange,TriggerRange,-(4096.0/2.0)*50.0, range1, range2);
            //         MyCuts.ApplyCuts(mywf,pm);
            GetEntry(ev);
            if (EventList[ev].getCut())
               continue;
            if (getWaveform(ev, pm)->getCut())
               continue;
            th->Fill(event.GetTimeStamp());
         }
         th->Scale(1.0 / 60.0);
         return th;
      }

      TH1D *TH1RateOfEvents() // evolution across time
      {
         int bins = ((lasttime - firsttime) / 60);
         TH1D *th = new TH1D(Form("Run%i_TH1RateOfEvents", Number), "Rate of events;Time;Rate (Hz)", bins, firsttime, lasttime);
         // th->SetErrorOption("S");
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            GetEntry(ev);
            th->Fill(event.GetTimeStamp());
         }
         th->Scale(1.0 / 60.0);
         return th;
      }
      TH1D *TH1TimeStamps() // evolution across time
      {
         GetEntry(0);
         long double FirstTS = event.GetTimeStamp();
         TH1D *th = new TH1D(Form("Run%i_TimeStamps", Number), ";Event;TimeStamp", NEvents, 0, NEvents);
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            GetEntry(ev); // cout << event.GetTimeStamp() << " " << FirstTS << " " << event.GetTimeStamp()-FirstTS << endl;
            th->SetBinContent(ev + 1, (event.GetTimeStamp() - FirstTS) * 8.e-9);
            //         th->SetBinContent(ev+1,(event.GetTimeStamp()));
         }
         return th;
      }

      TH1D *TH1PCTimeStamps() // evolution across time
      {
         GetEntry(0);
         long double FirstTS = event.GetPCTimeStamp();
         TH1D *th = new TH1D(Form("Run%i_PCTimeStamps", Number), ";Event;TimeStamp", NEvents, 0, NEvents);
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            GetEntry(ev);
            th->SetBinContent(ev + 1, event.GetPCTimeStamp() - FirstTS);
         }
         return th;
      }

      TH1D *TH1PCTimeStampsCorrected(double &length) // evolution across time
      {
         GetEntry(0);
         long double Shift = 0;
         long double OldTS = event.GetTimeStamp();
         long double NewTS;
         TH1D *th = new TH1D(Form("Run%i_PCTimeStampsCorrected", Number), ";Event;TimeStamp", NEvents, 0, NEvents);
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            GetEntry(ev);
            NewTS = event.GetTimeStamp();
            if (NewTS < OldTS)
               Shift += 2147483647;
            th->SetBinContent(ev + 1, 8.e-9 * (event.GetTimeStamp() + Shift));
            OldTS = NewTS;
         }
         cout << "Corrected run duration " << (th->GetBinContent(2) - th->GetBinContent(th->GetSize() - 2)) / 60.0 << " minutes." << endl;
         length = (th->GetBinContent(1) - th->GetBinContent(th->GetSize() - 2)) / 60.0;
         return th;
      }

      std::vector<double> myTimes;
      void FillMyTimes()
      {
         myTimes.resize(NEvents);
         GetEntry(0);
         long double Shift = 0;
         long double OldTS = event.GetTimeStamp();
         long double NewTS;
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            GetEntry(ev);
            NewTS = event.GetTimeStamp();
            if (NewTS < OldTS)
               Shift += 2147483647;
            myTimes[ev] = 8.e-9 * (event.GetTimeStamp() + Shift);
            OldTS = NewTS;
         }
      }
      TH2D *TH2DeltaTimes(int pm) // evolution across time
      {

         FillMyTimes();
         TH2D *th = new TH2D(Form("Run%i_DeltaTimes", Number), ";#Deltat (s);Amplitude (ADC)", 10000, 1.e-8, 4, 600, 0, 600);
         int lastev = 0;
         for (unsigned int ev = 1; ev < NEvents; ev++)
         {
            if (EventList[ev].getCut())
               continue;
            if (getWaveform(ev, pm)->getCut())
               continue;
            th->Fill(myTimes[ev] - myTimes[lastev], getWaveform(ev, pm)->getMaxAmplitudeRange());
            lastev = ev;
         }
         return th;
      }
      TProfile *TProfileEvolution(int pm, string variable, int bins)
      {
         double ChargeFactor = 1.e12;
         if (variable == "Charge_Q2MaxPeakRange_PE" || variable == "Charge_Q3MaxPeakRange_PE")
            ChargeFactor = 1.0 / 1.602e-19 / PMT_Gains[pm];

         cout << " TProfileEvolution Set bins " << bins << endl;
         TProfile *th = new TProfile(Form("Run%i_OpChannel%i_%s_%.0fV_Evolution_%s", Number, pm, PMT_SN[pm].c_str(), PMT_Voltages[pm], variable.c_str()), Form("%s;Time; Rate(Hz)", variable.c_str()), bins, firsttime, lasttime);
         //      th->SetErrorOption("S");
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            GetEntry(ev);
            if (EventList[ev].getCut())
               continue;
            if (getWaveform(ev, pm)->getCut())
               continue;
            if (variable == "MaxAmplitude")
               th->Fill(event.GetTimeStamp(), getWaveform(ev, pm)->getMaxAmplitude());
            if (variable == "pedestal")
               th->Fill(event.GetTimeStamp(), getWaveform(ev, pm)->getPedestalMean());
            if (variable == "MaxAmplitudeRange")
               th->Fill(event.GetTimeStamp(), getWaveform(ev, pm)->getMaxAmplitudeRange());
            if (variable == "ChargeRange")
               th->Fill(event.GetTimeStamp(), getWaveform(ev, pm)->getChargeRange() * ChargeFactor);
            if (variable == "ChargeRange_PE")
               th->Fill(event.GetTimeStamp(), getWaveform(ev, pm)->getChargeRange() * ChargeFactor);
            if (variable == "Charge_Q2MaxPeakRange")
               th->Fill(event.GetTimeStamp(), getWaveform(ev, pm)->getQ3MaxPeakRange() * ChargeFactor);
            if (variable == "Charge_Q2MaxPeakRange_PE")
               th->Fill(event.GetTimeStamp(), getWaveform(ev, pm)->getQ3MaxPeakRange() * ChargeFactor);
            if (variable == "Charge_Q3MaxPeakRange")
               th->Fill(event.GetTimeStamp(), getWaveform(ev, pm)->getQ3MaxPeakRange() * ChargeFactor);
            if (variable == "Charge_Q3MaxPeakRange_PE")
               th->Fill(event.GetTimeStamp(), getWaveform(ev, pm)->getQ3MaxPeakRange() * ChargeFactor);
         }
         return th;
      }

      TH1D *TH1Saturated(int pm)
      {
         TH1D *th = new TH1D(Form("Run%i_OpChannel%i_%s_%.0fV_Saturated", Number, pm, PMT_SN[pm].c_str(), PMT_Voltages[pm]), ";Rate of Events(Hz);Counts", 2, -0.5, 1.5);
         for (unsigned int ev = 1; ev < NEvents; ev++)
         {
            if (EventList[ev].getCut())
               continue;
            if (getWaveform(ev, pm)->getCut())
               continue;
            if (getWaveform(ev, pm)->getSaturation())
               th->Fill(1);
            else
               th->Fill(0);
         }
         return th;
      }

      TH1D *TH1AvWf(int pm)
      {
         cout << "Computing AverageWaveform for channel " << pm << ", just summing waveform. Applying waveform and trigger cuts, limited to NEvents set by user " << NEvents << ". " << endl;

         GetEntry(0);
         TH1D *th = event.GetWaveform(pm);
         th->Reset(); // th->Draw("HIST");
         th->SetName(Form("Run%i_OpChannel%i_AvWf", Number, pm));
         // th->SetTitle(Form("Run%i_OpChannel%i_%s_%.0fV_",Number,pm,PMT_SN[pm].c_str(),PMT_Voltages[pm]));
         int counter = 0;
         for (unsigned int ev = 1; ev < NEvents; ev++)
         {
            if (EventList[ev].getCut())
               continue;
            if (getWaveform(ev, pm)->getCut())
               continue;
            GetEntry(ev);
            TH1D *aux = event.GetWaveform(pm);
            th->Add(aux);
            counter++;
            delete aux;
         }
         th->Scale(1.0 / counter);
         cout << " Waveform averaged over channel " << pm << " using " << counter << " over all " << NEvents << " events. " << endl;
         //      th->Draw(); lets_pause();
         return th;
      }

      TH1D *TH1AvWfALL(int pm, int MaxNumberOfWaveforms = 200000)
      {
         cout << "Computing AverageWaveform for channel " << pm << ", just summing waveforms. Not applying any cut. " << endl;

         GetEntry(0);
         TH1D *th = event.GetWaveform(pm);
         th->Reset(); // th->Draw("HIST");
         th->SetName(Form("Run%i_OpChannel%i_AvWfALL", Number, pm));
         int counter = 0;
         if (MaxNumberOfWaveforms == -1)
            MaxNumberOfWaveforms = NEventsALL;
         for (unsigned int ev = 1; ev < NEventsALL && counter < MaxNumberOfWaveforms; ev++)
         {
            //         Waveform_t mywf=ProcessWaveform(ev,pm);
            //         MyCuts.ApplyCuts(mywf,pm);
            //         if(mywf.getCut()) continue;
            GetEntry(ev);
            TH1D *aux = event.GetWaveform(pm);
            th->Add(aux);
            counter++;
            delete aux;
            if (ev % 500 == 0)
               cout << "Added " << ev << " wvf." << endl;
         }
         th->Scale(1.0 / counter);
         cout << " Waveform averaged over channel " << pm << " using " << counter << " over all " << NEvents << " events. " << endl;
         return th;
      }

      TH2D *TH2AverageSignal(int pm, int MaxNumberOfWaveforms = 200000)
      {
         cout << "Computing TH2AverageSignal for channel " << pm << ", substracting the pedestal event by event, and shifting to the min bin within range. " << endl;
         cout << "Applying only waveform cuts, and trigger cuts" << endl;
         TH2D *th = new TH2D(Form("Run%i_OpChannel%i_AverageSignal", Number, pm), Form("Run%i_OpChannel%i_AverageSignal", Number, pm), NSamples, 0, NSamples * Sampling, ADCDynamicRange, 0, ADCDynamicRange);
         int counter = 0;
         int maxbin = (int)(range1 / Sampling);
         bool first = true;
         for (unsigned int ev = 0; ev < NEventsALL && counter < MaxNumberOfWaveforms; ev++)
         {
            Event_t evt = ProcessEvent(ev);
            MyCuts.ApplyMultiChannelCuts(evt);
            if (evt.getCut())
               continue;
            MyCuts.ApplyCuts(*evt.getChannel(pm), pm);
            if (evt.getChannel(pm)->getCut())
               continue;
            GetEntry(ev);
            TH1D *aux = event.GetScintProf(pm, evt.getChannel(pm)->getPedestalMean(), maxbin);
            for (int i = 0; i < aux->GetSize(); i++)
               th->Fill(aux->GetBinCenter(i + 1), aux->GetBinContent(i + 1));
            counter++;
            if (counter % 500 == 0)
               cout << "Added " << counter << " wvf." << endl;
            evt.Clear();
            delete aux;
         }
         th->Scale(1.0 / counter);
         cout << "Added " << counter << " events." << endl;
         return th;
      }

      TH1D *TH1ScintProf(int pm, int MaxNumberOfWaveforms, bool DoNotNormalize = false)
      {
         cout << "Computing TH1ScintProf for channel " << pm << ", substracting the pedestal event by event, and shifting to the min bin within range. " << endl;
         cout << "Applying waveform cuts and trigger cuts" << endl;
         GetEntry(0); // midas_chain->GetEntry(0);
         TH1D *th = event.GetWaveform(pm);
         th->Reset();
         th->Sumw2();
         th->SetName(Form("Run%i_OpChannel%i_ScintProf", Number, pm));
         int counter = 0;
         int maxbin = (int)(range1 / Sampling);
         bool first = true;
         //      TCanvas *c= new TCanvas("cscint");
         //      c->Divide(2);
         for (unsigned int ev = 0; ev < NEventsALL && counter < MaxNumberOfWaveforms; ev++)
         {
            Event_t evt = ProcessEvent(ev);
            MyCuts.ApplyMultiChannelCuts(evt);
            if (evt.getCut())
               continue;
            MyCuts.ApplyCuts(*evt.getChannel(pm), pm);
            if (evt.getChannel(pm)->getCut())
               continue;

            GetEntry(ev);
            TH1D *aux = event.GetScintProf(pm, evt.getChannel(pm)->getPedestalMean(), maxbin, range1, range2);
            th->Add(aux);
            counter++;
            //         c->cd(1);//gPad->SetLogy();//gPad->SetGridX();
            //         th->Draw("HIST");th->GetXaxis()->SetRangeUser(3.5e-6,5e-6);
            //         c->cd(2);//gPad->SetLogy();//gPad->SetGridX();
            //         aux->Draw("HIST"); aux->GetXaxis()->SetRangeUser(3.5e-6,5e-6);
            //         c->Update(); c->Modified(); cout << evt.getChannel(pm)->getFirstSampleBelowADCAmplitudeThresholdRange() << endl;
            //         lets_pause();

            if (counter % 500 == 0)
               cout << "Added " << counter << " wvf." << endl;
            delete aux;
         }

         if (!DoNotNormalize)
            th->Scale(1.0 / counter);
         // th->GetXaxis()->SetRangeUser(th->GetMinimum()-0.1,th->GetMaximum()*2);
         cout << "Added " << counter << " events." << endl;
         return th;
      }

      TH1D *TH1ScintProfBeginPeak(int pm, int MaxNumberOfWaveforms = 6000, bool DoNotNormalize = false)
      {

         cout << "Computing TH1ScintProfBeginPeak for channel " << pm << ", substracting the pedestal event by event, and shifting to the min bin within range. MaxNumberOfWaveforms: " << MaxNumberOfWaveforms << endl;
         cout << "Applying waveform cuts and trigger cuts" << endl;
         GetEntry(0); // midas_chain->GetEntry(0);
         TH1D *th = event.GetWaveform(pm);
         th->Reset();
         th->Sumw2();
         th->SetName(Form("Run%i_OpChannel%i_ScintProf", Number, pm));
         int counter = 0;
         int maxbin = (int)(range1 / Sampling);
         bool first = true;
         //      TCanvas *c= new TCanvas("cscint");
         //      c->Divide(2);
         for (unsigned int ev = 0; ev < NEventsALL && counter < MaxNumberOfWaveforms; ev++)
         {
            Event_t evt = ProcessEvent(ev);
            MyCuts.ApplyMultiChannelCuts(evt);
            if (evt.getCut())
               continue;
            MyCuts.ApplyCuts(*evt.getChannel(pm), pm);
            if (evt.getChannel(pm)->getCut())
               continue;
            GetEntry(ev);
            if (evt.getChannel(pm)->getTStartMaxPeakRange() == -1)
               continue;
            TH1D *aux = event.GetScintProfTimeShifted(pm, evt.getChannel(pm)->getPedestalMean(), evt.getChannel(pm)->getTStartMaxPeakRange(), maxbin);
            th->Add(aux);
            counter++; // aux->Draw("HIST"); lets_pause();

            //         c->cd(1);//gPad->SetLogy();//gPad->SetGridX();
            //         th->Draw("HIST");th->GetXaxis()->SetRangeUser(3.5e-6,5e-6);
            //         c->cd(2);//gPad->SetLogy();//gPad->SetGridX();
            //         aux->Draw("HIST"); aux->GetXaxis()->SetRangeUser(3.5e-6,5e-6);
            //         c->Update(); c->Modified(); cout << evt.getChannel(pm)->getFirstSampleBelowADCAmplitudeThresholdRange() << endl;
            //         lets_pause();
            if (counter % 1000 == 0)
               cout << "Added " << counter << " wvf." << endl;
            delete aux;
         }

         if (!DoNotNormalize)
            th->Scale(1.0 / counter);
         // th->GetXaxis()->SetRangeUser(th->GetMinimum()-0.1,th->GetMaximum()*2);
         cout << "Added " << counter << " events." << endl;
         return th;
      }

      TH1D *TH1ScintProfFirstSignalBin(int pm, int MaxNumberOfWaveforms = 6000, bool DoNotNormalize = false)
      {

         cout << "Computing TH1ScintProfFirstSignalBin for channel " << pm << ", substracting the pedestal event by event, and shifting to the min bin within range. MaxNumberOfWaveforms: " << MaxNumberOfWaveforms << endl;
         cout << "Applying waveform cuts and trigger cuts" << endl;
         GetEntry(0); // midas_chain->GetEntry(0);
         TH1D *th = event.GetWaveform(pm);
         th->Reset();
         th->Sumw2();
         th->SetName(Form("Run%i_OpChannel%i_ScintProf", Number, pm));
         int counter = 0;
         int maxbin = (int)(range1 / Sampling);
         bool first = true;
         bool debug = false;
         TCanvas *c;
         if (debug)
         {
            c = new TCanvas("cscint");
            c->Divide(2);
         }
         for (unsigned int ev = 0; ev < NEventsALL && counter < MaxNumberOfWaveforms; ev++)
         {
            Event_t evt = ProcessEvent(ev);
            MyCuts.ApplyMultiChannelCuts(evt);
            if (evt.getCut())
               continue;
            MyCuts.ApplyCuts(*evt.getChannel(pm), pm);
            if (evt.getChannel(pm)->getCut())
               continue;
            GetEntry(ev);
            if (evt.getChannel(pm)->getFirstSampleBelowADCAmplitudeThresholdRange() == -1)
               continue;
            TH1D *aux = event.GetScintProfTimeShifted(pm, evt.getChannel(pm)->getPedestalMean(), evt.getChannel(pm)->getFirstSampleBelowADCAmplitudeThresholdRange(), maxbin);
            th->Add(aux);
            counter++;
            if (debug)
            {
               c->cd(1); // gPad->SetLogy();//gPad->SetGridX();
               th->Draw("HIST");
               c->cd(2); // gPad->SetLogy();//gPad->SetGridX();
               aux->Draw("HIST");
               c->Update();
               c->Modified();
               cout << evt.getChannel(pm)->getFirstSampleBelowADCAmplitudeThresholdRange() << endl;
               lets_pause();
            }
            if (counter % 1000 == 0)
               cout << "Added " << counter << " wvf." << endl;
            delete aux;
         }

         if (!DoNotNormalize)
            th->Scale(1.0 / counter);
         // th->GetXaxis()->SetRangeUser(th->GetMinimum()-0.1,th->GetMaximum()*2);
         cout << "Added " << counter << " events." << endl;
         return th;
      }

      std::map<int, TH1D> SiPMResponse;
      void LoadSiPMResponse(int ch)
      {
         TFile *f = new TFile("AnalysisROOT/AvWf_SPE.root", "READ");
         SiPMResponse[ch] = *((TH1D *)f->Get(Form("Run4_ch%i_ADC%i_0V_ScintProfFirstSignalBin_0_", ch, ch + 1))->Clone(Form("SiPMResponse_ch%i", ch)));
         //       SiPMResponse[ch].Draw("HIST"); lets_pause();
         // if(!SiPMResponse[ch]) { cout << "Error, SiPM response not found!!!" << std::endl;throw std::exception();}
         f->Close();
      }
      TH1D *GetSiPMResponse(int ch)
      {
         if (SiPMResponse.find(ch) == SiPMResponse.end())
            LoadSiPMResponse(ch);
         return &SiPMResponse[ch];
      }
      TH1D *TH1ScintProfDeconvoluteSiPM(int pm, int MaxNumberOfWaveforms = 6000)
      {

         cout << "Computing TH1ScintProfDeconvoluteSiPM for channel " << pm << ", substracting the pedestal event by event, and shifting to the min bin within range. " << endl;
         cout << "Applying waveform cuts and trigger cuts" << endl;
         GetEntry(0); // midas_chain->GetEntry(0);
         TH1D *th = event.GetWaveform(pm);
         th->Reset();
         th->Sumw2();
         th->SetName(Form("Run%i_OpChannel%i_ScintProf", Number, pm));
         int counter = 0;
         int maxbin = (int)(range1 / Sampling);
         bool first = true;

         NoiseTools_t noise;
         TH1D *res = GetSiPMResponse(pm);
         res->Draw("HIST");
         lets_pause();
         noise.MovingAverageFilter(res, 15);
         int maxbin2 = res->FindFirstBinAbove(5);
         TSpectrum *s = new TSpectrum();
         int nbins = th->GetSize() - 2;
         Double_t source[nbins];
         Double_t response[nbins];
         for (int i = 0; i < nbins; i++)
            response[i] = 0.0;
         for (int i = 0; i < nbins; i++)
            source[i] = 0.0;
         for (int i = 0; i < res->GetSize() - 2; i++)
         {
            if (res->GetBinContent(maxbin2 - 2 + i) > 0.1)
            {
               response[i] = res->GetBinContent(maxbin2 - 2 + i);
            }
         }
         for (int i = 0; i < res->GetSize() - 2; i++)
         {
            if (res->GetBinContent(i + 1) > 0.1)
            {
               source[i + 200] = res->GetBinContent(i + 1);
            }
         }

         TCanvas *c = new TCanvas("cconv");
         c->Divide(2);
         for (unsigned int ev = 0; ev < NEventsALL && counter < MaxNumberOfWaveforms; ev++)
         {
            Event_t evt = ProcessEvent(ev);
            MyCuts.ApplyMultiChannelCuts(evt);
            if (evt.getCut())
               continue;
            MyCuts.ApplyCuts(*evt.getChannel(pm), pm);
            if (evt.getChannel(pm)->getCut())
               continue;

            GetEntry(ev);
            TH1D *aux = event.GetScintProfTimeShifted(pm, evt.getChannel(pm)->getPedestalMean(), evt.getChannel(pm)->getFirstSampleBelowADCAmplitudeThresholdRange(), maxbin);
            noise.MovingAverageFilter(aux, 20);
            //         aux->Draw("HIST"); //lets_pause();
            for (int i = 0; i < nbins; i++)
            {
               if (aux->GetBinContent(i + 1) > 0.1)
               {
                  source[i] = aux->GetBinContent(i + 1);
               }
               else
               {
                  source[i] = 0.0;
               }
            }
            //         cout << s->Deconvolution(source,response,256,1000,1,1) << endl;
            TH1F *d = (TH1F *)aux->Clone("DeconvResponse");
            d->SetTitle("Source");
            d->SetLineColor(kRed);
            TH1F *d2 = (TH1F *)aux->Clone("DeconvResponse");
            d2->SetTitle("Deconvolution");
            d2->SetLineColor(kGreen);
            //         for(unsigned int i = 0; i < nbins; i++) {d->SetBinContent(i + 1,source[i]); cout << i << " " << source[i] << " " << response[i] << endl;}
            //         for(unsigned int i = 0; i < nbins; i++) {d2->SetBinContent(i + 1,response[i]);}
            //         d->Draw("HIST L");
            //         d2->Draw("HIST SAME L");
            //         gPad->BuildLegend();
            //         lets_pause();

            auto result = s->Deconvolution(source, response, nbins, 1000, 1, 1);
            //   s->DeconvolutionRL(source,response,256,200,50,1.2);
            for (int i = 0; i < nbins; i++)
            {
               d2->SetBinContent(i + 1, source[i]);
            }

            //         cout << result << endl;

            th->Add(d2);
            counter++; // th->Draw("HIST"); lets_pause();
            c->cd(1);
            res->Draw("HIST");
            d->Draw("HIST L");
            d2->Draw("HIST SAME L");
            gPad->BuildLegend();
            c->cd(2);
            gPad->SetLogy();
            th->Draw("HIST");
            c->Modified();
            c->Update();
            //         if(counter%1000==0)
            lets_pause();

            if (counter % 100 == 0)
               cout << "Added " << counter << " wvf." << endl;
            delete aux;
            delete d2;
         }

         th->Scale(1.0 / counter);
         th->GetXaxis()->SetRangeUser(th->GetMinimum() - 0.1, th->GetMaximum() * 2);
         cout << "Added " << counter << " events." << endl;
         return th;
      }

      TH1D *TH1ScintProfMyDeconvolution(int pm, int MaxNumberOfWaveforms = 100000)
      {
         MaxNumberOfWaveforms = 100000;
         cout << "Computing TH1ScintProfDeconvoluteSiPM for channel " << pm << ", substracting the pedestal event by event, and shifting to the min bin within range. " << endl;
         cout << "Applying waveform cuts and trigger cuts" << endl;
         GetEntry(0); // midas_chain->GetEntry(0);
         TH1D *th = event.GetWaveform(pm);
         th->Reset();
         th->Sumw2();
         th->SetName(Form("Run%i_OpChannel%i_ScintProf", Number, pm));
         int counter = 0;
         int maxbin = (int)(range1 / Sampling);
         bool first = true;
         TH1D *res = GetSiPMResponse(pm);
         std::vector<double> myres;
         for (int i = res->GetMaximumBin(); i < res->GetSize() - 2; i++)
            myres.push_back(res->GetBinContent(i));
         Decon_t dd(myres);

         for (unsigned int ev = 0; ev < NEventsALL && counter < MaxNumberOfWaveforms; ev++)
         {
            Event_t evt = ProcessEvent(ev);
            MyCuts.ApplyMultiChannelCuts(evt);
            if (evt.getCut())
               continue;
            MyCuts.ApplyCuts(*evt.getChannel(pm), pm);
            if (evt.getChannel(pm)->getCut())
               continue;

            GetEntry(ev);
            TH1D *aux = event.GetScintProfTimeShifted(pm, evt.getChannel(pm)->getPedestalMean(), evt.getChannel(pm)->getFirstSampleBelowADCAmplitudeThresholdRange(), maxbin);
            TH1D *aux2 = dd.Deconvolute(aux);
            aux2->SetLineColor(kRed);
            //         aux->Draw("HIST"); //lets_pause();
            //         aux2->Draw("HIST"); lets_pause();
            //         for(unsigned int i = 0; i < nbins; i++) source[i]=aux->GetBinContent(i + 1);
            //         cout << s->Deconvolution(source,response,256,1000,1,1) << endl;
            //         TH1F *d = (TH1F*) aux->Clone("DeconvResponse"); d->SetTitle("Source");d->SetLineColor(kRed);
            //         TH1F *d2 = (TH1F*) aux->Clone("DeconvResponse"); d2->SetTitle("Deconvolution");d2->SetLineColor(kGreen);
            //         for(unsigned int i = 0; i < nbins; i++) {d->SetBinContent(i + 1,source[i]); cout << i << " " << source[i] << " " << response[i] << endl;}
            //         for(unsigned int i = 0; i < nbins; i++) {d2->SetBinContent(i + 1,response[i]);}
            //         d->Draw("HIST L");
            //         d2->Draw("HIST SAME L");
            //         gPad->BuildLegend();
            //         lets_pause();

            //         s->Deconvolution(source,response,nbins,5000,10,3);
            //   s->DeconvolutionRL(source,response,256,200,50,1.2);
            //         for(unsigned int i = 0; i < nbins; i++) {d2->SetBinContent(i + 1,source[i]);}

            //         res->Draw("HIST");
            //         d->Draw("HIST L");
            //         d2->Draw("HIST SAME L");
            //         gPad->BuildLegend();
            //         lets_pause();

            th->Add(aux2);
            counter++; // aux->Draw("HIST"); lets_pause();

            if (counter % 500 == 0)
               cout << "Added " << counter << " wvf." << endl;
            delete aux;
            delete aux2;
         }

         th->Scale(1.0 / counter);
         th->GetXaxis()->SetRangeUser(th->GetMinimum() - 0.1, th->GetMaximum() * 2);
         cout << "Added " << counter << " events." << endl;
         return th;
      }

      TH1D *TH1ScintProfTriggerCut(int pm, int MaxNumberOfWaveforms = -1)
      {
         cout << "Computing TH1ScintProf for channel " << pm << ", substracting the pedestal event by event, and shifting to the trigger time!. " << endl;
         cout << "Applying waveform cuts and trigger cuts" << endl;

         GetEntry(0);
         TH1D *th = event.GetWaveform(pm);
         th->Reset();
         th->Sumw2();
         th->SetName(Form("Run%i_OpChannel%i_ScintProfTriggerCut", Number, pm));
         if (MaxNumberOfWaveforms == -1 || MaxNumberOfWaveforms > (int)NEventsALL)
            MaxNumberOfWaveforms = NEventsALL;
         int counter = 0;
         int maxbin = (int)(range1 / Sampling);
         bool first = true;
         for (unsigned int ev = 0; ev < NEventsALL && counter < MaxNumberOfWaveforms; ev++)
         {
            Event_t evt = ProcessEvent(ev);
            MyCuts.ApplyMultiChannelCuts(evt);
            if (evt.getCut())
               continue;
            MyCuts.ApplyCuts(*evt.getChannel(pm), pm);
            if (evt.getChannel(pm)->getCut())
               continue;

            GetEntry(ev);
            //         cout << "Event " << ev << " passed cut: TS: " << evt.getTriggerTime() << endl;
            TH1D *aux = event.GetScintProfTimeShifted(pm, PedestalPeak[pm], evt.getTriggerTime(), maxbin); // aux->Draw("HIST"); lets_pause();
            th->Add(aux);
            counter++; // th->Draw("HIST"); lets_pause();
            if (counter % 500 == 0)
               cout << "Added " << counter << " wvf." << endl;
            delete aux;
         }
         th->Scale(1.0 / counter);
         th->GetXaxis()->SetRangeUser(th->GetMinimum() - 0.1, th->GetMaximum() * 2);
         cout << "Added " << counter << " events." << endl;
         return th;
      }

      TH1D *TH1ScintProf2(int pm, int MaxNumberOfWaveforms = -1)
      {
         cout << "Computing TH1ScintProf for channel " << pm << ", assuming a stable pedestal. Always substracting the average pedestal." << endl;
         cout << "Applying waveform cuts, not trigger cuts" << endl;
         GetEntry(0); // midas_chain->GetEntry(0);
         TH1D *th = event.GetWaveform(pm);
         th->Reset();
         th->Sumw2();
         th->SetName(Form("Run%i_OpChannel%i_ScintProf2", Number, pm));
         int counter = 0;
         int maxbin = (int)(range1 / Sampling);
         bool first = true;
         if (MaxNumberOfWaveforms == -1 || MaxNumberOfWaveforms > (int)NEventsALL)
            MaxNumberOfWaveforms = NEventsALL;
         for (int ev = 0; ev < MaxNumberOfWaveforms; ev++)
         {
            GetEntry(ev);
            waveana::Waveform_t myevt(event.GetTime(), event.GetAmp(pm), ParSet);
            MyCuts.ApplyCuts(myevt, pm);
            if (myevt.getCut())
               continue;
            TH1D *aux = event.GetScintProf(pm, PedestalPeak[pm], maxbin);
            th->Add(aux);
            counter++;
            if (counter % 500 == 0)
               cout << "Added " << counter << " wvf." << endl;
            delete aux;
         }
         // th->Draw(); lets_pause();
         th->Scale(1.0 / counter);
         th->GetXaxis()->SetRangeUser(th->GetMinimum() - 0.1, th->GetMaximum() * 2);
         //      th->Scale(1.0/th->GetBinContent(th->GetMaximumBin()));;
         // th->Draw(); lets_pause();
         cout << "Added " << counter << " events." << endl;
         //      ThreeExpoGausFitScintUtilsCh ChiaraFunction; ChiaraFunction.Fit(*th);
         return th;
      }

      TH1D *TH1ScintProfALL()
      {
         cout << "Computing TH1ScintProf for ALL channels." << endl;
         GetEntry(0);
         TH1D *th = event.GetWaveform(0);
         th->Reset();
         th->Sumw2(); // th->Draw(); lets_pause();
         th->SetName(Form("Run%i_ScintProfALL", Number));
         int counter = 0;
         int maxbin = (int)(range1 / Sampling);
         bool first = true;
         for (unsigned int ev = 0; ev < NEventsALL && counter < 50000000; ev++)
         {
            GetEntry(ev);
            for (auto pm : adcchannels)
            {
               waveana::Waveform_t myevt(event.GetTime(), event.GetAmp(pm), ParSet);
               MyCuts.ApplyCuts(myevt, pm);
               if (myevt.getCut())
                  continue;
               //         if(first) {th = event.GetWaveform(pm); maxbin = (int)(0.7e-6/sampling);/*maxbin=th->FindBin(getWaveform(ev,pm)->getPeakTime());*/ th->Reset(); first=false;}
               TH1D *aux = event.GetScintProf(pm, myevt.getPedestalMean(), maxbin);
               th->Add(aux);
               counter++;
               if (counter % 5000 == 0)
                  cout << "Added " << counter << " wvf." << endl;
               delete aux;
            }
         }
         th->Scale(1.0 / th->GetBinContent(th->GetMaximumBin()));
         ;
         cout << "Added " << counter << " events." << endl;
         //      ThreeExpoGausFitScintUtilsCh ChiaraFunction; ChiaraFunction.Fit(*th);
         return th;
      }
      void PlotFirstSampleBelowADCTriggerThresholdRange()
      {
         TH1F *th[64];
         for (auto pm : adcchannels)
         {
            th[pm] = TH1FirstSampleBelowADCTriggerThresholdRange(pm);
         }
         TCanvas *c = new TCanvas("c");
         c->cd(); // th[adcchannels[0]]->Draw("HIST");
         for (auto pm : adcchannels)
            th[pm]->Draw("HIST SAME");
         lets_pause();
      }
      void PlotFirstSampleAboveADCTriggerThresholdRange()
      {
         TH1F *th[64];
         for (auto pm : adcchannels)
         {
            th[pm] = TH1FirstSampleAboveADCTriggerThresholdRange(pm);
         }
         TCanvas *c = new TCanvas("c");
         c->cd(); // th[adcchannels[0]]->Draw("HIST");
         for (auto pm : adcchannels)
            th[pm]->Draw("HIST SAME");
         lets_pause();
      }
      void PlotPeakTimes()
      {

         TH1F *th[64];
         for (auto pm : adcchannels)
         {
            th[pm] = TH1PeakTime(pm);
         }
         TCanvas *c = new TCanvas("c");
         c->cd(); // th[adcchannels[0]]->Draw("HIST");
         for (auto pm : adcchannels)
            th[pm]->Draw("HIST SAME");
         lets_pause();
      }
      TH1F *TH1MaximumSampleRange(int pm)
      {
         TH1F *th = new TH1F(Form("Run%i_OpChannel%i_MaximumSampleRange", Number, pm), Form("Run%i_OpChannel%i_MaximumSampleRange;MaximumSample (ADC);Events", Number, pm), ADCDynamicRange, 0, ADCDynamicRange);
         for (unsigned int ev = 0; ev < NEvents; ev++)
            if (!getWaveform(ev, pm)->getCut())
            {
               th->Fill(getWaveform(ev, pm)->getMaximumSample_Range());
            }
         return th;
      }
      TH1F *TH1PreTriggerSTD(int pm)
      {
         TH1F *th = new TH1F(Form("Run%i_OpChannel%i_PreTriggerSTD", Number, pm), Form("Run%i_OpChannel%i_PreTriggerSTD;PreTriggerSTD (ADC);Events", Number, pm), 20000, 0, 1000);
         for (unsigned int ev = 0; ev < NEvents; ev++)
            if (!getWaveform(ev, pm)->getCut())
            {
               th->Fill(getWaveform(ev, pm)->getPreTriggerSTD());
               cout << getWaveform(ev, pm)->getPreTriggerSTD() << endl;
            }
         return th;
      }
      TH1F *TH1PreTriggerCharge(int pm)
      {
         TH1F *th = new TH1F(Form("Run%i_OpChannel%i_PreTriggerCharge", Number, pm), Form("Run%i_OpChannel%i_PreTriggerCharge;PreTriggerCharge (C);Events", Number, pm), 0, 0, 0);
         for (unsigned int ev = 0; ev < NEvents; ev++)
            if (!getWaveform(ev, pm)->getCut())
            {
               th->Fill(getWaveform(ev, pm)->getPreTriggerCharge());
            }
         return th;
      }
      TH1F *TH1PeakTime(int pm)
      {
         TH1F *th = new TH1F(Form("Run%i_OpChannel%i_PeakTime", Number, pm), Form("Run%i_OpChannel%i_PeakTime;PeakTime (s);Events", Number, pm), getWaveform(0, pm)->getNSamples(), 0, getWaveform(0, pm)->getNSamples() * getWaveform(0, pm)->getSampling());
         for (unsigned int ev = 0; ev < NEvents; ev++)
            if (!getWaveform(ev, pm)->getCut())
            {
               th->Fill(getWaveform(ev, pm)->getPeakTime());
            }
         return th;
      }
      TH1F *TH1FirstSampleBelowADCTriggerThresholdRange(int pm)
      {
         TH1F *th = new TH1F(Form("Run%i_OpChannel%i_FirstSampleBelowADCTriggerThresholdRange", Number, pm), Form("Run%i_OpChannel%i_FirstSampleBelow%.0fADC;FirstSampleBelow%.0fADC (s);Events", Number, pm, ParSet->getADCTriggerThreshold(), ParSet->getADCTriggerThreshold()), getWaveform(0, pm)->getNSamples(), 0, getWaveform(0, pm)->getNSamples() * getWaveform(0, pm)->getSampling());
         for (unsigned int ev = 0; ev < NEvents; ev++)
            if (!getWaveform(ev, pm)->getCut())
            {
               th->Fill(getWaveform(ev, pm)->getFirstSampleBelowADCTriggerThresholdRange());
            }
         return th;
      }
      TH1F *TH1FirstSampleAboveADCTriggerThresholdRange(int pm)
      {
         TH1F *th = new TH1F(Form("Run%i_OpChannel%i_FirstSampleAboveADCTriggerThresholdRange", Number, pm), Form("Run%i_OpChannel%i_FirstSampleAbove%.0fADC;FirstSampleAbove%.0fADC (s);Events", Number, pm, ParSet->getADCTriggerThreshold(), ParSet->getADCTriggerThreshold()), getWaveform(0, pm)->getNSamples(), 0, getWaveform(0, pm)->getNSamples() * getWaveform(0, pm)->getSampling());
         for (unsigned int ev = 0; ev < NEvents; ev++)
            if (!getWaveform(ev, pm)->getCut())
            {
               th->Fill(getWaveform(ev, pm)->getFirstSampleAboveADCTriggerThresholdRange());
            }
         return th;
      }
      TH1F *TH1ChargeInPEUnits(int pm)
      {
         TH1F *thaux = new TH1F(Form("Run%i_OpChannel%i_Charge_aux", Number, pm), Form("Run%i_OpChannel%i_%s_%.0fV_pedSTD<%.1fADC_aux;Charge (PE);Events", Number, pm, PMT_SN[pm].c_str(), PMT_Voltages[pm], PedestalSTDCUT[pm]), 5000, -2e6, 1e6);
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            //           double correction=(PedestalPeak[pm]-EventList[pm][i].getPedestalMean())*(2.0/4096)*Sampling*NSamples/50;
            if (!getWaveform(ev, pm)->getCut())
               thaux->Fill((getWaveform(ev, pm)->getCharge()) / 1.602e-19 / PMT_Gains[pm]);
         }
         Double_t xq2[2]; // position where to compute the quantiles in [0,1]
         Double_t yq2[2]; // array to contain the quantiles
         xq2[0] = 0.001;
         xq2[1] = 0.999;
         thaux->GetQuantiles(2, yq2, xq2);

         TH1F *th = new TH1F(Form("Run%i_OpChannel%i_ChargePE", Number, pm), Form("Run%i_OpChannel%i_%s_%.0fV_pedSTD<%.1fADC;Charge (PE);Events", Number, pm, PMT_SN[pm].c_str(), PMT_Voltages[pm], PedestalSTDCUT[pm]), 300, 2 * yq2[0], 1.2 * yq2[1]);
         th->SetLineColor(pm + 1);
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            //           double correction=(PedestalPeak[pm]-EventList[pm][i].getPedestalMean())*(2.0/4096)*Sampling*NSamples/50;
            if (!getWaveform(ev, pm)->getCut())
               th->Fill((getWaveform(ev, pm)->getCharge()) / 1.602e-19 / PMT_Gains[pm]);
         }

         delete thaux;
         return th;
      }

      void PlotPedestals()
      {

         TH1F *th[64];
         for (auto pm : adcchannels)
         {
            th[pm] = TH1Pedestal(pm);
         }
         TCanvas *c = new TCanvas("c");
         c->cd();
         for (auto pm : adcchannels)
            th[pm]->Draw("HIST SAME");
         lets_pause();
      }
      TH1F *TH1Pedestal(int pm)
      {
         TH1F *thaux = new TH1F(Form("Run%i_OpChannel%i_Ped_test", Number, pm), Form("Run%i_OpChannel%i_Ped;Pedestal (ADC);Events", Number, pm), 1000, 0, 0);
         for (unsigned int ev = 0; ev < NEvents; ev++)
            if (!getWaveform(ev, pm)->getCut())
            {
               thaux->Fill(getWaveform(ev, pm)->getPedestalMean());
            }

         Double_t xq2[2]; // position where to compute the quantiles in [0,1]
         Double_t yq2[2]; // array to contain the quantiles
         xq2[0] = 0.001;
         xq2[1] = 0.999;
         thaux->GetQuantiles(2, yq2, xq2);

         float min = yq2[0] - 60;
         float max = yq2[1] + 60;
         int NBins = (max - min) * (PedRange);
         TH1F *th = new TH1F(Form("Run%i_OpChannel%i_Ped", Number, pm), Form("Run%i_OpChannel%i_Ped;Pedestal (ADC);Events", Number, pm), NBins, min, max);
         th->SetLineColor(pm + 1);
         for (unsigned int ev = 0; ev < NEvents; ev++)
            if (!getWaveform(ev, pm)->getCut())
            {
               th->Fill(getWaveform(ev, pm)->getPedestalMean());
            }

         delete thaux;
         return th;
      }
      TH1F *TH1MaximumSample(int pm)
      {
         TH1F *th = new TH1F(Form("Run%i_OpChannel%i_MaximumSample", Number, pm), Form("Run%i_OpChannel%i_MaximumSample;MaximumSample (ADC);Events", Number, pm), ADCDynamicRange, 0, ADCDynamicRange);
         th->SetLineColor(pm + 1);
         for (unsigned int ev = 0; ev < NEvents; ev++)
            if (!getWaveform(ev, pm)->getCut())
            {
               th->Fill(getWaveform(ev, pm)->getMaximumSample());
            }
         return th;
      }
      TH2F *TH2PedestalVsPedSTD(int pm)
      {
         TH2F *th = new TH2F(Form("Run%i_OpChannel%i_PedestalVsPedSTD", Number, pm), Form("Run%i_OpChannel%i_PedestalVsPedSTD;PedSTD (ADC);Pedestal (ADC)", Number, pm), 200, 0, 20, 20000, 3500, 4100);
         th->SetLineColor(pm + 1);
         for (unsigned int ev = 0; ev < NEvents; ev++)
            if (!getWaveform(ev, pm)->getCut())
            {
               th->Fill(getWaveform(ev, pm)->getPedestalSTD(), getWaveform(ev, pm)->getPedestalMean());
            }
         return th;
      }
      void PlotPedSTDs()
      {
         std::cout << "Plotting PedSTDs..." << std::endl;
         TH1F *th[64];
         for (auto pm : adcchannels)
         {
            th[pm] = TH1PedSTD(pm);
         }
         TCanvas *c = new TCanvas("c");
         c->cd();
         for (auto pm : adcchannels)
            th[pm]->Draw("HIST SAME");
         lets_pause();
      }
      TH1F *TH1PedSTD(int pm)
      {
         TH1F *th = new TH1F(Form("Run%i_OpChannel%i_PedSTD", Number, pm), Form("Run%i_OpChannel%i_PedSTD;Pedestal STD (ADC);Events", Number, pm), 2000, 0, 100);
         th->SetLineColor(pm + 1);
         for (unsigned int ev = 0; ev < NEvents; ev++)
            if (!getWaveform(ev, pm)->getCut())
            {
               th->Fill(getWaveform(ev, pm)->getPedestalSTD());
            }
         return th;
      }

      void PlotChargeSPE(double rangemin = -2.e7, double rangemax = 1.e8)
      {
         TH1F *th[64];
         for (auto pm : adcchannels)
         {
            th[pm] = TH1ChargeSPE(pm, "Range", rangemin, rangemax);
         }
         TCanvas *c = new TCanvas("c");
         c->cd();
         for (auto pm : adcchannels)
            th[pm]->Draw("HIST SAME");
         lets_pause();
      }

      TH1F *TH1ChargeSPE(int pm, string option = "Range", double xmin = 0, double xmax = 0)
      {
         if (xmin == 0 && xmax == 0)
         {
            TH1F *thaux = new TH1F(Form("Run%i_OpChannel%i_Charge_aux", Number, pm), Form("Run%i_OpChannel%i_%s_%.0fV_pedSTD<%.1fADC_aux;Charge (/e/);Events", Number, pm, PMT_SN[pm].c_str(), PMT_Voltages[pm], PedestalSTDCUT[pm]), 5000, -2e8, 1e10);
            for (unsigned int ev = 0; ev < NEvents; ev++)
            {
               if (getWaveform(ev, pm)->getCut())
                  continue;
               if (option == "3bins")
                  thaux->Fill(getWaveform(ev, pm)->getChargeMaxPeakRange3bins() / 1.602e-19);
               if (option == "Range")
                  thaux->Fill(getWaveform(ev, pm)->getChargeRange() / 1.602e-19);
               if (option == "Peak")
                  thaux->Fill(getWaveform(ev, pm)->getChargeMaxPeakRange() / 1.602e-19);
            }

            Double_t xq2[2]; // position where to compute the quantiles in [0,1]
            Double_t yq2[2]; // array to contain the quantiles
            xq2[0] = 0.001;
            xq2[1] = 0.999;
            thaux->GetQuantiles(2, yq2, xq2);
            xmax = 1.15 * yq2[1];
            xmin = yq2[0] - 0.15 * yq2[1];
            delete thaux;
         }
         double BinWidth = Sampling * 1.0 * (2.0 / ADCDynamicRange) / 50.0 / 1.602e-19; // resolution of 1 tickxADC in # of electrons.
         int NBins = (xmax - xmin) / BinWidth;                                          // cout << "xmin, xmax: "<< xmin << " " << xmax<< " - NBins: " << NBins << endl; //maximum number of bins
         // cout << BinWidth << " " << NBins << " " << xmin << " " << xmax << endl; lets_pause();
         if (NBins > 600)
            NBins = NBins / (int)(1.0 * NBins / 300.0); // Rebin to 300 avoiding binning issues.
         TH1F *th = new TH1F(Form("Run%i_ch%i_SPE", Number, pm), Form("Run%i_ch%i_%s_%.0fV_SPE%s;Charge (/e/);Events", Number, pm, PMT_SN[pm].c_str(), PMT_Voltages[pm], option.c_str()), NBins, xmin, xmax);
         th->SetLineColor(pm + 1);
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            if (getWaveform(ev, pm)->getCut())
               continue;
            if (option == "3bins")
               th->Fill(getWaveform(ev, pm)->getChargeMaxPeakRange3bins() / 1.602e-19);
            if (option == "Range")
               th->Fill(getWaveform(ev, pm)->getChargeRange() / 1.602e-19);
            if (option == "Peak")
               th->Fill(getWaveform(ev, pm)->getChargeMaxPeakRange() / 1.602e-19);
         }
         return th;
      }
      TH1F *TH1MaxAmplitudeOutOfMaxPeak(int pm)
      {
         TH1F *th = new TH1F(Form("Run%i_OpChannel%i_MaxAmplitudeOutOfMaxPeak", Number, pm), Form("Run%i_OpChannel%i_%s_Amp;Amplitude (ADC);Events", Number, pm, PMT_SN[pm].c_str()), 1000, 0, ADCDynamicRange);
         th->SetLineColor(pm + 1);
         for (unsigned int ev = 0; ev < NEvents; ev++)
            if (!getWaveform(ev, pm)->getCut())
            {
               th->Fill(getWaveform(ev, pm)->getMaxAmplitudeOutOfMaxPeak());
            }
         th->SetLineWidth(2.0);
         return th;
      }
      int MaximumWaveformsToProcess = 6000;
      void SetMaximumWaveformsToProcess(int n)
      {
         if (n == -1)
            MaximumWaveformsToProcess = NEventsALL;
         else
            MaximumWaveformsToProcess = n;
      }
      HistogramCollection_t *GetHistCollectionByMode(string mode, double range1 = 0, double range2 = 0, std::vector<int> *SelectedChannels = NULL)
      {
         std::cout << "Creating HistogramCollection of " << mode << endl;
         HistogramCollection_t *HC = new HistogramCollection_t();
         int i = 1;
         for (auto pm : *SelectedChannels)
         {
            TH1 *th = NULL;
            if (mode == "Amp")
            {
               th = TH1Amp(pm);
            }
            if (mode == "AmpRange")
            {
               th = TH1AmpRange(pm, "ADC", new std::vector<float>({static_cast<float>(range1), static_cast<float>(range2)}));
            }
            if (mode == "AmpRange_PE")
            {
               th = TH1AmpRange(pm, "PE");
            }
            if (mode == "Amp2")
            {
               th = TH1Amp2(pm);
            }
            if (mode == "ChargePE")
            {
               th = TH1ChargeInPEUnits(pm);
            }
            if (mode == "Charge")
            {
               th = TH1Charge(pm, "", "pC", range1, range2);
            }
            if (mode == "ChargeRange")
            {
               th = TH1Charge(pm, "Range", "pC", range1, range2);
            }
            if (mode == "ChargeRangePE")
            {
               th = TH1Charge(pm, "Range", "PE", range1, range2);
            }
            if (mode == "Charge_Q2MaxPeakRange")
            {
               th = TH1Charge(pm, "Q2MaxPeakRange", "pC", range1, range2);
            }
            if (mode == "Charge_Q2MaxPeakRange_PE")
            {
               th = TH1Charge(pm, "Q2MaxPeakRange", "PE", range1, range2);
            }
            if (mode == "Charge_Q3MaxPeakRange")
            {
               th = TH1Charge(pm, "Q3MaxPeakRange", "pC", range1, range2);
            }
            if (mode == "Charge_Q3MaxPeakRange_PE")
            {
               th = TH1Charge(pm, "Q3MaxPeakRange", "PE", range1, range2);
            }
            if (mode == "Charge_MaxPeakRange_PE")
            {
               th = TH1Charge(pm, "MaxPeakRange", "PE", range1, range2);
            }
            if (mode == "Charge_MaxPeakRange")
            {
               th = TH1Charge(pm, "MaxPeakRange", "pC", range1, range2);
            }
            if (mode == "AvWf")
            {
               th = TH1AvWf(pm);
            }
            if (mode == "AmpVsCharge")
            {
               th = TH2AmpVsCharge(pm);
               HC->SetDrawOption("COLZ");
            }
            if (mode == "AmpVsCharge_Range")
            {
               th = TH2AmpVsCharge(pm, "Range");
               HC->SetDrawOption("COLZ");
            }
            if (mode == "AmpVsCharge_MaxPeak")
            {
               th = TH2AmpVsCharge(pm, "MaxPeak");
               HC->SetDrawOption("COLZ");
            }
            if (mode == "AmpVsCharge_MaxPeakRange")
            {
               th = TH2AmpVsCharge(pm, "MaxPeakRange");
               HC->SetDrawOption("COLZ");
            }
            if (mode == "AmpVsCharge_Q1MaxPeakRange")
            {
               th = TH2AmpVsCharge(pm, "Q1MaxPeakRange");
               HC->SetDrawOption("COLZ");
            }
            if (mode == "AmpVsCharge_Q2MaxPeakRange")
            {
               th = TH2AmpVsCharge(pm, "Q2MaxPeakRange");
               HC->SetDrawOption("COLZ");
            }
            if (mode == "AmpVsCharge_Q3MaxPeakRange")
            {
               th = TH2AmpVsCharge(pm, "Q3MaxPeakRange");
               HC->SetDrawOption("COLZ");
            }
            if (mode == "Range34VsRange56")
            {
               th = TH2Range34vsRange56(pm);
               HC->SetDrawOption("COLZ");
            }
            if (mode == "PeakTime")
            {
               th = TH1PeakTime(pm);
            }
            if (mode == "ScintProf")
            {
               th = TH1ScintProf(pm, MaximumWaveformsToProcess);
            }
            if (mode == "ScintProf_NotNormalize")
            {
               th = TH1ScintProf(pm, MaximumWaveformsToProcess, true);
            }
            if (mode == "ScintProf2")
            {
               th = TH1ScintProf2(pm, MaximumWaveformsToProcess);
            }
            if (mode == "ScintProfFirstSignalBin")
            {
               th = TH1ScintProfFirstSignalBin(pm, MaximumWaveformsToProcess);
            }
            if (mode == "ScintProfBeginPeak")
            {
               th = TH1ScintProfBeginPeak(pm, MaximumWaveformsToProcess);
            }
            if (mode == "ScintProfFirstSignalBin_NotNormalize")
            {
               th = TH1ScintProfFirstSignalBin(pm, MaximumWaveformsToProcess, true);
            }
            if (mode == "ScintProfDeconvoluteSiPM")
            {
               th = TH1ScintProfDeconvoluteSiPM(pm, MaximumWaveformsToProcess);
            }
            if (mode == "ScintProfMyDeconvolution")
            {
               th = TH1ScintProfMyDeconvolution(pm, MaximumWaveformsToProcess);
            }
            if (mode == "EvolutionPedestal")
            {
               th = TProfileEvolution(pm, "pedestal", 100);
               HC->SetDrawOption("E1");
            }
            if (mode == "EvolutionMaxAmplitude")
            {
               th = TProfileEvolution(pm, "MaxAmplitude", 100);
               HC->SetDrawOption("E1");
            }
            if (mode == "EvolutionMaxAmplitudeRange")
            {
               th = TProfileEvolution(pm, "MaxAmplitudeRange", 100);
               HC->SetDrawOption("E1");
            }
            if (mode == "MaxAmplitudeOutOfMaxPeak")
            {
               th = TH1MaxAmplitudeOutOfMaxPeak(pm);
            }
            if (mode == "FreqMinimum")
            {
               th = TH1FreqMinimum(pm);
            }
            if (mode == "Freq")
            {
               th = TH1Freq(pm);
            }
            if (mode == "ChargeSPERange")
            {
               th = TH1ChargeSPE(pm, "Range");
            }
            if (mode == "ChargeSPE3bins")
            {
               th = TH1ChargeSPE(pm, "3bins");
            }
            if (mode == "ChargeSPEPeak")
            {
               th = TH1ChargeSPE(pm, "Peak");
            }
            if (mode == "MaximumSampleRange")
            {
               th = TH1MaximumSampleRange(pm);
            }
            if (mode == "PreTriggerCharge")
            {
               th = TH1PreTriggerCharge(pm);
            }
            if (mode == "PreTriggerSTD")
            {
               th = TH1PreTriggerSTD(pm);
            }

            if (mode == "TH2AverageSignal")
            {
               th = TH2AverageSignal(pm);
               HC->SetDrawOption("COLZ");
            }
            if (mode == "TH2DeltaTimes")
            {
               th = TH2DeltaTimes(pm);
               HC->SetDrawOption("");
            }
            if (!th)
            {
               std::cout << "ERROR: Plot36(mode=" << mode << ") not implemented in this software version. " << std::endl;
               throw std::exception();
            }
            th->SetLineColor(i);
            if (i == 9)
            {
               i = 20;
            }
            else
            {
               i++;
            }
            HC->Add(th, Number, pm, PMT_SN[pm], PMT_Voltages[pm], mode, firsttime, "");
         }
         return HC;
      }

      void ExportScintProf(string dumpfile)
      {
         TH1 *th = TH1ScintProfALL();
         th->SetTitle(Form("Run%i", Number));
         th->SetName(Form("Run%i", Number));
         TFile *dump = new TFile(dumpfile.c_str(), "UPDATE");
         dump->cd();
         th->Write();
         dump->Close();
         delete th;
      }
      TH1F *TH1Amp(int pm)
      {
         TH1F *th = new TH1F(Form("Run%i_OpChannel%i_Amp", Number, pm), Form("Run%i_OpChannel%i_%s_Amp;Amplitude (ADC);Events", Number, pm, PMT_SN[pm].c_str()), 500, 0.5, 2000 + 0.5);
         th->SetLineColor(pm + 1);
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            if (EventList[ev].getCut())
               continue;
            if (!getWaveform(ev, pm)->getCut())
               th->Fill(getWaveform(ev, pm)->getMaxAmplitude());
         }
         th->SetLineWidth(2.0);
         return th;
      }
      TH1F *TH1Amp2(int pm)
      {
         TH1F *th = new TH1F(Form("Run%i_OpChannel%i_Amp", Number, pm), Form("Run%i_OpChannel%i_%s_Amp;Amplitude (ADC);Events", Number, pm, PMT_SN[pm].c_str()), ADCDynamicRange * 10, 0.5, ADCDynamicRange + 0.5);
         th->SetLineColor(pm + 1);
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            if (EventList[ev].getCut())
               continue;
            if (!getWaveform(ev, pm)->getCut())
               th->Fill(PedestalPeak[pm] - getWaveform(ev, pm)->getMinimumSample());
         }
         th->SetLineWidth(2.0);
         return th;
      }
      TH1F *TH1AmpRange2(int pm, string units = "ADC")
      {
         cout << "AmpRange2" << endl;

         double max = ADCDynamicRange + 0.5;
         double min = 0.5;
         double Nbins = max - min;

         if (units == "PE")
         {
            max = max / PMT_SPEAmp[pm];
            min = min / PMT_SPEAmp[pm];
         }

         TH1F *th = new TH1F(Form("Run%i_OpChannel%i_Amp", Number, pm), Form("Run%i_OpChannel%i_%s_Amp;Amplitude (ADC);Events", Number, pm, PMT_SN[pm].c_str()), Nbins, min, max);
         th->SetLineColor(pm + 1);
         th->SetLineWidth(2.0);

         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            if (EventList[ev].getCut())
               continue;
            if (!getWaveform(ev, pm)->getCut())
            {
               if (units == "ADC")
                  th->Fill(PedestalPeak[pm] - getWaveform(ev, pm)->getMinimumSampleRange());
               if (units == "PE")
                  th->Fill((PedestalPeak[pm] - getWaveform(ev, pm)->getMinimumSampleRange()) / PMT_SPEAmp[pm]);
            }
         }
         return th;
      }
      TH1F *TH1AmpRange(int pm, string units = "ADC", std::vector<float> *limits = NULL)
      {
         float min = 0.5;
         float max = ADCDynamicRange + 0.5;
         if (limits)
         {
            min = limits->at(0);
            max = limits->at(1);
         }
         if (!limits && units == "PE")
         {
            min = min / PMT_SPEAmp[pm];
            max = max / PMT_SPEAmp[pm];
         }
         //       if(!limits && units=="ADC") {min*=Sampling*2/ADCDynamicRange/50.0/1.602e-19/PMT_Gains[pm]; max*=Sampling*2/ADCDynamicRange/50.0/1.602e-19/PMT_Gains[pm];}

         TH1F *th = new TH1F(Form("Run%i_OpChannel%i_Amp", Number, pm), Form("Run%i_OpChannel%i_%s_Amp;Amplitude (%s);Events", Number, pm, PMT_SN[pm].c_str(), units.c_str()), ADCDynamicRange, min, max);
         th->SetLineColor(pm + 1);
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {

            if (EventList[ev].getCut())
               continue;
            if (getWaveform(ev, pm)->getCut())
               continue;
            if (units == "ADC")
               th->Fill(getWaveform(ev, pm)->getMaxAmplitudeRange());
            if (units == "PE")
               th->Fill(getWaveform(ev, pm)->getMaxAmplitudeRange() / PMT_SPEAmp[pm]);
         }
         th->SetLineWidth(2.0);
         return th;
      }
      TH2F *TH2AmpVsCharge(int pm, string option = "", string units = "PE")
      {
         double ChargeFactor;
         if (units == "pC")
            ChargeFactor = 1.e12;
         if (units == "PE")
            ChargeFactor = 1.0 / 1.602e-19 / PMT_Gains[pm];

         TH1F *thaux = new TH1F(Form("Run%i_OpChannel%i_aux", Number, pm), Form("Run%i_OpChannel%i_%s_%.0fV_aux;Charge (%s);Events", Number, pm, PMT_SN[pm].c_str(), PMT_Voltages[pm], units.c_str()), 5000, -0, 0);
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            //           double correction=(PedestalPeak[pm]-EventList[pm][i].getPedestalMean())*(2.0/ADCDynamicRange)*Sampling*NSamples/50;
            if (getEvent(ev)->getCut())
               continue;
            if (getWaveform(ev, pm)->getCut())
               continue;
            if (option == "")
               thaux->Fill(getWaveform(ev, pm)->getCharge() * ChargeFactor);
            if (option == "Range")
               thaux->Fill(getWaveform(ev, pm)->getChargeRange() * ChargeFactor);
            if (option == "MaxPeak")
               thaux->Fill(getWaveform(ev, pm)->getChargeMaxPeak() * ChargeFactor);
            if (option == "MaxPeakRange")
               thaux->Fill(getWaveform(ev, pm)->getChargeMaxPeakRange() * ChargeFactor);
            if (option == "Q1MaxPeakRange")
               thaux->Fill(getWaveform(ev, pm)->getQ1MaxPeakRange() * ChargeFactor);
            if (option == "Q2MaxPeakRange")
               thaux->Fill(getWaveform(ev, pm)->getQ2MaxPeakRange() * ChargeFactor);
            if (option == "Q3MaxPeakRange")
               thaux->Fill(getWaveform(ev, pm)->getQ3MaxPeakRange() * ChargeFactor);
         }
         Double_t xq2[2]; // position where to compute the quantiles in [0,1]
         Double_t yq2[2]; // array to contain the quantiles
         xq2[0] = 0.0005;
         xq2[1] = 0.9995;
         thaux->GetQuantiles(2, yq2, xq2);

         TH2F *th = new TH2F(Form("Run%i_OpChannel%i_AmpVsCharge_%s", Number, pm, option.c_str()), Form("Run%i_OpChannel%i_%s_%.0fV_AmpVsCharge_%s;Charge (%s);Amplitude (ADC)", Number, pm, PMT_SN[pm].c_str(), PMT_Voltages[pm], option.c_str(), units.c_str()), 300, yq2[0] - 0.15 * (yq2[1] - yq2[0]), yq2[1] + 0.15 * (yq2[1] - yq2[0]), 1000, 0, ADCDynamicRange);
         th->SetLineColor(pm + 1);
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            //           double correction=(PedestalPeak[pm]-EventList[pm][i].getPedestalMean())*(2.0/ADCDynamicRange)*Sampling*NSamples/50;
            if (getEvent(ev)->getCut())
               continue;
            if (getWaveform(ev, pm)->getCut())
               continue;
            if (option == "")
               th->Fill(getWaveform(ev, pm)->getCharge() * ChargeFactor, getWaveform(ev, pm)->getMaxAmplitude());
            if (option == "Range")
               th->Fill(getWaveform(ev, pm)->getChargeRange() * ChargeFactor, getWaveform(ev, pm)->getMaxAmplitudeRange());
            if (option == "MaxPeak")
               th->Fill(getWaveform(ev, pm)->getChargeMaxPeak() * ChargeFactor, getWaveform(ev, pm)->getMaxAmplitude());
            if (option == "MaxPeakRange")
               th->Fill(getWaveform(ev, pm)->getChargeMaxPeakRange() * ChargeFactor, getWaveform(ev, pm)->getMaxAmplitudeRange());
            if (option == "Q1MaxPeakRange")
               th->Fill(getWaveform(ev, pm)->getQ1MaxPeakRange() * ChargeFactor, getWaveform(ev, pm)->getMaxAmplitudeRange());
            if (option == "Q2MaxPeakRange")
               th->Fill(getWaveform(ev, pm)->getQ2MaxPeakRange() * ChargeFactor, getWaveform(ev, pm)->getMaxAmplitudeRange());
            if (option == "Q3MaxPeakRange")
               th->Fill(getWaveform(ev, pm)->getQ3MaxPeakRange() * ChargeFactor, getWaveform(ev, pm)->getMaxAmplitudeRange());
         }
         th->Draw();
         lets_pause();
         delete thaux;
         return th;
      }

      TH2F *TH2Range34vsRange56(int pm, string units = "pC")
      {
         double ChargeFactor;
         if (units == "pC")
            ChargeFactor = 1.e12;
         //~if(units=="PE") ChargeFactor = 1.0/1.602e-19/PMT_Gains[pm];
         //~TH2F *th = new TH2F(Form("Run%i_OpChannel%i_AmpVsCharge_%s",Number,pm,"ChargeCustomRanges",Form("Run%i_OpChannel%i_%s_%.0fV_AmpVsCharge_%s;Charge (%s);Amplitude (ADC)",Number,pm,PMT_SN[pm].c_str(),PMT_Voltages[pm],"ChargeCustomRanges",units.c_str()),300, 0,1 ,1000,0,1); th->SetLineColor(pm+1);
         TH2F *th = new TH2F(Form("Run%i_OpChannel%i_AmpVsCharge_%s", Number, pm, "ChargeCustomRanges"), Form("Run%i_OpChannel%i_%s_%.0fV_AmpVsCharge_%s;Charge in 1-4 #mus (%s);Charge in 80 ns (pC)", Number, pm, PMT_SN[pm].c_str(), PMT_Voltages[pm], "ChargeCustomRanges", units.c_str()), 1000, 0, 3000, 1000, 0, 3000);
         th->SetLineColor(pm + 1);
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            th->Fill(getWaveform(ev, pm)->getChargeRange56() * ChargeFactor, getWaveform(ev, pm)->getChargeRange34() * ChargeFactor);
            //~cout<<getWaveform(ev,pm)->getChargeRange34()*ChargeFactor<<endl;
            //~cout<<getWaveform(ev,pm)->getChargeRange56()*ChargeFactor<<endl;
            //~lets_pause();
         }
         th->Draw();
         lets_pause();
         return th;
      }

      TH2F *TH2F90(int pm, string units = "pC")
      {
         double ChargeFactor;
         if (units == "pC")
            ChargeFactor = 1.e12;
         //~if(units=="PE") ChargeFactor = 1.0/1.602e-19/PMT_Gains[pm];
         //~TH2F *th = new TH2F(Form("Run%i_OpChannel%i_AmpVsCharge_%s",Number,pm,"ChargeCustomRanges",Form("Run%i_OpChannel%i_%s_%.0fV_AmpVsCharge_%s;Charge (%s);Amplitude (ADC)",Number,pm,PMT_SN[pm].c_str(),PMT_Voltages[pm],"ChargeCustomRanges",units.c_str()),300, 0,1 ,1000,0,1); th->SetLineColor(pm+1);
         TH2F *th = new TH2F(Form("Run%i_OpChannel%i_AmpVsCharge_%s", Number, pm, "ChargeCustomRanges"), Form("Run%i_OpChannel%i_%s_%.0fV_AmpVsCharge_%s;Charge in 1-4 #mus (%s); F90", Number, pm, PMT_SN[pm].c_str(), PMT_Voltages[pm], "ChargeCustomRanges", units.c_str()), 1000, 0, 3000, 1000, 0, 1);
         th->SetLineColor(pm + 1);
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            th->Fill(getWaveform(ev, pm)->getChargeRange56() * ChargeFactor, getWaveform(ev, pm)->getChargeRange34() / getWaveform(ev, pm)->getChargeRange56());
            //~cout<<getWaveform(ev,pm)->getChargeRange34()*ChargeFactor<<endl;
            //~cout<<getWaveform(ev,pm)->getChargeRange56()*ChargeFactor<<endl;
            //~lets_pause();
         }
         th->Draw();
         lets_pause();
         return th;
      }

      TH1F *TH1Charge(int pm, string option = "", string units = "PE", double rangemin = 0, double rangemax = 0)
      {
         double ChargeFactor = 1.0;
         if (units == "pC")
            ChargeFactor = 1.e12;
         if (units == "PE")
            ChargeFactor = 1.0 / 1.602e-19 / PMT_Gains[pm];

         if (rangemin == 0 && rangemax == 0)
         {
            TH1F *thaux = new TH1F(Form("Run%i_OpChannel%i_aux", Number, pm), Form("Run%i_OpChannel%i_%s_%.0fV_aux;Charge (%s);Events", Number, pm, PMT_SN[pm].c_str(), PMT_Voltages[pm], units.c_str()), 5000, 0, 0);
            for (unsigned int ev = 0; ev < NEvents; ev++)
            {
               //           double correction=(PedestalPeak[pm]-EventList[pm][i].getPedestalMean())*(2.0/ADCDynamicRange)*Sampling*NSamples/50;
               if (getEvent(ev)->getCut())
                  continue;
               if (getWaveform(ev, pm)->getCut())
                  continue;
               if (option == "")
                  thaux->Fill(getWaveform(ev, pm)->getCharge() * ChargeFactor);
               if (option == "Range")
                  thaux->Fill(getWaveform(ev, pm)->getChargeRange() * ChargeFactor);
               if (option == "MaxPeak")
                  thaux->Fill(getWaveform(ev, pm)->getChargeMaxPeak() * ChargeFactor);
               if (option == "MaxPeakRange")
                  thaux->Fill(getWaveform(ev, pm)->getChargeMaxPeakRange() * ChargeFactor);
               if (option == "Q1MaxPeakRange")
                  thaux->Fill(getWaveform(ev, pm)->getQ1MaxPeakRange() * ChargeFactor);
               if (option == "Q2MaxPeakRange")
                  thaux->Fill(getWaveform(ev, pm)->getQ2MaxPeakRange() * ChargeFactor);
               if (option == "Q3MaxPeakRange")
                  thaux->Fill(getWaveform(ev, pm)->getQ3MaxPeakRange() * ChargeFactor);
            }
            // thaux->Draw("HIST"); lets_pause();
            Double_t xq2[2]; // position where to compute the quantiles in [0,1]
            Double_t yq2[2]; // array to contain the quantiles
            xq2[0] = 0.001;
            xq2[1] = 0.999;
            thaux->GetQuantiles(2, yq2, xq2);
            rangemin = yq2[0] - 0.15 * yq2[1];
            rangemax = 1.2 * yq2[1];
            delete thaux;
         }

         TH1F *th = new TH1F(Form("Run%i_OpChannel%i_Charge_%s", Number, pm, option.c_str()), Form("Run%i_OpChannel%i_%s_%.0fV_Charge_%s;Charge (%s);Counts", Number, pm, PMT_SN[pm].c_str(), PMT_Voltages[pm], option.c_str(), units.c_str()), 600, rangemin, rangemax);
         th->SetLineColor(pm + 1);
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            //           double correction=(PedestalPeak[pm]-EventList[pm][i].getPedestalMean())*(2.0/ADCDynamicRange)*Sampling*NSamples/50;
            if (getEvent(ev)->getCut())
               continue;
            if (getWaveform(ev, pm)->getCut())
               continue;
            if (option == "")
               th->Fill(getWaveform(ev, pm)->getCharge() * ChargeFactor);
            if (option == "Range")
               th->Fill(getWaveform(ev, pm)->getChargeRange() * ChargeFactor);
            if (option == "MaxPeak")
               th->Fill(getWaveform(ev, pm)->getChargeMaxPeak() * ChargeFactor);
            if (option == "MaxPeakRange")
               th->Fill(getWaveform(ev, pm)->getChargeMaxPeakRange() * ChargeFactor);
            if (option == "Q1MaxPeakRange")
               th->Fill(getWaveform(ev, pm)->getQ1MaxPeakRange() * ChargeFactor);
            if (option == "Q2MaxPeakRange")
               th->Fill(getWaveform(ev, pm)->getQ2MaxPeakRange() * ChargeFactor);
            if (option == "Q3MaxPeakRange")
               th->Fill(getWaveform(ev, pm)->getQ3MaxPeakRange() * ChargeFactor);
         }
         return th;
      }

      TH2F *TH2AmpVsCharge_FixedPedestal(int pm)
      {
         TH1F *thaux = new TH1F(Form("Run%i_OpChannel%i_Charge_aux", Number, pm), Form("Run%i_OpChannel%i_%s_%.0fV_pedSTD<%.1fADC_aux;Charge (PE);Events", Number, pm, PMT_SN[pm].c_str(), PMT_Voltages[pm], PedestalSTDCUT[pm]), 5000, -2e6, 1e6);
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            //           double correction=(PedestalPeak[pm]-EventList[pm][i].getPedestalMean())*(2.0/ADCDynamicRange)*Sampling*NSamples/50;
            if (!getWaveform(ev, pm)->getCut())
               thaux->Fill(getWaveform(ev, pm)->getSampling() * (getWaveform(ev, pm)->getNSamples() * PedestalPeak[pm] - getWaveform(ev, pm)->getSumOfSamples()) * 2 / ADCDynamicRange / 50 / 1.602e-19 / PMT_Gains[pm]);
         }
         Double_t xq2[2]; // position where to compute the quantiles in [0,1]
         Double_t yq2[2]; // array to contain the quantiles
         xq2[0] = 0.001;
         xq2[1] = 0.999;
         thaux->GetQuantiles(2, yq2, xq2);

         TH2F *th = new TH2F(Form("Run%i_OpChannel%i_ChargePEVsAmp", Number, pm), Form("Run%i_OpChannel%i_%s_%.0fV_pedSTD<%.1fADC;Charge (PE);Amplitude (ADC)", Number, pm, PMT_SN[pm].c_str(), PMT_Voltages[pm], PedestalSTDCUT[pm]), 300, 2 * yq2[0], 1.2 * yq2[1], 1000, 0, ADCDynamicRange);
         th->SetLineColor(pm + 1);
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            //           double correction=(PedestalPeak[pm]-EventList[pm][i].getPedestalMean())*(2.0/ADCDynamicRange)*Sampling*NSamples/50;
            double totalcharge = getWaveform(ev, pm)->getSampling() * (getWaveform(ev, pm)->getNSamples() * PedestalPeak[pm] - getWaveform(ev, pm)->getSumOfSamples()) * 2 / ADCDynamicRange / 50 / 1.602e-19 / PMT_Gains[pm];
            if (!getWaveform(ev, pm)->getCut())
               th->Fill(totalcharge, PedestalPeak[pm] - getWaveform(ev, pm)->getMinimumSample());
            /*
                     event.GetEntry(midas_chain,ev);
                     TCanvas *c = new TCanvas("PMT");
                     DrawEvent(*getWaveform(ev,pm),pm,ev,c);
                       lets_pause();*/
         }
         delete thaux;
         return th;
      }
      TH1F *TH1_AmplitudeRatio(int pm1, int pm2)
      {
         /*
        TH1F *thaux = new TH1F(Form("Run%i_OpChannel%i_Charge_aux",Number,pm1),Form("Run%i_OpChannel%i_%s_%.0fV_pedSTD<%.1fADC_aux;Charge (PE);Events",Number,pm1,PMT_SN[pm1].c_str(),PMT_Voltages[pm1],PedestalSTDCUT[pm1]),5000,-2e6,1e6);
        for(int ev=0;ev<NEvents;ev++)
        {
            if(EventList[ev].getCut()) continue;
            if(getWaveform(ev,pm1)->getCut()) continue;
            if(getWaveform(ev,pm2)->getCut()) continue;
            thaux->Fill((PedestalPeak[pm1]-getWaveform(ev,pm1)->getMinimumSample()) / (PedestalPeak[pm2]-getWaveform(ev,pm2)->getMinimumSample()));
        }
        Double_t xq2[2];  // position where to compute the quantiles in [0,1]
        Double_t yq2[2];  // array to contain the quantiles
        xq2[0] = 0.001;
        xq2[1] = 0.999;
        thaux->GetQuantiles(2,yq2,xq2);
 */
         TH1F *th = new TH1F(Form("Run%i_AmplitudeRatio", Number), Form("Run%i_AmplitudeRatio_channel%i_over_channel%i;Ratio;Counts", Number, pm1, pm2), 500, -10, 100); // th->SetLineColor(pm+1);
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {

            if (EventList[ev].getCut())
               continue;
            if (getWaveform(ev, pm1)->getCut())
               continue;
            if (getWaveform(ev, pm2)->getCut())
               continue;
            th->Fill((PedestalPeak[pm1] - getWaveform(ev, pm1)->getMinimumSample()) / (PedestalPeak[pm2] - getWaveform(ev, pm2)->getMinimumSample()));
         }
         // delete thaux;
         return th;
      }
      TH2F *TH2_Amplitudes(int pm1, int pm2)
      {

         TH1F *thaux = new TH1F("Aux1", "", 5000, -2e6, 1e6);
         TH1F *thaux2 = new TH1F("Aux2", "", 5000, -2e6, 1e6);
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            if (EventList[ev].getCut())
               continue;
            if (getWaveform(ev, pm1)->getCut())
               continue;
            if (getWaveform(ev, pm2)->getCut())
               continue;
            thaux->Fill((PedestalPeak[pm1] - getWaveform(ev, pm1)->getMinimumSample()));
            thaux2->Fill((PedestalPeak[pm2] - getWaveform(ev, pm2)->getMinimumSample()));
         }
         Double_t xq[2]; // position where to compute the quantiles in [0,1]
         Double_t yq[2]; // array to contain the quantiles
         xq[0] = 0.001;
         xq[1] = 0.999;
         thaux->GetQuantiles(2, yq, xq);
         Double_t xq2[2]; // position where to compute the quantiles in [0,1]
         Double_t yq2[2]; // array to contain the quantiles
         xq2[0] = 0.001;
         xq2[1] = 0.999;
         thaux->GetQuantiles(2, yq2, xq2);

         TH2F *th = new TH2F(Form("Run%i_Amplitudes", Number), Form("Run%i_AmplitudesCorrleation;Ch%i Amplitude (ADC);Ch%i Amplitude (ADC)", Number, pm1, pm2), 1000, 2 * yq[0], 1.2 * yq[1], 1000, 2 * yq2[0], 1.2 * yq2[1]); // th->SetLineColor(pm+1);
         for (unsigned int ev = 0; ev < NEvents; ev++)
         {

            if (EventList[ev].getCut())
               continue;
            if (getWaveform(ev, pm1)->getCut())
               continue;
            if (getWaveform(ev, pm2)->getCut())
               continue;
            th->Fill((PedestalPeak[pm1] - getWaveform(ev, pm1)->getMinimumSample()), (PedestalPeak[pm2] - getWaveform(ev, pm2)->getMinimumSample()));
         }
         delete thaux;
         delete thaux2;
         return th;
      }

      void PlotFreqMinimum(string ofile)
      {
         TGraphErrors *th[64];
         for (auto pm : adcchannels)
         {
            th[pm] = TGraphFreqMinimum(pm);
         }
         TCanvas *c = new TCanvas("c");
         TMultiGraph *m = new TMultiGraph("m", "m");
         for (auto pm : adcchannels)
            m->Add(th[pm]);
         c->cd();
         m->Draw("AP");
         gPad->BuildLegend();
         m->GetYaxis()->SetTitle("Frequency (Hz)");
         m->GetXaxis()->SetTitle("Amplitude threshold (ADC)");
         lets_pause();
         TFile *of = new TFile(ofile.c_str(), "RECREATE");
         of->cd();
         for (auto pm : adcchannels)
         {
            th[pm]->Write();
         }
         of->Close();
         delete of;
         for (auto pm : adcchannels)
            delete th[pm];
      }

      void PlotFreqAmp(string ofile)
      {
         TGraphErrors *th[64];
         for (auto pm : adcchannels)
         {
            th[pm] = TGraphFreqAmp(pm);
         }
         TCanvas *c = new TCanvas("c");
         TMultiGraph *m = new TMultiGraph("m", "m");
         for (auto pm : adcchannels)
            m->Add(th[pm]);
         c->cd();
         m->Draw("AP");
         gPad->BuildLegend();
         m->GetYaxis()->SetTitle("Frequency (Hz)");
         m->GetXaxis()->SetTitle("Amplitude threshold (ADC)");
         lets_pause();
         TFile *of = new TFile(ofile.c_str(), "RECREATE");
         of->cd();
         for (auto pm : adcchannels)
         {
            th[pm]->Write();
         }
         of->Close();
         delete of;
         for (auto pm : adcchannels)
            delete th[pm];
      }
      void PlotFreqAmpInverse(string ofile)
      {
         TGraphErrors *th[64];
         for (auto pm : adcchannels)
         {
            th[pm] = TGraphFreqAmpInverse(pm);
         }
         TCanvas *c = new TCanvas("c");
         TMultiGraph *m = new TMultiGraph("m", "m");
         for (auto pm : adcchannels)
            m->Add(th[pm]);
         c->cd();
         m->Draw("AP");
         gPad->BuildLegend();
         m->GetYaxis()->SetTitle("Frequency (Hz)");
         m->GetXaxis()->SetTitle("Amplitude (ADC)");
         lets_pause();
         TFile *of = new TFile(ofile.c_str(), "RECREATE");
         of->cd();
         for (auto pm : adcchannels)
         {
            th[pm]->Write();
         }
         of->Close();
         delete of;
         for (auto pm : adcchannels)
            delete th[pm];
      }

      void PlotFreqCharge()
      {
         TGraphErrors *th[64];
         for (auto pm : adcchannels)
         {
            th[pm] = TGraphFreqCharge(pm);
         }
         TCanvas *c = new TCanvas("c");
         TMultiGraph *m = new TMultiGraph("m", "m");
         for (auto pm : adcchannels)
            m->Add(th[pm]);
         c->cd();
         m->Draw("ALP");
         gPad->BuildLegend();
         m->GetYaxis()->SetTitle("Frequency (Hz)");
         m->GetXaxis()->SetTitle("Charge threshold (PEs)");
         lets_pause();
      }
      bool HasTPB(string name)
      {
         if (name == "FA0114")
            return true;
         if (name == "FA0110")
            return true;
         if (name == "FA0132")
            return true;
         if (name == "FA0139")
            return true;
         if (name == "FA0107")
            return true;
         if (name == "FA0137")
            return true;
         return false;
      }
      bool HasTPB(int pm)
      {
         return HasTPB(PMT_SN[pm]);
      }
      TGraphErrors *TGraphFreqAmp(int pm) // only have sense using random trigger!
      {
         int Npoints = 24;
         std::vector<double> Thres;
         Thres.resize(Npoints);
         Thres[0] = 2.6;
         double factor = 1.45;
         std::vector<double> Counts;
         Counts.resize(Npoints, 0.0);
         int eventcounter = 0;
         std::vector<double> ex;
         ex.resize(Npoints, 0.0);
         std::vector<double> ey;
         ey.resize(Npoints, 0.0);
         for (unsigned int j = 1; j < Thres.size(); j++)
         {
            Thres[j] = Thres[j - 1] * factor;
         }
         for (unsigned int ev = 0; ev < NEvents; ev++)
            if (!getWaveform(ev, pm)->getCut())
            {
               for (unsigned int j = 0; j < Thres.size(); j++)
                  if (getWaveform(ev, pm)->getMaxAmplitude() > Thres[j])
                  {
                     Counts[j]++;
                  }
               eventcounter++;
            }
         for (unsigned int j = 0; j < Thres.size(); j++)
         {
            Counts[j] = Counts[j] / (eventcounter * getWaveform(0, pm)->getSampling() * getWaveform(0, pm)->getNSamples());
            ex[j] = 1.0;
            ey[j] = TMath::Sqrt(Counts[j]) / (eventcounter * getWaveform(0, pm)->getSampling() * getWaveform(0, pm)->getNSamples());
            std::cout << "Threshold: " << Thres[j] << " , Freq: " << Counts[j] << " - error (" << ex[j] << " " << ey[j] << ")" << std::endl;
         }
         TGraphErrors *tg = new TGraphErrors(Thres.size(), &(Thres[0]), &(Counts[0]), &(ex[0]), &(ey[0]));
         tg->SetLineColor(pm + 1);
         //       tg->SetTitle(Form("Run%i_ch%i_%s_%.0fV;Frequency (Hz);Maximum amplitud per %.1fus event(ADC)",Number,pm,PMT_SN[pm].c_str(),PMT_Voltages[pm],getWaveform(0,pm)->getSampling()*getWaveform(0,pm)->getNSamples()*1.e6));
         if (HasTPB(PMT_SN[pm]))
            tg->SetTitle(Form("ch%i_%s_%.0fV_TPB", pm, PMT_SN[pm].c_str(), PMT_Voltages[pm]));
         else
            tg->SetTitle(Form("ch%i_%s_%.0fV", pm, PMT_SN[pm].c_str(), PMT_Voltages[pm]));
         tg->SetName(Form("Run%i_ch%i_%s", Number, pm, PMT_SN[pm].c_str()));
         tg->SetMarkerStyle(20);
         tg->SetMarkerColor(pm + 1);
         return tg;
      }
      TGraphErrors *TGraphFreqMinimum(int pm) // only have sense using random trigger!
      {
         int Npoints = 24;
         std::vector<double> Thres;
         Thres.resize(Npoints);
         Thres[0] = 2.6;
         double factor = 1.45;
         std::vector<double> Counts;
         Counts.resize(Npoints, 0.0);
         int eventcounter = 0;
         std::vector<double> ex;
         ex.resize(Npoints, 0.0);
         std::vector<double> ey;
         ey.resize(Npoints, 0.0);
         for (unsigned int j = 1; j < Thres.size(); j++)
         {
            Thres[j] = Thres[j - 1] * factor;
         }
         for (unsigned int ev = 0; ev < NEvents; ev++)
            if (!getWaveform(ev, pm)->getCut())
            {
               for (unsigned int j = 0; j < Thres.size(); j++)
                  if (PedestalPeak[pm] - getWaveform(ev, pm)->getMinimumSample() > Thres[j])
                  {
                     Counts[j]++;
                  }
               eventcounter++;
            }
         for (unsigned int j = 0; j < Thres.size(); j++)
         {
            Counts[j] = Counts[j] / (eventcounter * getWaveform(0, pm)->getSampling() * getWaveform(0, pm)->getNSamples());
            ex[j] = 1.0;
            ey[j] = TMath::Sqrt(Counts[j]) / (eventcounter * getWaveform(0, pm)->getSampling() * getWaveform(0, pm)->getNSamples());
            std::cout << "Threshold: " << Thres[j] << " , Freq: " << Counts[j] << " - error (" << ex[j] << " " << ey[j] << ")" << std::endl;
         }
         TGraphErrors *tg = new TGraphErrors(Thres.size(), &(Thres[0]), &(Counts[0]), &(ex[0]), &(ey[0]));
         tg->SetLineColor(pm + 1);
         //       tg->SetTitle(Form("Run%i_ch%i_%s_%.0fV;Frequency (Hz);Maximum amplitud per %.1fus event(ADC)",Number,pm,PMT_SN[pm].c_str(),PMT_Voltages[pm],getWaveform(0,pm)->getSampling()*getWaveform(0,pm)->getNSamples()*1.e6));
         if (HasTPB(PMT_SN[pm]))
            tg->SetTitle(Form("ch%i_%s_%.0fV_TPB", pm, PMT_SN[pm].c_str(), PMT_Voltages[pm]));
         else
            tg->SetTitle(Form("ch%i_%s_%.0fV", pm, PMT_SN[pm].c_str(), PMT_Voltages[pm]));
         tg->SetName(Form("Run%i_ch%i_%s", Number, pm, PMT_SN[pm].c_str()));
         tg->SetMarkerStyle(20);
         tg->SetMarkerColor(pm + 1);
         return tg;
      }
      TH1F *TH1FreqMinimum(int pm) // only have sense using random trigger!
      {
         const Int_t Npoints = 24;
         std::vector<double> Thres;
         Thres.resize(Npoints);
         Thres[0] = 2.6;
         double factor = 1.45;
         std::vector<double> Counts;
         Counts.resize(Npoints, 0.0);
         int eventcounter = 0;
         std::vector<double> ex;
         ex.resize(Npoints, 0.0);
         std::vector<double> ey;
         ey.resize(Npoints, 0.0);
         Double_t edges[Npoints + 1];
         for (unsigned int j = 1; j < Thres.size(); j++)
         {
            Thres[j] = Thres[j - 1] * factor;
         }
         for (unsigned int j = 1; j < Npoints + 1; j++)
         {
            edges[j] = Thres[j - 1] * factor;
         }
         edges[0] = 2.6;

         TH1F *h = new TH1F(Form("ch%i_%s_%.0fV", pm, PMT_SN[pm].c_str(), PMT_Voltages[pm]), Form("ch%i_%s_%.0fV;Amplitude (ADC); Frequency (Hz)", pm, PMT_SN[pm].c_str(), PMT_Voltages[pm]), Npoints, edges);
         h->Sumw2();

         for (unsigned int ev = 0; ev < NEvents; ev++)
            if (!getWaveform(ev, pm)->getCut())
            {
               h->Fill(PedestalPeak[pm] - getWaveform(ev, pm)->getMinimumSample());
               eventcounter++;
            }
         h->Scale(1.0 / (eventcounter * getWaveform(0, pm)->getSampling() * getWaveform(0, pm)->getNSamples()));

         if (HasTPB(PMT_SN[pm]))
            h->SetTitle(Form("ch%i_%s_%.0fV_TPB", pm, PMT_SN[pm].c_str(), PMT_Voltages[pm]));
         else
            h->SetTitle(Form("ch%i_%s_%.0fV", pm, PMT_SN[pm].c_str(), PMT_Voltages[pm]));
         h->SetName(Form("Run%i_ch%i_%s", Number, pm, PMT_SN[pm].c_str()));
         return h;
      }

      TH1F *TH1Freq(int pm) // only have sense using random trigger!
      {
         const Int_t Npoints = 24;
         std::vector<double> Thres;
         Thres.resize(Npoints);
         Thres[0] = 2.6;
         double factor = 1.45;
         std::vector<double> Counts;
         Counts.resize(Npoints, 0.0);
         int eventcounter = 0;
         std::vector<double> ex;
         ex.resize(Npoints, 0.0);
         std::vector<double> ey;
         ey.resize(Npoints, 0.0);
         Double_t edges[Npoints + 1];
         for (unsigned int j = 1; j < Thres.size(); j++)
         {
            Thres[j] = Thres[j - 1] * factor;
         }
         for (unsigned int j = 1; j < Npoints + 1; j++)
         {
            edges[j] = Thres[j - 1] * factor;
         }
         edges[0] = 2.6;

         TH1F *h = new TH1F(Form("ch%i_%s_%.0fV", pm, PMT_SN[pm].c_str(), PMT_Voltages[pm]), Form("ch%i_%s_%.0fV;Amplitude (ADC); Frequency (Hz)", pm, PMT_SN[pm].c_str(), PMT_Voltages[pm]), Npoints, edges);
         h->Sumw2();

         for (unsigned int ev = 0; ev < NEvents; ev++)
            if (!getWaveform(ev, pm)->getCut())
            {
               h->Fill(getWaveform(ev, pm)->getMaxAmplitude());
               eventcounter++;
            }
         h->Scale(1.0 / (eventcounter * getWaveform(0, pm)->getSampling() * getWaveform(0, pm)->getNSamples()));

         if (HasTPB(PMT_SN[pm]))
            h->SetTitle(Form("ch%i_%s_%.0fV_TPB", pm, PMT_SN[pm].c_str(), PMT_Voltages[pm]));
         else
            h->SetTitle(Form("ch%i_%s_%.0fV", pm, PMT_SN[pm].c_str(), PMT_Voltages[pm]));
         h->SetName(Form("Run%i_ch%i_%s", Number, pm, PMT_SN[pm].c_str()));
         return h;
      }

      TGraphErrors *TGraphFreqAmpInverse(int pm) // only have sense using random trigger!
      {
         int Npoints = 24;
         std::vector<double> Thres;
         Thres.resize(Npoints);
         Thres[0] = 2.6;
         double factor = 1.45;
         std::vector<double> Counts;
         Counts.resize(Npoints, 0.0);
         int eventcounter = 0;
         std::vector<double> ex;
         ex.resize(Npoints, 0.0);
         std::vector<double> ey;
         ey.resize(Npoints, 0.0);
         for (unsigned int j = 1; j < Thres.size(); j++)
         {
            Thres[j] = Thres[j - 1] * factor;
         }
         for (unsigned int ev = 0; ev < NEvents; ev++)
            if (!getWaveform(ev, pm)->getCut())
            {
               for (unsigned int j = 0; j < Thres.size(); j++)
                  if (getWaveform(ev, pm)->getMaxAmplitude() < Thres[j])
                  {
                     Counts[j]++;
                  }
               eventcounter++;
            }
         for (unsigned int j = 0; j < Thres.size(); j++)
         {
            Counts[j] = Counts[j] / (eventcounter * getWaveform(0, pm)->getSampling() * getWaveform(0, pm)->getNSamples());
            ex[j] = 1.0;
            ey[j] = TMath::Sqrt(Counts[j]) / (eventcounter * getWaveform(0, pm)->getSampling() * getWaveform(0, pm)->getNSamples());
            std::cout << "Threshold: " << Thres[j] << " , Freq: " << Counts[j] << " - error (" << ex[j] << " " << ey[j] << ")" << std::endl;
         }
         TGraphErrors *tg = new TGraphErrors(Thres.size(), &(Thres[0]), &(Counts[0]), &(ex[0]), &(ey[0]));
         tg->SetLineColor(pm + 1);
         //       tg->SetTitle(Form("Run%i_ch%i_%s_%.0fV;Frequency (Hz);Maximum amplitud per %.1fus event(ADC)",Number,pm,PMT_SN[pm].c_str(),PMT_Voltages[pm],getWaveform(0,pm)->getSampling()*getWaveform(0,pm)->getNSamples()*1.e6));
         if (HasTPB(PMT_SN[pm]))
            tg->SetTitle(Form("ch%i_%s_%.0fV_TPB", pm, PMT_SN[pm].c_str(), PMT_Voltages[pm]));
         else
            tg->SetTitle(Form("ch%i_%s_%.0fV", pm, PMT_SN[pm].c_str(), PMT_Voltages[pm]));
         tg->SetName(Form("Run%i_ch%i_%s", Number, pm, PMT_SN[pm].c_str()));
         tg->SetMarkerStyle(20);
         tg->SetMarkerColor(pm + 1);
         return tg;
      }
      TGraphErrors *TGraphFreqCharge(int pm)
      {
         std::vector<double> Thres;
         Thres.resize(23);
         Thres[0] = 0.8;
         double factor = 1.52;
         std::vector<double> Counts;
         Counts.resize(23, 0.0);
         int eventcounter = 0;
         for (unsigned int j = 1; j < Thres.size(); j++)
         {
            Thres[j] = Thres[j - 1] * factor;
         }
         for (unsigned int ev = 0; ev < NEvents; ev++)
            if (!getWaveform(ev, pm)->getCut())
            {
               for (unsigned int j = 0; j < Thres.size(); j++)
                  if (getWaveform(ev, pm)->getCharge() / 1.602e-19 / PMT_Gains[pm] > Thres[j])
                  {
                     Counts[j]++;
                  }
               eventcounter++;
            }
         for (unsigned int j = 0; j < Thres.size(); j++)
         {
            Counts[j] = Counts[j] / (eventcounter * getWaveform(0, pm)->getSampling() * getWaveform(0, pm)->getNSamples());
            std::cout << Thres[j] << "  " << Counts[j] << std::endl;
         }
         TGraphErrors *tg = new TGraphErrors(Thres.size(), &(Thres[0]), &(Counts[0]));
         tg->SetLineColor(pm + 1);
         //       tg->SetTitle(Form("Run%i_ch%i_%s_%.0fV;Frequency (Hz);Maximum charge per %.1fus event(PE)",Number,pm,PMT_SN[pm].c_str(),PMT_Voltages[pm],getWaveform(0,pm)->getSampling()*getWaveform(0,pm)->getNSamples()*1.e6));
         tg->SetTitle(Form("Run%i_ch%i_%s_%.0fV_Charge", Number, pm, PMT_SN[pm].c_str(), PMT_Voltages[pm]));
         tg->SetMarkerStyle(20);
         tg->SetMarkerColor(pm + 1);
         return tg;
      }
      /*
          void MakeSPESpec()
          {

              TFile fout(Form("Run%08i_SPESpec.root",Number), "RECREATE");
              TCanvas *c = new TCanvas("c");
              c->Divide(6,6);
              TH1F *th[64];
              for(auto i: adcchannels)
              {
                std::cout << "Analysing PMT " << i << std::endl;
                th[i]=TGraphSPE(i);
                fout.cd();
           th[i]->Write();

                c->cd(i+1);
                th[i]->Draw("HIST");
                gPad->SetLogy();c->cd(i+1)->Modified();c->cd(i+1)->Update();
              }

              fout.Close();
              lets_pause();

          }
          TH1F * TGraphSPE(int pm)
          {
             TH1F *th = new TH1F(Form("Run%i_OpChannel%i_SPE",Number,pm),Form("Run%i_OpChannel%i_SPE;Charge (/e/);Counts",Number,pm),300,-2e7,5e7); th->SetLineColor(pm+1);
             MakeChargeSpec(midas_chain,*th,pm,NSamples);
             return th;
          }*/
      void SetPMTs(std::vector<string> names) // Manually set the PTMs of interest
      {
         PMT_SN = names;
         PMT_Voltages.resize(PMT_SN.size(), 0);
         PMT_Gains.resize(PMT_SN.size(), 1); // adcchannels.clear();
                                             //       for(auto pm: adcchannelsALL){if(PMT_SN[pm]!="") {adcchannels.push_back(pm);}}
         cout << PMT_SN.size() << endl;
      }

      void SetVoltages(std::vector<double> volt) // Manually set the PTM voltages
      {
         PMT_Voltages = volt;
      }
      void SetGains(std::vector<double> mg) // Manually set the PTMs gains
      {
         PMT_Gains = mg;
      }
      void SetGains(std::map<string, double> mg) // Manually set the PTMs gains using a pmt map: PMT_SN -> Gain
      {
         for (auto pm : adcchannels)
            PMT_Gains[pm] = mg[PMT_SN[pm]];
      }
      void SetSPEAmps(std::vector<double> mg) // Manually set the PTMs SPE Amplitudes
      {
         PMT_SPEAmp = mg;
      }

      void SelectChannels(std::vector<int> mychannels)
      {
         adcchannels.clear();
         for (auto pm : mychannels)
         {
            adcchannels.push_back(pm);
            cout << "Selected channel " << pm << " " << PMT_SN[pm] << " to perform the analysis." << endl;
         }
      }
      void RemoveChannels(std::vector<int> channels)
      {
         std::vector<int> aux = adcchannels;
         for (auto pm : channels)
         {
            cout << "Removing channel " << pm << " " << PMT_SN[pm] << " to perform the analysis." << endl;
            adcchannels.erase(std::find(adcchannels.begin(), adcchannels.end(), pm));
         }
      }
      void PrintRunReport(string outfile)
      {
         std::cout << "Creating Report... " << std::endl;
         if (EventList.size() == 0)
            Process();
         TFile *ofile = new TFile(outfile.c_str(), "RECREATE");
         TH1RateOfEvents()->Write();
         for (auto pm : adcchannels)
         {
            ofile->cd();
            TH1Pedestal(pm)->Write();
            //         TH1RateOfEventsPerPMT(pm)->Write();
            //         TH1Saturated(pm)->Write();
            //         TH1ChargeInPEUnits(pm)->Write();
            //         TH1AmpRange2(pm)->Write();
            //         TH2AmpVsCharge(pm,"MaxPeakRange")->Write();
            //         TH2AmpVsCharge(pm,"Q1MaxPeakRange")->Write();
            //         TH2AmpVsCharge(pm,"Q2MaxPeakRange")->Write();
            //         TH2AmpVsCharge(pm,"Q3MaxPeakRange")->Write();
            TH1PeakTime(pm)->Write();
            //         TH1PedSTD(pm)->Write();
            //         TH1ScintProf(pm,-1)->Write();
            //         TH1ScintProf2(pm,50000)->Write();
            //         TH1AvWf(pm)->Write();
            //         TH1AvWfALL(pm)->Write();
            //         TH1AvWf(pm)->Write();
            //         TH1ScintProfTriggerCut(pm)->Write();
            // TProfile_ChargeVsTime(pm)->Write();
            // TProfileEvolution(pm,"pedestal")->Write();
         }
         ofile->Close();
         std::cout << "Report dumped to " << outfile << std::endl;
      }
      void Close()
      {
         EventList.clear();
         cout << "Run " << Number << " closed. " << endl;
      }
      void LoopWaveformsPMT(int pm, int ev0 = 0, string BoxOption = "", std::vector<double> *ranges = NULL, int Rebin = 1)
      {
         for (unsigned int ev = ev0; ev < NEvents; ev++)
         {
            if (EventList[ev].getCut())
               continue;
            if (getWaveform(ev, pm)->getCut())
               continue;
            cout << "Event " << ev << endl;
            GetEntry(ev);
            TCanvas *c = new TCanvas("PMT");
            DrawEvent(getWaveform(ev, pm), pm, ev, c, ranges, BoxOption, Rebin);
            lets_pause();
         }
      }
      void DrawEventPMT(int pm, int ev, string BoxOption = "", std::vector<double> *ranges = NULL, int Rebin = 1)
      {
         TCanvas *c = new TCanvas("PMT");
         DrawEvent(getWaveform(ev, pm), pm, ev, c, ranges, BoxOption, Rebin);
         lets_pause();
      }
      void DrawEvent(int ev, std::vector<int> Channels, string BoxOption = "", std::vector<double> *ranges = NULL, int Rebin = 1, bool Overlap = false)
      {
         lets_pause();
         std::cout << "Drawing waveform" << std::endl;
         TH1F *th[64];
         std::vector<double> *ranges2 = NULL;
         TCanvas *c = new TCanvas("c");
         c->DivideSquare(Channels.size(), 0.002, 0.002);
         if (ev > (int)NEventsALL)
         {
            std::cout << "ERROR: You are trying to draw event " << ev << ", but there is only " << NEvents << " . Exiting ... " << std::endl;
            throw std::exception();
         }
         cout << "Looking at " << ev << endl;

         GetEntry(ev);
         int canvas = 0;
         double max = 0;
         double min = 10000;
         for (auto pm : Channels)
         {
            if (min > getWaveform(ev, pm)->getMinimumSample())
               min = getWaveform(ev, pm)->getMinimumSample();
         }
         for (auto pm : Channels)
         {
            if (max < getWaveform(ev, pm)->getMaximumSample())
               max = getWaveform(ev, pm)->getMaximumSample();
         }
         double delta = (max - min) * 0.15;

         cout << "Lets loop on channels" << endl;
         int i = 0;
         for (auto pm : Channels)
         {
            cout << "drawing channel " << pm << endl;
            int color = pm + 1;
            if (!Overlap)
            {
               c->cd(canvas + 1);
               cout << "MIN " << min << ", MAX " << max << " delta " << delta << endl;
               color = 1;
            }
            if (ranges == NULL)
            {
               ranges2 = new std::vector<double>();
               *ranges2 = {0.0, Sampling * NSamples, min - delta, max + delta};
            }
            else
            {
               ranges2 = ranges;
            }
            if (PMT_SN[pm] == "")
            {
               cout << "Entro" << pm << " vs " << NChannels << endl;
               continue;
            }
            else
            {
               std::cout << "Drawing PMT" << pm << " - evt " << ev << " - timestamp " << Form("%.0f %s", event.GetTimeStamp(), TDatime(event.GetTimeStamp()).AsString()) << " " << event.GetTimeStamp_nsec() << "ns" << std::endl;
               TH1D *h = new TH1D();
               // DrawEvent( waveana::Waveform_t(event.GetTime(),event.GetAmp(pm),PedRange,TriggerRange,-(ADCDynamicRange/2.0)*50.0, range1, range2),pm,ev,gPad,ranges,BoxOption);  canvas+=1;}
               lets_pause();
               DrawEvent(getWaveform(ev, pm), pm, ev, gPad, ranges2, BoxOption, Rebin, color);
               canvas += 1;
            }
            if (!Overlap)
            {
               c->cd(i + 1)->Modified();
               c->cd(i + 1)->Update();
            }
            i++;
         }
         lets_pause();
      }

      void LoopWaveforms(int myev = 0, string BoxOption = "", std::vector<double> *ranges = NULL, int Rebin = 1, int CutOption = -1, bool donotloop = false, std::vector<int> *SelectedChannels = NULL, bool Overlap = false)
      {
         /*CutOption:
         -1 (default) : Loop over all events, skipping cut events.
         n (channel #): Loop over all events, skipping cut events, and events with waveform of channel n cut.
         -2 : Loop over all events (not skipping any event), you will see a text pannel on cut waveforms though.
         -3 : Loop over cut events. */
         std::cout << "Drawing waveform" << std::endl;
         TH1F *th[64];
         std::vector<double> *ranges2 = NULL;
         TCanvas *c = new TCanvas("c");
         if (!Overlap)
            c->DivideSquare(adcchannels.size(), 0.002, 0.002); // c->Divide(6,6,0.01,0.01);
         if (myev > (int)NEventsALL)
         {
            std::cout << "ERROR: You are trying to draw event " << myev << ", but there is only " << NEvents << " . Exiting ... " << std::endl;
            throw std::exception();
         }
         int lastev = NEventsALL;
         if (donotloop)
            lastev = myev + 1;
         std::vector<int> Channels = adcchannels;
         if (SelectedChannels)
            Channels = *SelectedChannels;
         for (unsigned int ev = myev; (int)ev < lastev; ev++)
         {
            cout << "Looking at " << ev << endl;

            GetEntry(ev);
            int canvas = 0;
            double max = 0;
            double min = 10000;

            Event_t EVT = ProcessEvent(ev);
            for (auto pm : Channels)
            {
               MyCuts.ApplyCuts(*EVT.getChannel(pm), pm, TH1getWaveform(ev, pm));
            }
            MyCuts.ApplyMultiChannelCuts(EVT);
            cout << "Cuts applied" << endl;

            if (CutOption > -2)
            {
               if (EVT.getCut())
                  continue; // std::cout << "Entro, cut: " << EVT.getChannel(CutOption)->getCut() << std::endl;
               if (CutOption != -1)
                  if (EVT.getChannel(CutOption)->getCut())
                     continue;
            }
            if (CutOption == -3)
            {
               if (!EVT.getCut())
                  continue;
            }
            for (auto pm : Channels)
            {
               if (min > EVT.getChannel(pm)->getMinimumSample())
                  min = EVT.getChannel(pm)->getMinimumSample();
            }
            for (auto pm : Channels)
            {
               if (max < EVT.getChannel(pm)->getMaximumSample())
                  max = EVT.getChannel(pm)->getMaximumSample();
            }
            double delta = (max - min) * 0.15;

            cout << "Lets loop on channels" << endl;

            for (auto pm : Channels)
            {
               cout << "drawing channel " << pm << endl;

               int color = pm + 1;
               if (!Overlap)
               {
                  c->cd(canvas + 1);
                  cout << "MIN " << min << ", MAX " << max << " delta " << delta << endl;
                  color = 1;
               }

               if (ranges == NULL)
               {
                  ranges2 = new std::vector<double>();
                  *ranges2 = {0.0, Sampling * NSamples, min - delta, max + delta};
               }
               else
               {
                  ranges2 = ranges;
               }

               if (PMT_SN[pm] == "")
               {
                  cout << "Entro" << pm << " vs " << NChannels << endl;
                  continue;
               }
               else
               {
                  std::cout << "Drawing PMT" << pm << " - evt " << ev << " - timestamp " << Form("%.0f %s", event.GetTimeStamp(), TDatime(event.GetTimeStamp()).AsString()) << " " << event.GetTimeStamp_nsec() << "ns" << std::endl;
                  TH1D *h = new TH1D();
                  // DrawEvent( waveana::Waveform_t(event.GetTime(),event.Vamp[pm],PedRange,TriggerRange,-(ADCDynamicRange/2.0)*50.0, range1, range2),pm,ev,gPad,ranges,BoxOption);  canvas+=1;}
                  //  lets_pause();
                  DrawEvent(EVT.getChannel(pm), pm, ev, gPad, ranges2, BoxOption, Rebin, color);
                  canvas += 1;
               } // peta aquí
               /*          gPad->SetTickx(2);
                gPad->SetTicky(2);
                th[pm]=event.GetWaveform(pm);th[pm]->SetLineColor(1);
                th[pm]->Draw("HIST"); if(getWaveform(ev,pm)->getMaxAmplitude()>400) th[pm]->SetLineColor(2);
                //gPad->SetLogy();
                th[pm]->GetXaxis()->SetLabelSize(0.08);
                th[pm]->GetYaxis()->SetLabelSize(0.08);
      */
               c->cd(pm + 1)->Modified();
               c->cd(pm + 1)->Update();
            }
            int print = lets_pause(ev);
            if (print == 1)
               c->Print(Form("Run%iEvent%i.png", Number, ev));
	    if (print==2) break;
         }
      }
      int GEOGetPMTBin(string SN, TH2 *h)
      {
         if (SN == "FA0104")
            return h->GetBin(9, 7);
         if (SN == "FA0105")
            return h->GetBin(9, 1);
         if (SN == "FA0106")
            return h->GetBin(7, 11);
         if (SN == "FA0107")
            return h->GetBin(11, 9);
         if (SN == "FA0110")
            return h->GetBin(5, 9);
         if (SN == "FA0111")
            return h->GetBin(9, 11);
         if (SN == "FA0112")
            return h->GetBin(13, 5);
         if (SN == "FA0113")
            return h->GetBin(13, 15);
         if (SN == "FA0114")
            return h->GetBin(1, 9);
         if (SN == "FA0115")
            return h->GetBin(5, 13);
         if (SN == "FA0116")
            return h->GetBin(11, 3);
         if (SN == "FA0119")
            return h->GetBin(3, 11);
         if (SN == "FA0121")
            return h->GetBin(1, 13);
         if (SN == "FA0122")
            return h->GetBin(7, 1);
         if (SN == "FA0124")
            return h->GetBin(5, 7);
         if (SN == "FA0129")
            return h->GetBin(9, 5);
         if (SN == "FA0130")
            return h->GetBin(7, 7);
         if (SN == "FA0132")
            return h->GetBin(7, 9);
         if (SN == "FA0133")
            return h->GetBin(13, 11);
         if (SN == "FA0134")
            return h->GetBin(15, 7);
         if (SN == "FA0135")
            return h->GetBin(1, 7);
         if (SN == "FA0136")
            return h->GetBin(7, 15);
         if (SN == "FA0137")
            return h->GetBin(15, 9);
         if (SN == "FA0139")
            return h->GetBin(9, 9);
         if (SN == "FA0146")
            return h->GetBin(3, 15);
         if (SN == "FA0147")
            return h->GetBin(9, 15);
         if (SN == "FA0148")
            return h->GetBin(3, 9);
         if (SN == "FA0149")
            return h->GetBin(5, 3);
         if (SN == "FA0150")
            return h->GetBin(11, 13);
         if (SN == "FA0151")
            return h->GetBin(13, 9);
         if (SN == "FA0153")
            return h->GetBin(11, 7);
         if (SN == "FA0155")
            return h->GetBin(7, 5);
         if (SN == "FA0156")
            return h->GetBin(15, 13);
         if (SN == "FA0157")
            return h->GetBin(15, 3);
         if (SN == "FC0004")
            return h->GetBin(1, 3);
         if (SN == "FC0005")
            return h->GetBin(3, 5);
         return 0;
      }
      void PMTposition(string SN, TCanvas *c)
      {
         if (SN == "FA0146")
            c->cd(2);
         if (SN == "FA0136")
            c->cd(4);
         if (SN == "FA0147")
            c->cd(5);
         if (SN == "FA0113")
            c->cd(7);
         if (SN == "FA0121")
            c->cd(9);
         if (SN == "FA0115")
            c->cd(11);
         if (SN == "FA0150")
            c->cd(14);
         if (SN == "FA0156")
            c->cd(16);
         if (SN == "FA0119")
            c->cd(18);
         if (SN == "FA0106")
            c->cd(20);
         if (SN == "FA0111")
            c->cd(21);
         if (SN == "FA0133")
            c->cd(23);
         if (SN == "FA0114")
            c->cd(25);
         if (SN == "FA0148")
            c->cd(26);
         if (SN == "FA0110")
            c->cd(27);
         if (SN == "FA0132")
            c->cd(28);
         if (SN == "FA0139")
            c->cd(29);
         if (SN == "FA0107")
            c->cd(30);
         if (SN == "FA0151")
            c->cd(31);
         if (SN == "FA0137")
            c->cd(32);
         if (SN == "FA0135")
            c->cd(33);
         if (SN == "FA0124")
            c->cd(35);
         if (SN == "FA0130")
            c->cd(36);
         if (SN == "FA0104")
            c->cd(37);
         if (SN == "FA0153")
            c->cd(38);
         if (SN == "FA0134")
            c->cd(40);
         if (SN == "FC0005")
            c->cd(42);
         if (SN == "FA0155")
            c->cd(44);
         if (SN == "FA0129")
            c->cd(45);
         if (SN == "FA0112")
            c->cd(47);
         if (SN == "FC0004")
            c->cd(49);
         if (SN == "FA0149")
            c->cd(51);
         if (SN == "FA0116")
            c->cd(54);
         if (SN == "FA0157")
            c->cd(56);
         if (SN == "FA0122")
            c->cd(60);
         if (SN == "FA0105")
            c->cd(61);
      }
      void LoopWaveformsGEO(unsigned int myev = 0, string BoxOption = "", std::vector<double> *ranges = NULL, int Rebin = 1, int CutOption = -1)
      {
         /*CutOption:
         -1 (default) : Loop over all events, skipping cut events.
         n (channel #): Loop over all events, skipping cut events, and events with waveform of channel n cut.
         -2 : Loop over all events (not skipping any event), you will see a text pannel on cut waveforms though. */
         std::cout << "Drawing waveform" << std::endl;
         TH1F *th[64];
         TCanvas *c = new TCanvas("c");
         c->Divide(8, 8);

         if (myev > NEventsALL)
         {
            std::cout << "ERROR: You are trying to draw event " << myev << ", but there is only " << NEvents << " . Exiting ... " << std::endl;
            throw std::exception();
         }
         for (unsigned int ev = myev; ev < NEventsALL; ev++)
         {
            cout << "Looking at " << ev << endl;
            //        if(EventList[17][ev].getPeakTime()<3.2e-6 || EventList[17][ev].getPeakTime()>3.4e-6) continue;
            GetEntry(ev);
            GetEntry(ev);
            int canvas = 0;

            double max = 0;

            Event_t *EVT;
            if (ev < NEvents)
            {
               if (CutOption != -2)
               {
                  if (EventList[ev].getCut())
                     continue;
                  if (CutOption != -1)
                     if (getWaveform(ev, CutOption)->getCut())
                        continue;
               }
               for (auto pm : adcchannels)
               {
                  if (max < getWaveform(ev, pm)->getMaxAmplitude())
                     max = getWaveform(ev, pm)->getMaxAmplitude();
               }
               EVT = getEvent(ev);
            }
            else
            {
               //           Event_t EVT;
               std::cout << "creating evt:" << std::endl;
               EVT = new Event_t();
               for (auto pm : adcchannels)
               {
                  std::cout << "adding ch " << pm << std::endl;
                  EVT->AddWaveform(waveana::Waveform_t(event.GetTime(), event.GetAmp(pm), ParSet), pm);
               }
               std::cout << "creating , applying cuts:" << std::endl;
               MyCuts.ApplyMultiChannelCuts(*EVT);
               std::cout << "evt created, applied cuts" << std::endl;
               if (CutOption != -2)
               {
                  if (EVT->getCut())
                     continue;
                  if (CutOption != -1)
                     if (EVT->getChannel(CutOption)->getCut())
                        continue;
               }
               for (auto pm : adcchannels)
               {
                  if (max < EVT->getChannel(pm)->getMaxAmplitude())
                     max = EVT->getChannel(pm)->getMaxAmplitude();
               }
            }

            for (auto pm : adcchannels)
            {
               PMTposition(PMT_SN[pm], c);
               if (ranges == NULL)
               {
                  ranges = new std::vector<double>();
                  *ranges = {0.0, Sampling * NSamples, 4000 - 1.15 * max, 4050};
               }
               if (PMT_SN[pm] == "")
               {
                  cout << "Entro" << pm << " vs " << NChannels << endl;
                  continue;
               }
               else
               {
                  std::cout << "Drawing PMT" << pm << " - evt " << ev << " - timestamp " << event.GetTimeStamp() << std::endl;
                  TH1D *h = new TH1D();
                  // DrawEvent( waveana::Waveform_t(event.GetTime(),event.Vamp[pm],PedRange,TriggerRange,-(ADCDynamicRange/2.0)*50.0, range1, range2),pm,ev,gPad,ranges,BoxOption);  canvas+=1;}
                  DrawEvent(EVT->getChannel(pm), pm, ev, gPad, NULL, BoxOption, Rebin);
                  canvas += 1;
               }

               /*          gPad->SetTickx(2);
                         gPad->SetTicky(2);
                         th[pm]=event.GetWaveform(pm);th[pm]->SetLineColor(1);
                         th[pm]->Draw("HIST"); if(getWaveform(ev,pm)->getMaxAmplitude()>400) th[pm]->SetLineColor(2);
                         //gPad->SetLogy();
                         th[pm]->GetXaxis()->SetLabelSize(0.08);
                         th[pm]->GetYaxis()->SetLabelSize(0.08);
               */
               c->cd(pm + 1)->Modified();
               c->cd(pm + 1)->Update();
            }
            int print = lets_pause(ev);
            if (print == 1)
               c->Print(Form("Run%iEvent%i.png", Number, ev));
         }
      }
      void LoopWaveformsAll() // All waveforms in the same canvas
      {

         std::vector<TH1D *> th;
         th.resize(NChannels);
         TCanvas *c = new TCanvas("c");
         c->Divide(2, 1);

         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            std::cout << "Evento " << ev << std::endl;
            GetEntry(ev);
            for (auto pm : adcchannels)
            {
               c->cd(pm + 1);
               gPad->SetTickx(2);
               gPad->SetTicky(2);
               th[pm] = event.GetWaveform(pm);
               th[pm]->SetLineColor(1);
               th[pm]->Draw("HIST");
               if (getWaveform(ev, pm)->getMaxAmplitude() > 400)
                  th[pm]->SetLineColor(2);
               // gPad->SetLogy();
               th[pm]->GetXaxis()->SetLabelSize(0.08);
               th[pm]->GetYaxis()->SetLabelSize(0.08);
               c->cd(pm + 1)->Modified();
               c->cd(pm + 1)->Update();
            }
            lets_pause();
            for (auto pm : adcchannels)
            {
               th[pm]->Clear();
               th[pm]->Delete();
            }
         }
      }

      void DrawEvent(waveana::Waveform_t *myevt, int pm, int ev, TVirtualPad *c, std::vector<double> *ranges = NULL, string BoxOption = "", int Rebin = 1, int color = 1)
      {
         TGaxis::SetMaxDigits(3);
         c->cd();
         gStyle->SetOptStat(0);
         string tpb = "";
         if (HasTPB(PMT_SN[pm]))
            tpb = "- TPB";
         TH1D *h = TH1getWaveform(ev, pm);
         if (Rebin != 1)
         {
            h->Rebin(Rebin);
            h->Scale(1. / Rebin);
         }
         /*      if(BoxOption.find("o") !=std::string::npos)//correcting overshooting!
               {
                  ana::AnalyzerS2_t ana(this);
                  h=ana.CorrectOvershooting(ev,pm);
               }
         */
         h->SetTitle(Form("Run %i - ADC Channel %i - Event %i %s", Number, pm, ev, tpb.c_str()));
         h->SetLineColor(color);
         h->Draw("HIST");
         cout << h->GetEntries() << endl;
         h->GetXaxis()->SetLabelSize(0.07);
         h->GetYaxis()->SetLabelSize(0.07);
         h->GetXaxis()->SetTitleSize(0.07);
         h->GetXaxis()->SetTitleSize(0.07);
         if (adcchannels.size() < 10)
         {
            h->GetXaxis()->SetLabelSize(0.05);
            h->GetYaxis()->SetLabelSize(0.05);
            h->GetXaxis()->SetTitleSize(0.05);
            h->GetYaxis()->SetTitleSize(0.05);
         }
         TF1 *fped = new TF1("pedestal", Form("%.5f", myevt->getPedestalMean()), -10, 10);
         fped->SetLineColor(2);

         TPaveText *pt = new TPaveText(.35, .15, .90, .15 + BoxOption.length() * 0.055, "NDC");

         if (BoxOption.find("r") != std::string::npos)
         {
            if (BoxOption.find("q") != std::string::npos)
               pt->AddText(Form("Q_{Range} = %.2fpC = %.0fPEs", 1.0e12 * myevt->getChargeRange(), myevt->getChargeRange() / 1.602e-19 / PMT_Gains[pm]));
            if (BoxOption.find("a") != std::string::npos)
               pt->AddText(Form("Amp_{Range} = %.1fADC = %.1fPEs", myevt->getMaxAmplitudeRange(), myevt->getMaxAmplitudeRange() / PMT_SPEAmp[pm]));
         }
         else
         {
            if (BoxOption.find("q") != std::string::npos)
               pt->AddText(Form("Q_{Tot} = %.2fpC = %.0fPEs", 1.0e12 * myevt->getCharge(), myevt->getCharge() / 1.602e-19 / PMT_Gains[pm]));
            if (BoxOption.find("a") != std::string::npos)
               pt->AddText(Form("Amp = %.1fADC", myevt->getMaxAmplitude()));
         }
         if (BoxOption.find("p") != std::string::npos)
            pt->AddText(Form("Ped = %.1f #pm %.1fADC", myevt->getPedestalMean(), myevt->getPedestalSTD()));
         //   pt->AddText(Form("PedSTD %.1fADC",myevt.getPedestalSTD()));
         //   pt->AddText(Form("Width <%ins",16));
         //   pt->AddText(Form("PedSTD %.1fADC, MaxAmp = %.1fADC",myevt.getPedestalSTD(),myevt.getMaxAmplitude()));
         //   pt->AddText(Form("Q_{Ped} = %.2fpC - Q_{Trigger} = %.2fpC",1.e12*myevt.getPedestalCharge(),1.e12*myevt.getTriggerCharge()));
         //   pt->AddText(Form("SigToNoise_{Q} = %.2f, SigToNoise_{Amp} = %.1f",myevt.getSignalToNoise(),myevt.getSignalToNoiseAmplitude()));
         if (BoxOption.length() != 0)
            pt->Draw();
         fped->Draw("SAME");
         // lets_pause();

         TPaveText *cut = new TPaveText(.25, .25, .75, .75, "NDC");
         cut->AddText(Form("CUT"));
         cut->SetFillStyle(0);
         cut->SetBorderSize(0);
         cut->SetLineWidth(0);
         ((TText *)cut->GetListOfLines()->Last())->SetTextColor(kRed);
         if (myevt->getCut())
            cut->Draw();
         Float_t ymax = h->GetMaximum();
         Float_t ymin = h->GetMinimum();
         TLine *line1 = new TLine(event.GetTime()->at(0) + PedRange * myevt->getSampling(), myevt->getPedestalMean() - 5, event.GetTime()->at(0) + PedRange * myevt->getSampling(), myevt->getPedestalMean() + 5);
         TLine *line2 = new TLine(range1, ymin, range1, ymax);
         TLine *line3 = new TLine(range2, ymin, range2, ymax);
         line1->SetLineColor(kRed);
         line2->SetLineWidth(2);
         line2->SetLineColor(kGreen);
         line3->SetLineColor(kGreen);
         line1->Draw();
         line2->Draw();
         line3->Draw();
         myevt->Print(); // cout << ranges->at(0) << " " << ranges->at(1)<< " " << ranges->at(2)<< " " << ranges->at(3) << endl;
         if (ranges)
         {
            cout << "ESTOY!!!" << endl;
            h->GetXaxis()->SetRangeUser(ranges->at(0), ranges->at(1));
            h->GetYaxis()->SetRangeUser(ranges->at(2), ranges->at(3));
         }
         h->GetYaxis()->UnZoom();
      }

      void Plot36(string mode, string dumpfile = "", int logoption = 0, bool pause = false, string fileoption = "RECREATE", double range1 = 0, double range2 = 0, std::vector<int> *SelectedChannels = NULL)
      {
         std::cout << "Plotting " << mode << endl;
         Run_t *run = this;
         std::vector<int> ThisChannels = run->adcchannels;
         if (SelectedChannels)
            ThisChannels = *SelectedChannels;

         HistogramCollection_t *HC;
         HC = run->GetHistCollectionByMode(mode, range1, range2, &ThisChannels);
         std::vector<TH1 *> th;
         th = HC->GetByOpChannel64();

         TCanvas *c = new TCanvas("c");
         float txtsize = 0.06;
         c->DivideSquare(run->adcchannels.size());
         int canvas = 1;
         TGaxis::SetMaxDigits(3);
         for (auto pm : ThisChannels)
         {
            c->cd(canvas);
            th[pm]->Draw(HC->GetDrawOption().c_str());
            // if(mode=="ScintProf"){th[pm]->GetXaxis()->SetRangeUser(0.e-6,8.5e-6);gStyle->SetOptFit(1);gStyle->SetOptStat(0);th[pm]->GetListOfFunctions()->FindObject("scint_time")->Draw("SAME");}
            // if(mode=="AmpVsCharge") {th[pm]->Draw("COL");th[pm]->GetYaxis()->SetRangeUser(0,80);}
            //   gPad->SetTickx(2);
            //   gPad->SetTicky(2);
            if (logoption == 2 || logoption == 3)
               gPad->SetLogy();
            if (logoption == 1 || logoption == 3)
               gPad->SetLogx();
            gPad->SetGridx();
            gPad->SetGridy();
            th[pm]->GetXaxis()->SetLabelSize(txtsize);
            th[pm]->GetXaxis()->SetNdivisions(5);
            th[pm]->GetYaxis()->SetLabelSize(txtsize);
            th[pm]->GetXaxis()->SetTitleSize(txtsize);
            th[pm]->GetYaxis()->SetTitleSize(txtsize);
            c->cd(canvas)->Modified();
            c->cd(canvas)->Update();
            canvas++;
         }
         if (pause)
            lets_pause();
         if (dumpfile != "")
         {
            TFile *dump = new TFile(dumpfile.c_str(), fileoption.c_str());
            dump->cd();
            for (auto pm : ThisChannels)
            {
               th[pm]->Write();
            }
            dump->Close();
            std::cout << "Dumped " << mode << " to: " << dumpfile << std::endl;
         }
         for (auto pm : ThisChannels)
            delete th[pm];
      }

      void DumpVariableToNtuple(string filename, bool debug = false)
      {
         TFile *ofile = new TFile(filename.c_str(), "RECREATE");
         ofile->cd();
         TNtuple *nt = new TNtuple("ntuple", "ntuple", Form("evt:ch:time:Amp:AmpPE:Q1:Q2:Q3:QPeak:TStartQpeak:TEndQPeak:TDiff:QTotal"));
         TNtuple *nt_Q = new TNtuple("charge", "charge", "QFixRange1:QFixRange2:QFixRange3:QFixRange4:QPeakRange1:QPeakRange2:QPeakRange3:QPeakRange4");

         for (unsigned int ev = 0; ev < NEvents; ev++)
         {
            GetEntry(ev);
            if (EventList[ev].getCut())
               continue;
            for (auto pm : adcchannels)
            {
               if (getWaveform(ev, pm)->getCut())
                  continue;
               auto vQFixed = getWaveform(ev, pm)->getQFixed_ranges();
               auto vQPeak = getWaveform(ev, pm)->getQPeak_ranges();

               nt->Fill(ev, pm, event.GetTimeStamp(), getWaveform(ev, pm)->getMaxAmplitudeRange(), getWaveform(ev, pm)->getMaxAmplitudeRange() / PMT_SPEAmp[pm], getWaveform(ev, pm)->getQ1MaxPeakRange() * 1e12, getWaveform(ev, pm)->getQ2MaxPeakRange() * 1e12, getWaveform(ev, pm)->getQ3MaxPeakRange() * 1e12, getWaveform(ev, pm)->getChargeMaxPeakRange() * 1e12, getWaveform(ev, pm)->getTStartMaxPeakRange(), getWaveform(ev, pm)->getTEndMaxPeakRange(), getWaveform(ev, pm)->getTEndMaxPeakRange() - getWaveform(ev, pm)->getTStartMaxPeakRange(), getWaveform(ev, pm)->getCharge() * 1e12);
               nt_Q->Fill(
                   vQFixed[0] * 1e12,
                   vQFixed[1] * 1e12,
                   vQFixed[2] * 1e12,
                   vQFixed[3] * 1e12,
                   vQPeak[0] * 1e12,
                   vQPeak[1] * 1e12,
                   vQPeak[2] * 1e12,
                   vQPeak[3] * 1e12);
               if (debug)
               {
                  cout << ev << " " << pm << " " << event.GetTimeStamp() << endl;
                  cout << getWaveform(ev, pm)->getMaxAmplitudeRange() << endl;
                  cout << getWaveform(ev, pm)->getMaxAmplitudeRange() / PMT_SPEAmp[pm] << endl;
                  cout << getWaveform(ev, pm)->getQ1MaxPeakRange() * 1e12 << endl;
                  cout << getWaveform(ev, pm)->getQ2MaxPeakRange() * 1e12 << endl;
                  cout << getWaveform(ev, pm)->getQ3MaxPeakRange() * 1e12 << endl;
                  cout << getWaveform(ev, pm)->getChargeMaxPeakRange() * 1e12 << endl;
                  cout << getWaveform(ev, pm)->getTStartMaxPeakRange() << endl;
                  cout << getWaveform(ev, pm)->getTEndMaxPeakRange() << endl;
                  lets_pause();
               };
            }
         }
         nt->Write("ntuple");
         nt_Q->Write("charge");
         ofile->Close();
         cout << "ntuple dumped to " << filename << endl;
      }
      
      void autofit(TH1F* h, bool print, string output)
      { // Función para realizar ajustes automaticos
         int npeaks = 0; cout << "N of peaks?: "; cin >> npeaks;
         double_t par[100];
         double min_fit, max_fit;
         cout << "Min value for fit: "; cin >> min_fit;
         cout << "Max value for fit: "; cin >> max_fit;
         // Use TSpectrum to find the peak candidates
         h->GetXaxis()->SetRangeUser(min_fit, max_fit);
         TSpectrum *s = new TSpectrum(2*npeaks,1);
         Int_t nfound = s->Search(h,2,"goff",0.01);
         printf("Found %d candidate peaks to fit\n",nfound);
         h->GetXaxis()->SetRange(0,0);
         Double_t *xpeaks;
         xpeaks = s->GetPositionX();
      
         for (int p=0;p<npeaks;p++)
         {
            Double_t xp = xpeaks[p];
            Int_t bin = h->GetXaxis()->FindBin(xp);
            Double_t yp = h->GetBinContent(bin);
            //if (yp-TMath::Sqrt(yp) < fline->Eval(xp)) continue;
            par[3*p] = yp; // "height"
            par[3*p+1] = xp; // "mean"
            par[3*p+2] = 0.05; // "sigma"
            // cout << Form("sigma_%i: ", p); cin >> par[3*npeaks+2];
            //~ p++;
         }
	
         //int_t max_bin[npeaks], j=0;
         string init = Form("gaus(0)");
         string added;
         for (int i=1; i<npeaks; i++)
         {
            added = Form("+gaus(%i)",i*3);
            init = init + added;   
         }
	      //int Nbins = h->GetNbinsX();
	
         TF1 *fb = new TF1("fb",init.c_str(),min_fit,max_fit);
         fb->SetParameters(par);

         h->Fit("fb","MER+","HIST SAME",min_fit, max_fit); gPad->Update();
         lets_pause();
      
         cout << "Printing mean values: " << endl;
         for (int p=0; p<npeaks;p++)
         {
            cout << Form("mu_%i: ", p) << fb->GetParameter(3*p+1) << " +- " << fb->GetParError(3*p+1) << endl;
         }
         cout << "Printing sigma values: " << endl;
         for (int p=0; p<npeaks;p++)
         {
            cout << Form("sigma_%i: ", p) << fb->GetParameter(3*p+2) << " +- " << fb->GetParError(3*p+2) << endl;
         }
         
         lets_pause();

         // ------------ Volcamos a txt ---------------
         if (print == true)
         {
            string filename = output+"_GAUSS_FIT";
            ifstream ifile;
            ifile.open("fit_data/"+filename+".txt");
            string str = "fit_data/"+filename+".txt";
            char c[str.size() + 1];
            strcpy(c, str.c_str());
            if(ifile){remove(c);}
            for (int p=0; p<npeaks;p++)
            {
               double mu = fb->GetParameter(3*p+1), Dmu = fb->GetParError(3*p+1);
               double sigma = fb->GetParameter(3*p+2), Dsigma = fb->GetParError(3*p+2);
               double amp = fb->GetParameter(p), Damp = fb->GetParError(p);
         
               FILE *f = fopen(Form("fit_data/%s.txt",filename.c_str()), "a");
               if (f == NULL)
               {
                  printf("Error opening file!\n");
                  exit(1);
               }
               fprintf(f, "%i\t%f\t%f\t%f\t%f\n", p, mu, Dmu, sigma, Dsigma);
               fclose(f);
            }
            cout << "\nFile with name '"; cout << filename; cout << ".txt' has been generated.\n\n";

            if (npeaks>1)
            {
               string filename = output+"_GAIN";
               ifstream ifile;
               ifile.open("fit_data/"+filename+".txt");
               string str = "fit_data/"+filename+".txt";
               char c[str.size() + 1];
               strcpy(c, str.c_str());
               if(ifile){remove(c);}
               for (int p=1; p<npeaks;p++)
               {  
                  double mu0 = fb->GetParameter(3*(p-1)+1), Dmu0 = fb->GetParError(3*(p-1)+1);
                  double sigma0 = fb->GetParameter(3*(p-1)+2), Dsigma0 = fb->GetParError(3*(p-1)+2);
                  double mu = fb->GetParameter(3*p+1), Dmu = fb->GetParError(3*p+1);
                  double sigma = fb->GetParameter(3*p+2), Dsigma = fb->GetParError(3*p+2);
                  
                  double gain = (mu-mu0)*1e-12/1.602e-19;
                  double Dgain = pow((Square(Dmu)+Square(Dmu0)),0.5)*1e-12/1.602e-19;
                  double sn0 = (mu-mu0)/sigma0;
                  double Dsn0 = sn0*pow(((Square(Dmu)+Square(Dmu0))/Square(mu-mu0))+Square(Dsigma0/sigma0),0.5);
                  double sn1 = (mu-mu0)/sigma;
                  double Dsn1 = sn1*pow(((Square(Dmu)+Square(Dmu0))/Square(mu-mu0))+Square(Dsigma/sigma),0.5);
                  double snc = (mu-mu0)/pow((Square(sigma)+Square(sigma0)),0.5);
                  double Dsnc = snc*pow((Square(Dgain/gain)+Square(sigma0*Dsigma0/(Square(sigma)+Square(sigma0)))+Square(sigma*Dsigma/(Square(sigma)+Square(sigma0)))),0.5);

                  FILE *f = fopen(Form("fit_data/%s.txt",filename.c_str()), "a");
                  if (f == NULL)
                  {
                     printf("Error opening file!\n");
                     exit(1);
                  }
                  fprintf(f, "%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", p, gain, Dgain, sn0, Dsn0, sn1, Dsn1, snc, Dsnc);
                  fclose(f);
               }
               cout << "\nFile with name '"; cout << filename; cout << ".txt' has been generated.\n\n";
            }
            
         }

      }
      void HistEventCounter(TH1F* h)
      { // Función contar eventos en un histograma
         cout << "\nEnter rangevalues for event counting\n";
         int peaks = 1; cout << "Number of ranges?: "; cin >> peaks;

         for (int x = 1; x <= peaks; x++)
         {
            double min_value = 0; cout << "Min value?: "; cin >> min_value;
            Int_t min_bin = h->GetXaxis()->FindBin(min_value);
            double max_value = 0; cout << "Max value?: "; cin >> max_value;
            Int_t max_bin = h->GetXaxis()->FindBin(max_value);

            int counts = 0;
            for (int p=min_bin;p<max_bin;p++)
            {
               Double_t yp = h->GetBinContent(p);
               counts = counts + yp;
            }
            cout << "Total number of events between "; cout << min_value; cout << " and "; cout << max_value; cout << " = "; cout << counts; cout << "\n";
         }
             
      }
   };
}
#endif
