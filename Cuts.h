namespace ana
{
   class Cuts_t
   {
   public:
      bool SaturationCut = false;
      bool CutPedSTD = false;
      std::vector<double> PedestalSTDCUT;
      bool CutBaselineStability = false;
      double BaselineStability_sigma;
      std::vector<double> BaselineStability_PedPeak;
      std::vector<double> BaselineStability_PedSigma;
      bool CutMaxAmplitude = false;
      double MaxAmplitude_minAmp, MaxAmplitude_maxAmp;
      bool CutMaxAmplitudeRange = false;
      double MaxAmplitudeRange_minAmp, MaxAmplitudeRange_maxAmp;
      bool CutRemoveOutOfRangePeaksAbove = false;
      double RemoveOutOfRangePeaksAbove_maxAmp;
      bool CutRemovePeaksOutOfMaxPeakAbove = false;
      double RemovePeaksOutOfMaxPeakAbove_maxAmp;
      bool CutPeakTimeRange = false;
      double PeakTimeRange_minTime;
      double PeakTimeRange_maxTime;
      bool CutS2Amplitude = false;
      double S2Amplitude_minAmp;
      bool CutCharge = false;
      double Charge_minCharge;
      double Charge_maxCharge;
      bool CutChargeRange = false;
      double ChargeRange_minCharge;
      double ChargeRange_maxCharge;
      bool CutMinimumMaximumSample = false;
      double MinimumMaximumSample_Sample;

      bool CutTriggerAmplitudeRange = false;
      double TriggerAmplitudeRange_minAmp;
      double TriggerAmplitudeRange_maxAmp;
      std::vector<int> TriggerAmplitudeRange_TriggerChannels;
      bool CutTriggerFirstSampleBelowADCTriggerThresholdRange = false;
      double TriggerFirstSampleBelowADCTriggerThresholdRange_ADCThreshold;
      std::vector<int> TriggerFirstSampleBelowADCTriggerThresholdRange_TriggerChannels;

      bool CutNsamplesBelowThresholdRebin = false;
      int NsamplesBelowThresholdRebin_nsamples, NsamplesBelowThresholdRebin_ADCThreshold, NsamplesBelowThresholdRebin_Rebin;

      bool CutTriggerMinimumThresholdRange = false;
      double TriggerMinimumThresholdRange_maxSample;
      std::vector<int> TriggerMinimumThresholdRange_TriggerChannels;
      bool CutThresholdRange = false;
      double ThresholdRange_maxSample;
      bool CutThreshold = false;
      double Threshold_maxSample;
      bool CutPedestalStatus = false;
      bool CutOvershootingRange = false;
      int OvershootingRange_Sample;
      bool CutTriggerWaveformCuts = false;
      std::vector<int> TriggerWaveformCuts_TriggerChannels;
      bool CutVariableVector = false;
      string VariableVector_var;
      std::map<int, std::pair<double, double>> VariableVector_map;
      bool CutVariable = false;
      std::vector<string> Variable_var;
      std::vector<double> Variable_min, Variable_max;

      void SetCutPedestalStatus()
      {
         CutPedestalStatus = true;
      }
      void SetCutTriggerMinimumThresholdRange(double maxSample, std::vector<int> TriggerChannels)
      {
         CutTriggerMinimumThresholdRange = true;
         TriggerMinimumThresholdRange_maxSample = maxSample;
         TriggerMinimumThresholdRange_TriggerChannels = TriggerChannels;
      }
      void SetCutThresholdRange(double maxSample)
      {
         CutThresholdRange = true;
         ThresholdRange_maxSample = maxSample;
      }
      void SetCutThreshold(double maxSample)
      {
         CutThreshold = true;
         Threshold_maxSample = maxSample;
      }

      void SetCutPedSTD(std::vector<double> p)
      {
         CutPedSTD = true;
         PedestalSTDCUT = p;
      }
      void SetCutNsamplesBelowThresholdRebin(int nsamples, int ADCThreshold, int Rebin)
      {
         CutNsamplesBelowThresholdRebin = true;
         NsamplesBelowThresholdRebin_nsamples = nsamples;
         NsamplesBelowThresholdRebin_ADCThreshold = ADCThreshold;
         NsamplesBelowThresholdRebin_Rebin = Rebin;
      }
      void SetCutMinimumMaximumSample(double Sample)
      {
         CutMinimumMaximumSample = true;
         MinimumMaximumSample_Sample = Sample;
      }
      void SetCutBaselineStability(double sigma, std::vector<double> PedPeak, std::vector<double> PedSigma)
      {
         CutBaselineStability = true;
         BaselineStability_sigma = sigma;
         BaselineStability_PedPeak = PedPeak;
         BaselineStability_PedSigma = PedSigma;
      }
      void SetSaturationCut()
      {
         SaturationCut = true;
      }
      void SetCutMaxAmplitude(double minAmp, double maxAmp)
      {
         CutMaxAmplitude = true;
         MaxAmplitude_minAmp = minAmp;
         MaxAmplitude_maxAmp = maxAmp;
      }
      void SetCutMaxAmplitudeRange(double minAmp, double maxAmp)
      {
         CutMaxAmplitudeRange = true;
         MaxAmplitudeRange_minAmp = minAmp;
         MaxAmplitudeRange_maxAmp = maxAmp;
      }
      void SetCutRemoveOutRangePeaksAbove(double maxAmp)
      {
         CutRemoveOutOfRangePeaksAbove = true;
         RemoveOutOfRangePeaksAbove_maxAmp = maxAmp;
      }
      void SetCutRemovePeaksOutOfMaxPeakAbove(double maxAmp) // cut events with secondary signals above maxAmp
      {
         CutRemovePeaksOutOfMaxPeakAbove = true;
         RemovePeaksOutOfMaxPeakAbove_maxAmp = maxAmp;
      }
      void SetCutPeakTimeRange(double minTime, double maxTime) // cut events with peak time out of the range
      {
         CutPeakTimeRange = true;
         PeakTimeRange_minTime = minTime;
         PeakTimeRange_maxTime = maxTime;
      }
      void SetCutS2Amplitude(double minAmp)
      {
         CutS2Amplitude = true;
         S2Amplitude_minAmp = minAmp;
      }
      void SetCutTriggerAmplitudeRange(double minAmp, double maxAmp, std::vector<int> TriggerChannels)
      {
         CutTriggerAmplitudeRange = true;
         TriggerAmplitudeRange_minAmp = minAmp;
         TriggerAmplitudeRange_maxAmp = maxAmp;
         TriggerAmplitudeRange_TriggerChannels = TriggerChannels;
      }
      void SetCutTriggerFirstSampleBelowADCTriggerThresholdRange(double ADCThreshold, std::vector<int> TriggerChannels)
      {
         CutTriggerFirstSampleBelowADCTriggerThresholdRange = true;
         TriggerFirstSampleBelowADCTriggerThresholdRange_ADCThreshold = ADCThreshold;
         TriggerFirstSampleBelowADCTriggerThresholdRange_TriggerChannels = TriggerChannels;
      }
      void SetCutCharge(double minCharge, double maxCharge)
      {
         CutCharge = true;
         Charge_minCharge = minCharge;
         Charge_maxCharge = maxCharge;
      }
      void SetCutChargeRange(double minCharge, double maxCharge)
      {
         CutChargeRange = true;
         ChargeRange_minCharge = minCharge;
         ChargeRange_maxCharge = maxCharge;
      }
      void SetCutOvershootingRange(int Sample)
      {
         CutOvershootingRange = true;
         OvershootingRange_Sample = Sample;
      }
      void SetCutTriggerWaveformCuts(std::vector<int> TriggerChannels)
      {
         CutTriggerWaveformCuts = true;
         TriggerWaveformCuts_TriggerChannels = TriggerChannels;
      }
      void SetCutVariableVector(string var, std::map<int, std::pair<double, double>> map)
      {
         CutVariableVector = true;
         VariableVector_var = var;
         VariableVector_map = map;
      }
      void SetCutVariable(string var, double min, double max)
      {
         CutVariable = true;
         Variable_var.push_back(var);
         Variable_min.push_back(min);
         Variable_max.push_back(max);
      }

      void ResetCuts()
      {
         SaturationCut = false;
         CutPedSTD = false;
         CutMaxAmplitude = false;
         CutMaxAmplitudeRange = false;
         CutRemoveOutOfRangePeaksAbove = false;
         CutRemovePeaksOutOfMaxPeakAbove = false;
         CutPeakTimeRange = false;
         CutTriggerAmplitudeRange = false;
         CutChargeRange = false;
         CutBaselineStability = false;
         CutS2Amplitude = false;
         CutCharge = false;
         CutChargeRange = false;
         CutMinimumMaximumSample = false;
         CutTriggerAmplitudeRange = false;
         CutTriggerFirstSampleBelowADCTriggerThresholdRange = false;
         CutNsamplesBelowThresholdRebin = false;
         CutTriggerMinimumThresholdRange = false;
         CutThresholdRange = false;
         CutThreshold = false;
         CutPedestalStatus = false;
         CutOvershootingRange = false;
         CutTriggerWaveformCuts = false;
         CutVariableVector = false;
         CutVariable = false;
      }

      void ApplyCuts(waveana::Waveform_t &myevt, int pm, TH1D *h = NULL)
      {
         myevt.setCut(false);

         if (CutPedSTD)
            if (myevt.getPedestalSTD() > PedestalSTDCUT[pm])
               myevt.setCut(true);
         if (CutBaselineStability)
            if (myevt.getPedestalMean() > BaselineStability_PedPeak[pm] + BaselineStability_sigma * BaselineStability_PedSigma[pm] || myevt.getPedestalMean() < BaselineStability_PedPeak[pm] - BaselineStability_sigma * BaselineStability_PedSigma[pm])
               myevt.setCut(true);
         if (SaturationCut)
            if (myevt.getSaturation())
               myevt.setCut(true);
         if (CutMaxAmplitude)
            if (myevt.getMaxAmplitude() < MaxAmplitude_minAmp || myevt.getMaxAmplitude() > MaxAmplitude_maxAmp)
               myevt.setCut(true);
         if (CutMaxAmplitudeRange)
            if (myevt.getMaxAmplitudeRange() < MaxAmplitudeRange_minAmp || myevt.getMaxAmplitudeRange() > MaxAmplitudeRange_maxAmp)
               myevt.setCut(true);
         if (CutRemoveOutOfRangePeaksAbove)
            if (myevt.getMaxAmplitudeOutOfRange() > RemoveOutOfRangePeaksAbove_maxAmp)
               myevt.setCut(true);
         if (CutRemovePeaksOutOfMaxPeakAbove)
            if (myevt.getMaxAmplitudeOutOfMaxPeak() > RemovePeaksOutOfMaxPeakAbove_maxAmp)
               myevt.setCut(true);
         if (CutPeakTimeRange)
            if (myevt.getPeakTime() < PeakTimeRange_minTime || myevt.getPeakTime() > PeakTimeRange_maxTime)
               myevt.setCut(true);
         if (CutS2Amplitude)
            if (myevt.getS2AverageAmplitude() < S2Amplitude_minAmp)
               myevt.setCut(true);
         if (CutCharge)
            if (myevt.getCharge() < Charge_minCharge || myevt.getCharge() > Charge_maxCharge)
               myevt.setCut(true);
         if (CutChargeRange)
            if (myevt.getChargeRange() < ChargeRange_minCharge || myevt.getChargeRange() > ChargeRange_maxCharge)
               myevt.setCut(true);
         if (CutMinimumMaximumSample)
            if (myevt.getMaximumSample() <= MinimumMaximumSample_Sample)
               myevt.setCut(true);
         if (CutNsamplesBelowThresholdRebin)
            if (!CheckCutNsamplesBelowThresholdRebin(myevt, h))
               myevt.setCut(true);
         if (CutThresholdRange)
            if (myevt.getMinimumSampleRange() > ThresholdRange_maxSample)
            {
               myevt.setCut(true);
            }
         if (CutThreshold)
            if (myevt.getMinimumSample() > Threshold_maxSample)
            {
               myevt.setCut(true);
            }
         if (CutPedestalStatus)
            if (!myevt.getPedestalStatus())
            {
               myevt.setCut(true);
            }
         if (CutOvershootingRange)
            if (myevt.getMaximumSample_Range() - myevt.getPedestalMean() > OvershootingRange_Sample)
            {
               myevt.setCut(true);
            }

         if (CutVariableVector)
         {
            if (VariableVector_var == "PreTriggerSTD")
            {
               if (myevt.getPreTriggerSTD() < VariableVector_map[pm].first || myevt.getPreTriggerSTD() > VariableVector_map[pm].second)
                  myevt.setCut(true);
            }
            else if (VariableVector_var == "MaxAmplitudeRange")
            {
               if (myevt.getMaxAmplitudeRange() < VariableVector_map[pm].first || myevt.getMaxAmplitudeRange() > VariableVector_map[pm].second)
                  myevt.setCut(true);
            }
            else if (VariableVector_var == "ChargeMaxPeakRange")
            {
               if (myevt.getChargeMaxPeakRange() < VariableVector_map[pm].first || myevt.getChargeMaxPeakRange() > VariableVector_map[pm].second)
                  myevt.setCut(true);
            }
            else if (VariableVector_var == "Q3MaxPeakRange")
            {
               if (myevt.getQ3MaxPeakRange() < VariableVector_map[pm].first || myevt.getQ3MaxPeakRange() > VariableVector_map[pm].second)
                  myevt.setCut(true);
            }
            else if (VariableVector_var == "TEndMaxPeakRange")
            {
               if (myevt.getTEndMaxPeakRange() < VariableVector_map[pm].first || myevt.getTEndMaxPeakRange() > VariableVector_map[pm].second)
                  myevt.setCut(true);
            }
            else
            {
               cout << "Error in cut, var not defined" << endl;
               throw std::exception();
            }
         }
         if (CutVariable)
         {
            for (unsigned int i = 0; i < Variable_var.size(); i++)
            {
               if (Variable_var[i] == "PreTriggerSTD")
               {
                  if (myevt.getPreTriggerSTD() < Variable_min[i] || myevt.getPreTriggerSTD() > Variable_max[i])
                     myevt.setCut(true);
               }
               else if (Variable_var[i] == "MaxAmplitudeRange")
               {
                  if (myevt.getMaxAmplitudeRange() < Variable_min[i] || myevt.getMaxAmplitudeRange() > Variable_max[i])
                     myevt.setCut(true);
               }
               else if (Variable_var[i] == "ChargeMaxPeakRange")
               {
                  if (myevt.getChargeMaxPeakRange() < Variable_min[i] || myevt.getChargeMaxPeakRange() > Variable_max[i])
                     myevt.setCut(true);
               }
               else if (Variable_var[i] == "Q3MaxPeakRange")
               {
                  if (myevt.getQ3MaxPeakRange() < Variable_min[i] || myevt.getQ3MaxPeakRange() > Variable_max[i])
                     myevt.setCut(true);
               }
               else if (Variable_var[i] == "TEndMaxPeakRange")
               {
                  if (myevt.getTEndMaxPeakRange() < Variable_min[i] || myevt.getTEndMaxPeakRange() > Variable_max[i])
                     myevt.setCut(true);
               }
               else
               {
                  cout << "Error in cut, var not defined" << endl;
                  throw std::exception();
               }
            }
         }
      }
      bool CheckCutNsamplesBelowThresholdRebin(waveana::Waveform_t &myevt, TH1D *hh)
      {

         TH1D *h = (TH1D *)hh->Clone("htemp");
         h->Rebin(NsamplesBelowThresholdRebin_Rebin);
         h->Scale(1.0 / NsamplesBelowThresholdRebin_Rebin);
         int counter = 0;
         for (int i = 1; i <= h->GetSize(); i++)
         {
            if (h->GetBinContent(i) < NsamplesBelowThresholdRebin_ADCThreshold)
               counter++;
            else
               counter = 0;
            if (counter >= NsamplesBelowThresholdRebin_nsamples)
            {
               delete h;
               return true;
            }
         }
         delete h;
         return false;
      }

      void ApplyMultiChannelCuts(Event_t &EVT)
      {

         bool CutTriggerAmplitudeRange_ThisEVT = false;
         int pm = 0;
         for (auto myevt : EVT.waveform)
         {
            myevt.setCut(false);
            if (CutPedSTD)
               if (myevt.getPedestalSTD() > PedestalSTDCUT[pm])
               {
                  myevt.setCut(true);
               }
            if (SaturationCut)
               if (myevt.getSaturation())
               {
                  myevt.setCut(true);
               }
            if (CutMaxAmplitude)
               if (myevt.getMaxAmplitude() < MaxAmplitude_minAmp || myevt.getMaxAmplitude() > MaxAmplitude_maxAmp)
               {
                  myevt.setCut(true);
               }
            if (CutMaxAmplitudeRange)
               if (myevt.getMaxAmplitudeRange() < MaxAmplitudeRange_minAmp || myevt.getMaxAmplitudeRange() > MaxAmplitudeRange_maxAmp)
               {
                  myevt.setCut(true);
               }
            if (CutRemoveOutOfRangePeaksAbove)
               if (myevt.getMaxAmplitudeOutOfRange() > RemoveOutOfRangePeaksAbove_maxAmp)
               {
                  myevt.setCut(true);
               }
            if (CutRemovePeaksOutOfMaxPeakAbove)
               if (myevt.getMaxAmplitudeOutOfMaxPeak() > RemovePeaksOutOfMaxPeakAbove_maxAmp)
               {
                  myevt.setCut(true);
               }
            if (CutPeakTimeRange)
               if (myevt.getPeakTime() < PeakTimeRange_minTime || myevt.getPeakTime() > PeakTimeRange_maxTime)
               {
                  myevt.setCut(true);
               }
            if (CutS2Amplitude)
               if (myevt.getS2AverageAmplitude() < S2Amplitude_minAmp)
               {
                  myevt.setCut(true);
               }
            if (CutCharge)
               if (myevt.getCharge() < Charge_minCharge || myevt.getCharge() > Charge_maxCharge)
               {
                  myevt.setCut(true);
               }
            if (CutChargeRange)
               if (myevt.getChargeRange() < ChargeRange_minCharge || myevt.getChargeRange() > ChargeRange_maxCharge)
               {
                  myevt.setCut(true);
               }
            if (CutThresholdRange)
               if (myevt.getMinimumSampleRange() > ThresholdRange_maxSample)
               {
                  myevt.setCut(true);
               }
            if (CutThreshold)
               if (myevt.getMinimumSample() > Threshold_maxSample)
               {
                  myevt.setCut(true);
               }
            if (CutPedestalStatus)
               if (!myevt.getPedestalStatus())
               {
                  myevt.setCut(true);
               }
            if (CutOvershootingRange)
               if (myevt.getMaximumSample_Range() - myevt.getPedestalMean() > OvershootingRange_Sample)
               {
                  myevt.setCut(true);
               }
            if (CutVariableVector)
            {
               if (VariableVector_var == "PreTriggerSTD")
               {
                  if (myevt.getPreTriggerSTD() < VariableVector_map[pm].first || myevt.getPreTriggerSTD() > VariableVector_map[pm].second)
                     myevt.setCut(true);
               }
               if (VariableVector_var == "MaxAmplitudeRange")
               {
                  if (myevt.getMaxAmplitudeRange() < VariableVector_map[pm].first || myevt.getMaxAmplitudeRange() > VariableVector_map[pm].second)
                     myevt.setCut(true);
               }
               if (VariableVector_var == "ChargeMaxPeakRange")
               {
                  if (myevt.getChargeMaxPeakRange() < VariableVector_map[pm].first || myevt.getChargeMaxPeakRange() > VariableVector_map[pm].second)
                     myevt.setCut(true);
               }
            }

            pm++;
         }

         if (CutTriggerWaveformCuts)
         {
            for (auto pm : TriggerWaveformCuts_TriggerChannels)
               if (EVT.getChannel(pm)->getCut())
               {
                  EVT.setCut(true);
                  break;
               }
         }

         if (CutTriggerMinimumThresholdRange)
         {
            int c = 0;
            if (EVT.getChannel(pm)->getMinimumSampleRange() < TriggerMinimumThresholdRange_maxSample)
            {
               c++;
            }
            if (c == 0)
               EVT.setCut(true);
         }

         if (CutTriggerAmplitudeRange)
            for (auto pm : TriggerAmplitudeRange_TriggerChannels)
               if (EVT.getChannel(pm)->getMaxAmplitudeRange() < TriggerAmplitudeRange_minAmp || EVT.getChannel(pm)->getMaxAmplitudeRange() > TriggerAmplitudeRange_maxAmp)
               {
                  EVT.setCut(true);
                  break;
               }

         if (CutTriggerFirstSampleBelowADCTriggerThresholdRange)
         {
            if (CheckCutTriggerFirstSampleBelowADCTriggerThresholdRange(EVT))
               EVT.setCut(true);
            else
               EVT.setTriggerTime(EVT.getChannel(TriggerFirstSampleBelowADCTriggerThresholdRange_TriggerChannels[0])->getFirstSampleBelowADCTriggerThresholdRange());
         }
      }

      bool CheckCutTriggerFirstSampleBelowADCTriggerThresholdRange(Event_t &EVT)
      {
         unsigned int t = 0;
         for (auto pm : TriggerFirstSampleBelowADCTriggerThresholdRange_TriggerChannels)
         {
            if (EVT.getChannel(pm)->getMinimumSampleRange() < TriggerFirstSampleBelowADCTriggerThresholdRange_ADCThreshold)
               t++;
         }
         if (t == TriggerFirstSampleBelowADCTriggerThresholdRange_TriggerChannels.size())
            return false; // we keep it, dont cut it
         else
            return true; // lets remove it!
      }
   };
}
