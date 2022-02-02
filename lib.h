namespace ana
{
  class WF_t
  {
  public:
    TFile *f = NULL;
    TTree *t = NULL;
    std::vector<short> *ws = NULL;
    std::vector<double> *wd = NULL;
    int EventNumber;
    ULong64_t TriggerTimeStamp;
    ULong64_t PCTimeStamp;

    string name;
    WF_t() {}
    std::vector<double> Amp;
    std::vector<double> Time;
    int Number;
    int nsamples;
    double sampling = 4.e-9;

    void SetSampling(double s) { sampling = s; }
    double GetSampling() { return sampling; }
    bool IsOscilloscope = false;

    double GetTimeStamp() { return static_cast<double>(TriggerTimeStamp); }
    double GetPCTimeStamp() { return static_cast<double>(PCTimeStamp); }

    void Set(string file, string n, bool IsOsc)
    {

      IsOscilloscope = IsOsc;
      f = TFile::Open(file.c_str(), "READ");
      if (f)
      {
        t = (TTree *)f->Get("IR02");
        if (IsOscilloscope)
          t->SetBranchAddress("ADC", &wd);
        else
          t->SetBranchAddress("ADC", &ws);
        if (IsOscilloscope)
          t->SetBranchAddress("Sampling", &sampling);
        t->SetBranchAddress("EventNumber", &EventNumber);
        t->SetBranchAddress("TriggerTimeStamp", &TriggerTimeStamp);
        t->SetBranchAddress("PCTimeStamp", &PCTimeStamp);
        std::cout << "Branches set for " << file << std::endl;
        name = n;
      }
    }
    void GetEntry(int i)
    {
      t->GetEntry(i);
      Number = i;
      if (IsOscilloscope)
        nsamples = wd->size();
      else
        nsamples = ws->size();
      Amp.resize(nsamples);
      Time.resize(nsamples, 0.0);
      if (IsOscilloscope)
        for (int j = 0; j < nsamples; j++)
        {
          Amp[j] = 1.0 * static_cast<double>(wd->at(j));
        }
      else
        for (int j = 0; j < nsamples; j++)
        {
          Amp[j] = 1.0 * static_cast<double>(ws->at(j));
        }
      for (int j = 0; j < nsamples; j++)
      {
        Time[j] = j * sampling;
      }
    }

    TH1D *GetTH1(int pm)
    {
      TH1D *h = new TH1D(Form("ADCChannel%i_Event%i", pm, Number), Form("ADCChannel%i_Event%i;Time (s); ADC", pm, Number), nsamples, 0, nsamples * sampling);
      for (int j = 0; j < nsamples; j++)
        h->SetBinContent(j + 1, Amp[j]);
      return h;
    }

    TH1D *GetScintProf(int pm, double ped, int PeakBin, double r1 = 0, double r2 = 0)
    {
      TH1D *wf2 = new TH1D(Form("ADCChannel%i_Event%i_ScintProf2", pm, Number), Form("ADCChannel%i_Event%i;Time (s); ADC", pm, Number), nsamples, 0, nsamples * sampling);
      for (int i = 0; i < nsamples; i++)
      {
        wf2->SetBinContent(i + 1, ped - Amp[i]);
      }
      if (r1 || r2)
        wf2->GetXaxis()->SetRangeUser(r1, r2);
      int binmax = wf2->GetMaximumBin();
      TH1D *wf = new TH1D(Form("ADCChannel%i_Event%i_ScintProf", pm, Number), Form("ADCChannel%i_Event%i;Time (s); ADC", pm, Number), nsamples, 0, nsamples * sampling);
      for (int i = 0; i < nsamples; i++)
        if (i + 1 - binmax + PeakBin > 0 && i + 1 - binmax + PeakBin < nsamples + 1)
        {
          wf->SetBinContent(i + 1 - binmax + PeakBin, wf2->GetBinContent(i + 1));
        }
      delete wf2;
      return wf;
    }

    TH1D *GetScintProfTimeShifted(int pm, double ped, double TimeStamp, int PeakBin)
    {
      TH1D *wf2 = new TH1D(Form("ADCChannel%i_Event%i_ScintProf2", pm, Number), Form("ADCChannel%i_Event%i;Time (s); ADC", pm, Number), nsamples, 0, nsamples * sampling);
      for (int i = 0; i < nsamples; i++)
      {
        wf2->SetBinContent(i + 1, ped - Amp[i]);
      }
      TAxis *xaxis = wf2->GetXaxis();
      Int_t binmax = xaxis->FindBin(TimeStamp);
      TH1D *wf = new TH1D(Form("ADCChannel%i_Event%i_ScintProf", pm, Number), Form("ADCChannel%i_Event%i;Time (s); ADC", pm, Number), nsamples, 0, nsamples * sampling);
      for (int i = 0; i < nsamples; i++)
        if (i + 1 - binmax + PeakBin > 0 && i + 1 - binmax + PeakBin < nsamples + 1)
        {
          wf->SetBinContent(i + 1 - binmax + PeakBin, wf2->GetBinContent(i + 1));
        }
      delete wf2;
      return wf;
    }

    string GetName() { return name; }
  };
  class EventReader_t
  {
    std::vector<WF_t *> WF;
    bool TSReady = false;
    std::vector<double> TS;

  public:
    std::vector<double> Vtime;
    EventReader_t()
    {
      //    WF = new std::vector<WF_t>();
    }
    void SetFile(string file, string name, bool IsOsc = false)
    {
      WF_t *w = new WF_t();
      w->Set(file, name, IsOsc);
      WF.push_back(w);
    }
    void GetEntry(int i)
    {
      for (unsigned int j = 0; j < WF.size(); j++)
      {
        WF[j]->GetEntry(i);
      }
    }
    void FillTS()
    {
      TS.resize(GetNEvents());
      WF[0]->t->GetEntry(0);
      TS[0] = WF[0]->GetTimeStamp();
      double shift = 0;
      for (int i = 1; i < GetNEvents(); i++)
      {
        WF[0]->t->GetEntry(i);
        if (WF[0]->GetTimeStamp() < TS[i - 1])
        {
          shift += 40;
          TS[i] = WF[0]->GetTimeStamp() + shift;
        }
      }
      TSReady = true;
    }

    void Draw(int ev)
    {
      TCanvas *c = new TCanvas("c");
      c->Divide(3, 2);
      GetEntry(ev);
      for (unsigned int i = 0; i < WF.size(); i++)
      {
        c->cd(i + 1);
        TH1D *h = WF[i]->GetTH1(i);
        h->SetTitle(WF[i]->GetName().c_str());
        h->SetName(WF[i]->GetName().c_str());
        h->Draw("HIST");
      }
    }
    int GetNSamples()
    {
      return WF[0]->nsamples;
    }
    int GetNEvents()
    {
      return WF[0]->t->GetEntries();
    }
    std::vector<double> *GetAmp(int ch)
    {
      return &(WF[ch]->Amp);
    }
    std::vector<double> *GetTime()
    {
      return &(WF[0]->Time);
    }
    double GetTimeStamp() { return WF[0]->GetTimeStamp(); }
    double GetPCTimeStamp() { return WF[0]->GetPCTimeStamp(); }
    double GetTimeStamp_nsec() { return WF[0]->GetTimeStamp(); }
    double GetSampling() { return WF[0]->GetSampling(); }

    TH1D *GetWaveform(int pm)
    {
      return WF[pm]->GetTH1(pm);
    }
    TH1D *GetScintProf(int pm, double ped, int PeakBin, double r1 = 0, double r2 = 0)
    {
      return WF[pm]->GetScintProf(pm, ped, PeakBin, r1, r2);
    }
    TH1D *GetScintProfTimeShifted(int pm, double ped, double TimeStamp, int PeakBin)
    {
      return WF[pm]->GetScintProfTimeShifted(pm, ped, TimeStamp, PeakBin);
    }
  };
}
