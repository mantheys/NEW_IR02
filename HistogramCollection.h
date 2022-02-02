#ifndef HistogramCollection_h
#define HistogramCollection_h

#include "HistogramName.h"
namespace ana
{
  class HistogramCollection_t
  {
  public:
    string filein = "";
    string fileout = "";
    std::vector<TH1 *> hist;
    std::vector<ana::HistogramName_t *> histname;
    string drawoption = "HIST";

  public:
    HistogramCollection_t()
    {
      hist.clear();
      histname.clear();
    }

    void Add(TH1 *h, int r, int oc, string sn, double volt, string m, int firstt, string o)
    {
      histname.push_back(new ana::HistogramName_t(r, oc, sn, volt, m, firstt, o));
      h->SetName(histname[histname.size() - 1]->histname.c_str());
      h->SetTitle(histname[histname.size() - 1]->histname.c_str());
      hist.push_back(h);
    }
    void SetDrawOption(string option) { drawoption = option; }
    string GetDrawOption() { return drawoption; }

    void GetFromFile(string f)
    {
      filein = f;
      TFile *f1 = new TFile(filein.c_str(), "READ");
      TIter next(f1->GetListOfKeys());
      TKey *key;
      while ((key = (TKey *)next()))
      {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TH1"))
          continue;
        TH1 *h = (TH1 *)key->ReadObj();
        hist.push_back(h);
        histname.push_back(new ana::HistogramName_t(h->GetName()));
      }
      std::cout << "Loading Histogram Collection of size " << hist.size() << " TH1 from " << filein << std::endl;
    }
    void DumpToFile(string f, string fileoption = "RECREATE")
    {
      fileout = f;
      TFile *dump = new TFile(fileout.c_str(), fileoption.c_str());
      for (auto th : hist)
      {
        th->Write();
      }
      dump->Close();
      std::cout << "Dumped to: " << fileout << std::endl;
    }
    std::vector<TH1 *> GetByOpChannel(int pm)
    {
      std::vector<TH1 *> v;
      for (unsigned int i = 0; i < hist.size(); i++)
      {
        if (histname[i]->GetOpChannel() == pm)
        {
          v.push_back(hist[i]);
        }
      }
      return v;
    }
    std::vector<TH1 *> GetByPMTSN(string SN)
    {
      std::vector<TH1 *> v;
      for (unsigned int i = 0; i < hist.size(); i++)
      {
        if (histname[i]->GetPMTSN() == SN)
        {
          v.push_back(hist[i]);
        }
      }
      return v;
    }
    std::vector<TH1 *> GetByRun(int r)
    {
      std::vector<TH1 *> v;
      for (unsigned int i = 0; i < hist.size(); i++)
      {
        if (histname[i]->GetRun() == r)
        {
          v.push_back(hist[i]);
        }
      }
      return v;
    }
    std::vector<TH1 *> GetByOpChannel64()
    {
      std::vector<TH1 *> v;
      v.resize(64);
      for (unsigned int i = 0; i < hist.size(); i++)
      {
        v[histname[i]->GetOpChannel()] = hist[i];
      }
      return v;
    }

    void Plot(int logoption = 0)
    {
      TCanvas *c = new TCanvas("c");
      float txtsize = 0.06;
      c->DivideSquare(hist.size());
      int canvas = 1;
      TGaxis::SetMaxDigits(3);
      for (auto h : hist)
      {
        c->cd(canvas);
        h->Draw(GetDrawOption().c_str());
        if (logoption == 2 || logoption == 3)
          gPad->SetLogy();
        if (logoption == 1 || logoption == 3)
          gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        h->SetLineColor(2);
        h->GetXaxis()->SetLabelSize(txtsize);
        h->GetXaxis()->SetNdivisions(5);
        h->GetYaxis()->SetLabelSize(txtsize);
        h->GetXaxis()->SetTitleSize(txtsize);
        h->GetYaxis()->SetTitleSize(txtsize);
        c->cd(canvas)->Modified();
        c->cd(canvas)->Update();
        canvas++;
      }
      lets_pause();
    }

    void PlotAll(int logoption = 0)
    {
      TCanvas *c = new TCanvas("c");
      float txtsize = 0.06;
      TGaxis::SetMaxDigits(3);
      int i = 0;
      for (auto h : hist)
      {
        c->cd();
        h->Draw(Form("%s SAME", GetDrawOption().c_str()));
        if (logoption == 2 || logoption == 3)
          gPad->SetLogy();
        if (logoption == 1 || logoption == 3)
          gPad->SetLogx();
        gPad->SetGridx();
        gPad->SetGridy();
        h->SetLineColor(i + 1);
        i++;
        h->GetXaxis()->SetLabelSize(txtsize);
        h->GetXaxis()->SetNdivisions(5);
        h->GetYaxis()->SetLabelSize(txtsize);
        h->GetXaxis()->SetTitleSize(txtsize);
        h->GetYaxis()->SetTitleSize(txtsize);
        c->Modified();
        c->Update();
      }
      lets_pause();
    }
  };
}

#endif
