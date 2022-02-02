
#ifndef Decon_h
#define Decon_h

namespace ana
{
  class Decon_t
  {
  public:
    std::vector<double> response;
    int threshold = 19;

    Decon_t(std::vector<double> res) : response(res) {}
    TH1D *Deconvolute(TH1D *h)
    {
      TH1D *h2 = (TH1D *)h->Clone("temp");
      TH1D *hres = (TH1D *)h->Clone("hres");
      hres->Reset();
      while (1)
      {
        h2->GetXaxis()->UnZoom();
        int maxbin = h2->FindFirstBinAbove(threshold); // cout << maxbin << endl;
        if (maxbin == -1)
          break;
        double N = h2->GetBinContent(maxbin) / response[0];
        hres->Fill(h2->GetBinCenter(maxbin), N);
        for (unsigned int i = 0; i < response.size(); i++)
          h2->SetBinContent(maxbin + i, h2->GetBinContent(maxbin + i) - N * response[i]);
      }

      h2->Delete();
      return hres;
    }
  };
}

#endif
