namespace ana
{
  class NoiseTools_t
  {
  public:
    NoiseTools_t() {} /*
      double Averate(std::vector<double> v)
      {
         double av=0;
         for(auto i : v) av+=i;
         return
      }*/
    void MovingAverageFilter(TH1D *h, int n)
    {
      //     TCanvas * c1= new TCanvas("temp");
      //     h->Draw("HIST");
      // TH1D *h2 = (TH1D*)h->Clone("h2"); h2->Reset();
      std::vector<double> v;
      v.resize(n);
      for (int j = 0; j < n; j++)
        v[j] = h->GetBinContent(j + 1);

      for (int i = 1; i < h->GetSize() - n; i++)
      {
        h->SetBinContent(i, accumulate(v.begin(), v.end(), 0.0) / v.size());
        v.erase(v.begin());
        v.push_back(h->GetBinContent(i + n - 1));
      }
      for (int i = h->GetSize() - n; i < h->GetSize(); i++)
        h->SetBinContent(i, 0);
      // h2->SetLineColor(2);
      // h2->Draw("SAME");
      // lets_pause();
      // h2->Delete();
      // c1->Delete();
    }
  };
}
