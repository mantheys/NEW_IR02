#include "TSystem.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TKey.h"
#include "TMath.h"
#include "TImage.h"
#include "TMinuit.h"
#include "TPolyMarker.h"
#include "TPolyMarker.h"
#include "TGaxis.h"
#include "TNtupleD.h"

// para el mapa de runes de las tomas de medidas
#include "json_run_map.h"

std::vector<double> GetVector(TH1D *hh)
{
  std::vector<double> a;
  for (int i = 1; i < hh->GetSize(); i++)
  {
    a.push_back(hh->GetBinContent(i));
  }
  return a;
}
double VectorDistance(double ped, std::vector<double> pedestales)
{
  double sum = 0;
  for (unsigned int i = 0; i < pedestales.size(); i++)
    sum += TMath::Abs(ped - pedestales[i]);
  return sum;
}

void printpoints(TGraphErrors *tg)
{

  for (int i = 0; i < tg->GetN(); i++)
  {
    double x, y, ye;
    tg->GetPoint(i, x, y);
    ye = tg->GetErrorY(i);
    cout << tg->GetName() << "\t" << i << "\t" << x << "\t" << y << "\t" << ye << endl;
  }
}

struct physvalue
{
  double value;
  double error;
};

physvalue diff(physvalue val1, physvalue val2)
{
  physvalue dif;
  dif.value = val1.value - val2.value;
  dif.error = TMath::Sqrt(val1.error * val1.error + val2.error * val2.error);
  return dif;
}
string discrep(physvalue val1, physvalue val2)
{
  if (val2.value == 0 || val1.value == 0)
    return "";
  physvalue dif = diff(val2, val1);
  return Form("%.1f ( %.1f )", dif.value / dif.error, 100.0 * val2.value / val1.value - 100.0);
}
double getdiscrep(physvalue val1, physvalue val2)
{
  physvalue dif = diff(val2, val1);
  return dif.value / dif.error;
}

string printgain(physvalue val)
{
  if (val.value == 0)
    return "";
  if (val.error > 19e6)
  {
    return Form("%.0f0\u00b1%.0f0", 1e-7 * val.value, 1e-7 * val.error);
  }
  if (val.error > 1.9e6)
  {
    return Form("%.0f\u00b1%.0f", 1e-6 * val.value, 1e-6 * val.error);
  }
  return Form("%.1f\u00b1%.1f", 1e-6 * val.value, 1e-6 * val.error);
}
string printgain2(physvalue val)
{
  if (val.value == 0)
    return "";
  if (val.error > 1.9)
    return Form("%.0f+%.0f", 1e-6 * val.value, 1e-6 * val.error);
  return Form("%.1f+%.1f", 1e-6 * val.value, 1e-6 * val.error);
}

string printphysvalue(physvalue val)
{
  return Form("%.1f\u00b1%.1f", val.value, val.error);
}

void plotdiscrepancies(std::vector<TH1D *> h, std::vector<TH1D *> hh, TGraphErrors *tg, std::vector<TGraphErrors *> tg2)
{

  std::map<int, physvalue> map1;
  std::vector<std::map<int, physvalue>> map2;
  map2.resize(tg2.size());

  double x, y, ye;
  for (int i = 0; i < tg->GetN(); i++)
  {
    tg->GetPoint(i, x, y);
    ye = tg->GetErrorY(i);
    map1[(int)x].value = y;
    map1[(int)x].error = ye;
  }
  for (unsigned int k = 0; k < tg2.size(); k++)
    for (int i = 0; i < tg2[k]->GetN(); i++)
    {
      tg2[k]->GetPoint(i, x, y);
      ye = tg2[k]->GetErrorY(i);
      map2[k][(int)x].value = y;
      map2[k][(int)x].error = ye;
    }

  for (unsigned int k = 0; k < tg2.size(); k++)
    for (std::map<int, physvalue>::iterator it = map2[k].begin(); it != map2[k].end(); ++it)
    {
      // std::cout << it->first << " => " << it->second.value << '\n';
      physvalue empty = {0, 0};
      std::map<int, physvalue>::iterator it1;
      it1 = map1.find(it->first);
      if (it1 == map1.end())
        map1.insert(std::pair<int, physvalue>(it->first, empty));
    }

  cout << Form("Voltage\tGain_%s(10^6)", tg->GetName());
  for (unsigned int k = 0; k < tg2.size(); k++)
    cout << Form("\t Gain_%s(10^6)\t\u0394(\u03C3)", tg2[k]->GetName());
  cout << endl;

  for (std::map<int, physvalue>::iterator it = map1.begin(); it != map1.end(); ++it)
  {

    for (unsigned int k = 0; k < tg2.size(); k++)
      if (it->second.value != 0 && map2[k][it->first].value != 0)
        h[k]->Fill(getdiscrep(it->second, map2[k][it->first]));

    for (unsigned int k = 0; k < tg2.size(); k++)
      if (it->second.value != 0 && map2[k][it->first].value != 0)
        hh[k]->Fill(100.0 * map2[k][it->first].value / it->second.value - 100);
    std::cout << it->first << "\t" << printgain(it->second);
    for (unsigned int k = 0; k < tg2.size(); k++)
      std::cout << "\t" << printgain(map2[k][it->first]) << "\t" << discrep(it->second, map2[k][it->first]);
    cout << endl;
  }
}

int autorangemin(TH1 *myhist)
{
  for (Int_t i = 0; i < myhist->GetNbinsX(); i++)
    if (myhist->GetBinContent(i) > 1)
      return i;
  return -1;
}
int autorangemax(TH1 *myhist)
{
  for (Int_t i = myhist->GetNbinsX(); i > 0; i--)
    if (myhist->GetBinContent(i) > 1)
      return i;
  return -1;
}

void setlogscale(TCanvas *c1, TMultiGraph *mg, Double_t xlow = 0, Double_t xup = 0, Double_t LabelStep = 100)
{

  double ylow = ((TGraphErrors *)mg->GetListOfGraphs()->At(0))->GetY()[0];
  double yup = ((TGraphErrors *)mg->GetListOfGraphs()->At(0))->GetY()[0];
  for (int i = 0; i < mg->GetListOfGraphs()->GetSize(); i++)
    if (((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetY()[0] < ylow)
      ylow = ((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetY()[0];
  for (int i = 0; i < mg->GetListOfGraphs()->GetSize(); i++)
    if (((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetY()[((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetN() - 1] > yup)
      yup = ((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetY()[((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetN() - 1];

  ylow = 0.8 * ylow;
  yup = 2 * yup;
  std::cout << "setlogscale::Setting Ylevels: " << ylow << " " << yup << std::endl;

  mg->GetYaxis()->SetRangeUser(ylow, yup);

  xlow = ((TGraphErrors *)mg->GetListOfGraphs()->At(0))->GetX()[0];
  xup = ((TGraphErrors *)mg->GetListOfGraphs()->At(0))->GetX()[0];
  for (int i = 0; i < mg->GetListOfGraphs()->GetSize(); i++)
    if (((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetX()[0] < xlow)
      xlow = ((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetX()[0];
  for (int i = 0; i < mg->GetListOfGraphs()->GetSize(); i++)
    if (((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetX()[((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetN() - 1] > xup)
      xup = ((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetX()[((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetN() - 1];
  xlow = xlow - LabelStep;
  xup = xup + LabelStep;
  // c1->SetGrid();
  c1->SetLogy();
  c1->SetLogx();
  gPad->SetLogx();
  gPad->SetLogy();
  TStyle *plain = new TStyle("Plain", "Plain Style (no colors/fill areas)");
  plain->SetCanvasBorderMode(0);
  plain->SetPadBorderMode(0);
  plain->SetPadColor(0);
  plain->SetCanvasColor(0);
  plain->SetStatColor(0);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0000); // no mostrar cuadro fit (0011, etc)

  // Lines needed to have the value in Log Scale
  // xlow   = 1000;
  // xup    = 1800;
  Int_t nbinsx = Int_t(xup - xlow);
  // cout << " ho visto xmin" << endl;
  Int_t step_normal = LabelStep;
  Int_t factor_large = 1;
  //  Int_t    factor_large    = ;
  Double_t fraction_normal = 0.015;
  Double_t fraction_large = 0.035;

  Double_t uxmin = c1->GetUxmin();
  Double_t uxmax = c1->GetUxmax();
  Double_t uymin = c1->GetUymin();
  Double_t uymax = c1->GetUymax();

  // x-axis ticks
  //----------------------------------------------------------------------------
  if (c1->GetLogx())
  {

    Int_t ixlow = Int_t(xlow);
    Int_t ixup = Int_t(xup);
    Int_t step_large = factor_large * step_normal;

    mg->GetXaxis()->SetNdivisions(0);

    Double_t ymin = (c1->GetLogy()) ? pow(10, uymin) : uymin;
    Double_t ymax = (c1->GetLogy()) ? pow(10, uymax) : uymax;

    Double_t bottom_tick_normal = fraction_normal * (uymax - uymin);
    Double_t bottom_tick_large = fraction_large * (uymax - uymin);

    Double_t top_tick_normal = fraction_normal * (uymax - uymin);
    Double_t top_tick_large = fraction_large * (uymax - uymin);

    if (c1->GetLogy())
    {

      bottom_tick_normal = pow(10, uymin + bottom_tick_normal) - pow(10, uymin);
      bottom_tick_large = pow(10, uymin + bottom_tick_large) - pow(10, uymin);

      TLine tick;

      // for (Int_t i=xlow; i<=xup; i+=100) {
      for (Int_t i = xlow + LabelStep; i < xup; i += LabelStep)
      {

        Double_t xx = i;

        tick.DrawLine(xx, ymin, xx, ymin + ((i - ixlow + 50) % step_large == 0 ? bottom_tick_large : bottom_tick_normal));
      }

      // x-axis labels
      //--------------------------------------------------------------------------

      Double_t ylatex = ymin - bottom_tick_normal;

      if ((ixup - ixlow) % step_large >= step_normal)
      {

        TLatex *latex = new TLatex(xup, ylatex, Form("%.0f", xup));

        latex->SetTextAlign(23);
        latex->SetTextFont(42);
        latex->SetTextSize(0.035);
        latex->Draw("same");
      }

      for (Int_t i = ixlow + LabelStep; i < ixup; i += step_large)
      {

        Double_t xx = i;
        TLatex *latex = new TLatex(xx, ylatex, Form("%.0f", xx));

        latex->SetTextAlign(23);
        latex->SetTextFont(42);
        latex->SetTextSize(0.035);
        latex->Draw("same");
      }
    }
  } // end initial if
}

void setlogscale(TVirtualPad *c1, TMultiGraph *mg, Double_t xlow = 0, Double_t xup = 0, Double_t LabelStep = 100)
{

  double ylow = ((TGraphErrors *)mg->GetListOfGraphs()->At(0))->GetY()[0];
  double yup = ((TGraphErrors *)mg->GetListOfGraphs()->At(0))->GetY()[0];
  // for(unsigned int i=0;i<mg->GetListOfGraphs()->GetSize();i++) if (((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetY()[0]<ylow) ylow=((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetY()[0];
  // for(unsigned int i=0;i<mg->GetListOfGraphs()->GetSize();i++)if (((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetY()[((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetN()-1]>yup) yup=((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetY()[((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetN()-1];

  for (int i = 0; i < mg->GetListOfGraphs()->GetSize(); i++)
    for (int j = 0; j < ((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetN(); j++)
      if (((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetY()[j] < ylow)
        ylow = ((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetY()[j];
  for (int i = 0; i < mg->GetListOfGraphs()->GetSize(); i++)
    for (int j = 0; j < ((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetN(); j++)
      if (((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetY()[j] > yup)
        yup = ((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetY()[j];

  ylow = 0.8 * ylow;
  yup = 2 * yup;
  std::cout << "setlogscale::Setting Ylevels: " << ylow << " " << yup << std::endl;

  mg->GetYaxis()->SetRangeUser(ylow, yup);
  mg->GetYaxis()->SetLabelSize(0.07);

  xlow = ((TGraphErrors *)mg->GetListOfGraphs()->At(0))->GetX()[0];
  xup = ((TGraphErrors *)mg->GetListOfGraphs()->At(0))->GetX()[0];
  for (int i = 0; i < mg->GetListOfGraphs()->GetSize(); i++)
    if (((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetX()[0] < xlow)
      xlow = ((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetX()[0];
  for (int i = 0; i < mg->GetListOfGraphs()->GetSize(); i++)
    if (((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetX()[((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetN() - 1] > xup)
      xup = ((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetX()[((TGraphErrors *)mg->GetListOfGraphs()->At(i))->GetN() - 1];
  xlow = xlow - LabelStep;
  xup = xup + LabelStep;
  // c1->SetGrid();
  c1->SetLogy();
  c1->SetLogx();
  gPad->SetLogx();
  gPad->SetLogy();
  TStyle *plain = new TStyle("Plain", "Plain Style (no colors/fill areas)");
  plain->SetCanvasBorderMode(0);
  plain->SetPadBorderMode(0);
  plain->SetPadColor(0);
  plain->SetCanvasColor(0);
  plain->SetStatColor(0);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0000); // no mostrar cuadro fit (0011, etc)

  // Lines needed to have the value in Log Scale
  // xlow   = 1000;
  // xup    = 1800;
  Int_t nbinsx = Int_t(xup - xlow);
  // cout << " ho visto xmin" << endl;
  Int_t step_normal = LabelStep;
  Int_t factor_large = 1;
  //  Int_t    factor_large    = ;
  Double_t fraction_normal = 0.015;
  Double_t fraction_large = 0.035;

  Double_t uxmin = c1->GetUxmin();
  Double_t uxmax = c1->GetUxmax();
  Double_t uymin = c1->GetUymin();
  Double_t uymax = c1->GetUymax();

  // x-axis ticks
  //----------------------------------------------------------------------------
  if (c1->GetLogx())
  {

    Int_t ixlow = Int_t(xlow);
    Int_t ixup = Int_t(xup);
    Int_t step_large = factor_large * step_normal;

    mg->GetXaxis()->SetNdivisions(0);

    Double_t ymin = (c1->GetLogy()) ? pow(10, uymin) : uymin;
    Double_t ymax = (c1->GetLogy()) ? pow(10, uymax) : uymax;

    Double_t bottom_tick_normal = fraction_normal * (uymax - uymin);
    Double_t bottom_tick_large = fraction_large * (uymax - uymin);

    Double_t top_tick_normal = fraction_normal * (uymax - uymin);
    Double_t top_tick_large = fraction_large * (uymax - uymin);

    if (c1->GetLogy())
    {

      bottom_tick_normal = pow(10, uymin + bottom_tick_normal) - pow(10, uymin);
      bottom_tick_large = pow(10, uymin + bottom_tick_large) - pow(10, uymin);

      TLine tick;

      // for (Int_t i=xlow; i<=xup; i+=100) {
      for (Int_t i = xlow + LabelStep; i < xup; i += LabelStep)
      {

        Double_t xx = i;

        tick.DrawLine(xx, ymin, xx, ymin + ((i - ixlow + 50) % step_large == 0 ? bottom_tick_large : bottom_tick_normal));
      }

      // x-axis labels
      //--------------------------------------------------------------------------

      Double_t ylatex = ymin - bottom_tick_normal;

      if ((ixup - ixlow) % step_large >= step_normal)
      {

        TLatex *latex = new TLatex(xup, ylatex, Form("%.0f", xup));

        latex->SetTextAlign(23);
        latex->SetTextFont(42);
        latex->SetTextSize(0.07);
        latex->Draw("same");
      }

      for (Int_t i = ixlow + LabelStep; i < ixup; i += step_large)
      {

        Double_t xx = i;
        TLatex *latex = new TLatex(xx, ylatex, Form("%.0f", xx));

        latex->SetTextAlign(23);
        latex->SetTextFont(42);
        latex->SetTextSize(0.07);
        latex->Draw("same");
      }
    }
  } // end initial if
}

void setlogscale(TCanvas *c1, TGraphErrors *mg, Double_t xlow = 0, Double_t xup = 0, Double_t LabelStep = 100)
{

  double ylow = mg->GetY()[0];
  double yup = mg->GetY()[0];
  for (int i = 0; i < mg->GetN(); i++)
    if (mg->GetY()[i] < ylow)
      ylow = mg->GetY()[i];
  for (int i = 0; i < mg->GetN(); i++)
    if (mg->GetY()[i] > yup)
      yup = mg->GetY()[i];

  ylow = 0.8 * ylow;
  yup = 2 * yup;
  // mg->GetYaxis()->SetRangeUser(ylow,yup);
  std::cout << "setlogscale1::Setting Ylevels: " << ylow << " " << yup << std::endl;
  xlow = mg->GetX()[0];
  xup = mg->GetX()[0];
  for (int i = 0; i < mg->GetN(); i++)
    if (mg->GetX()[i] < xlow)
      xlow = mg->GetX()[i];
  for (int i = 0; i < mg->GetN(); i++)
    if (mg->GetX()[i] > xup)
      xup = mg->GetX()[i];

  xlow = LabelStep * ((int)(xlow / LabelStep - 1));
  xup = LabelStep * ((int)(xup / LabelStep + 1));
  std::cout << "setlogscale1::Setting XLevels: " << xlow << " " << xup << std::endl;

  // c1->SetGrid();
  c1->SetLogy();
  c1->SetLogx();
  gPad->SetLogx();
  gPad->SetLogy();
  TStyle *plain = new TStyle("Plain", "Plain Style (no colors/fill areas)");
  plain->SetCanvasBorderMode(0);
  plain->SetPadBorderMode(0);
  plain->SetPadColor(0);
  plain->SetCanvasColor(0);
  plain->SetStatColor(0);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0000); // no mostrar cuadro fit (0011, etc)

  // Lines needed to have the value in Log Scale
  xlow = 1000;
  xup = 1800;
  Int_t nbinsx = Int_t(xup - xlow);
  cout << "nbinsx " << nbinsx << endl;
  // cout << " ho visto xmin" << endl;
  Int_t step_normal = LabelStep;
  Int_t factor_large = 1;
  //  Int_t    factor_large    = ;
  Double_t fraction_normal = 0.015;
  Double_t fraction_large = 0.035;

  Double_t uxmin = c1->GetUxmin();
  Double_t uxmax = c1->GetUxmax();
  Double_t uymin = c1->GetUymin();
  Double_t uymax = c1->GetUymax();

  // x-axis ticks
  //----------------------------------------------------------------------------
  if (c1->GetLogx())
  {
    cout << "in X " << endl;

    Int_t ixlow = Int_t(xlow);
    Int_t ixup = Int_t(xup);
    Int_t step_large = factor_large * step_normal;

    mg->GetXaxis()->SetNdivisions(0);

    Double_t ymin = (c1->GetLogy()) ? pow(10, uymin) : uymin;
    Double_t ymax = (c1->GetLogy()) ? pow(10, uymax) : uymax;

    Double_t bottom_tick_normal = fraction_normal * (uymax - uymin);
    Double_t bottom_tick_large = fraction_large * (uymax - uymin);

    Double_t top_tick_normal = fraction_normal * (uymax - uymin);
    Double_t top_tick_large = fraction_large * (uymax - uymin);

    if (c1->GetLogy())
    {

      bottom_tick_normal = pow(10, uymin + bottom_tick_normal) - pow(10, uymin);
      bottom_tick_large = pow(10, uymin + bottom_tick_large) - pow(10, uymin);

      TLine tick;

      // for (Int_t i=xlow; i<=xup; i+=100) {
      for (Int_t i = xlow + LabelStep; i < xup; i += LabelStep)
      {

        Double_t xx = i;

        tick.DrawLine(xx, ymin, xx, ymin + ((i - ixlow + 50) % step_large == 0 ? bottom_tick_large : bottom_tick_normal));
        cout << "Draw line " << endl;
      }

      // x-axis labels
      //--------------------------------------------------------------------------

      Double_t ylatex = ymin - bottom_tick_normal;

      if ((ixup - ixlow) % step_large >= step_normal)
      {

        TLatex *latex = new TLatex(xup, ylatex, Form("%.0f", xup));

        latex->SetTextAlign(23);
        latex->SetTextFont(42);
        latex->SetTextSize(0.035);
        latex->Draw("same");
      }

      for (Int_t i = ixlow + LabelStep; i < ixup; i += step_large)
      {

        Double_t xx = i;
        TLatex *latex = new TLatex(xx, ylatex, Form("%.0f", xx));

        latex->SetTextAlign(23);
        latex->SetTextFont(42);
        latex->SetTextSize(0.035);
        latex->Draw("same");
      }
    }
  } // end initial if
}

void setlogscale(TVirtualPad *c1, TGraphErrors *mg, Double_t xlow = 0, Double_t xup = 0, Double_t LabelStep = 100)
{
  lets_pause();

  double ylow = mg->GetY()[0];
  double yup = mg->GetY()[0];
  for (int i = 0; i < mg->GetN(); i++)
    if (mg->GetY()[i] < ylow)
      ylow = mg->GetY()[i];
  for (int i = 0; i < mg->GetN(); i++)
    if (mg->GetY()[i] > yup)
      yup = mg->GetY()[i];

  ylow = 0.7 * ylow;
  yup = 1.5 * yup;
  std::cout << "setlogscale2::Setting YLevels: " << ylow << " " << yup << std::endl;

  xlow = mg->GetX()[0];
  xup = mg->GetX()[0];
  for (int i = 0; i < mg->GetN(); i++)
    if (mg->GetX()[i] < xlow)
      xlow = mg->GetX()[i];
  for (int i = 0; i < mg->GetN(); i++)
    if (mg->GetX()[i] > xup)
      xup = mg->GetX()[i];

  xlow = LabelStep * ((int)(xlow / LabelStep - 1));
  xup = LabelStep * ((int)(xup / LabelStep + 1));

  // lets_pause();
  c1->SetLogy();
  c1->SetLogx();
  gPad->SetLogx();
  gPad->SetLogy();

  mg->GetYaxis()->SetRangeUser(ylow, yup);
  // mg->GetYaxis()->UnZoom();
  std::cout << "setlogscale2::Setting XLevels: " << xlow << " " << xup << std::endl; // lets_pause();
  // c1->SetGrid();
  TStyle *plain = new TStyle("Plain", "Plain Style (no colors/fill areas)");
  plain->SetCanvasBorderMode(0);
  plain->SetPadBorderMode(0);
  plain->SetPadColor(0);
  plain->SetCanvasColor(0);
  plain->SetStatColor(0);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0000); // no mostrar cuadro fit (0011, etc)

  // Lines needed to have the value in Log Scale
  Int_t nbinsx = Int_t(xup - xlow);
  cout << "nbinsx " << nbinsx << endl;
  // cout << " ho visto xmin" << endl;
  Int_t step_normal = LabelStep;
  Int_t factor_large = 1;
  //  Int_t    factor_large    = ;
  Double_t fraction_normal = 0.015;
  Double_t fraction_large = 0.035;

  Double_t uxmin = c1->GetUxmin();
  Double_t uxmax = c1->GetUxmax();
  Double_t uymin = c1->GetUymin();
  Double_t uymax = c1->GetUymax();

  // x-axis ticks
  //----------------------------------------------------------------------------
  if (c1->GetLogx())
  {
    cout << "in X " << endl;

    Int_t ixlow = Int_t(xlow);
    Int_t ixup = Int_t(xup);
    Int_t step_large = factor_large * step_normal;

    mg->GetXaxis()->SetNdivisions(0);

    Double_t ymin = (c1->GetLogy()) ? pow(10, uymin) : uymin;
    Double_t ymax = (c1->GetLogy()) ? pow(10, uymax) : uymax;

    Double_t bottom_tick_normal = fraction_normal * (uymax - uymin);
    Double_t bottom_tick_large = fraction_large * (uymax - uymin);

    Double_t top_tick_normal = fraction_normal * (uymax - uymin);
    Double_t top_tick_large = fraction_large * (uymax - uymin);

    if (c1->GetLogy())
    {

      bottom_tick_normal = pow(10, uymin + bottom_tick_normal) - pow(10, uymin);
      bottom_tick_large = pow(10, uymin + bottom_tick_large) - pow(10, uymin);

      TLine tick;

      // for (Int_t i=xlow; i<=xup; i+=100) {
      for (Int_t i = xlow + LabelStep; i < xup; i += LabelStep)
      {

        Double_t xx = i;

        tick.DrawLine(xx, ymin, xx, ymin + ((i - ixlow + 50) % step_large == 0 ? bottom_tick_large : bottom_tick_normal));
        cout << "Draw line " << endl;
      }

      // x-axis labels
      //--------------------------------------------------------------------------

      Double_t ylatex = ymin - bottom_tick_normal;

      if ((ixup - ixlow) % step_large >= step_normal)
      {

        TLatex *latex = new TLatex(xup, ylatex, Form("%.0f", xup));

        latex->SetTextAlign(23);
        latex->SetTextFont(42);
        latex->SetTextSize(0.035);
        latex->Draw("same");
      }

      for (Int_t i = ixlow + LabelStep; i < ixup; i += step_large)
      {

        Double_t xx = i;
        TLatex *latex = new TLatex(xx, ylatex, Form("%.0f", xx));

        latex->SetTextAlign(23);
        latex->SetTextFont(42);
        latex->SetTextSize(0.035);
        latex->Draw("same");
      }
    }
  } // end initial if
}

bool lets_check_overflow(string histname, string histfile)
{

  TFile data1(histfile.c_str(), "READ"); //,"RECREATE");

  TH1F *hist = (TH1F *)data1.Get(histname.c_str());

  TAxis *xaxis = hist->GetXaxis();
  float vh[4097];
  float vx[4097];
  float adcfirst = 0;
  float adclast = 0;
  int nloop = 0;
  float vhmax = 0;
  float bmx = 0;

  float second_max = 0;
  float h_second_max = 0;
  float bin_mean;

  float ntot = 0;
  for (unsigned int ic = 1; ic <= 4096; ic++)
  {
    //         for(unsigned int ic=1; ic<=2048; ic++){
    vh[ic] = hist->GetBinContent(ic); // altura del bin
    vx[ic] = xaxis->GetBinCenter(ic); // valor x del bin
    ntot += vh[ic];

    if (vh[ic] > vhmax)
    {
      vhmax = vh[ic];
      bmx = vx[ic];
    }
  }
  float trigger = 2; // vhmax*0.005;
  float ini_second_max = 4096;

  for (unsigned int ic = 1; ic <= 4096; ic++)
  {

    if (vh[ic] > trigger && nloop == 0)
    {
      nloop++;
      // First bin with contents
      adcfirst = vx[ic]; // primer bin con señal
      float aux = 2 * bmx - adcfirst;
      ini_second_max = (int)aux;

      float aux2 = 2.05 * bmx - adcfirst;
      h_second_max = vh[(int)aux2];
      second_max = (int)aux2;

      // cout << "ini second max " << ini_second_max << "bmx y adc" << bmx << " " << adcfirst<< endl;
    }
    // Maximum of the pedestal
    if (ic > ini_second_max)
    {

      if (vh[ic] > h_second_max)
      {
        h_second_max = vh[ic];
        second_max = vx[ic];
      }

      //}
    }
  }
  data1.Close();

  if (bmx == 4095 || second_max == 4095)
  {
    cout << "   ------------------------  OVERFLOW -in-" << histname << "----------------------" << endl;
    return true;
  }
  else
  {
    return false;
  }
}

void print_pdf(int i, int tanda, string voltage_str, float voltage, float maxvolt, string histfile, string backupfile)
{
  string canvasname = "SPE HI-RES V" + voltage_str;
  TFile data1(histfile.c_str(), "READ"); //,"RECREATE");
  TH1F *hist;

  TCanvas *c1 = new TCanvas(canvasname.c_str(), canvasname.c_str());
  c1->SetTitle(canvasname.c_str());
  c1->Divide(3, 2);

  for (unsigned int ch = 1; ch <= 5; ch++)
  {
    string histn = "pmt" + std::to_string(ch + tanda * 5) + "_hi" + voltage_str;
    cout << "Abro histograma " << histn << endl;
    hist = (TH1F *)data1.Get(histn.c_str());
    hist->SetTitle(histn.c_str());
    hist->GetXaxis()->SetRange(autorangemin(hist) * 0.9, autorangemax(hist) * 1.1);
    c1->cd(ch);
    hist->Draw("HIST");
    gPad->SetLogy(); // pm2->Draw();
    c1->cd(ch)->Modified();
    c1->cd(ch)->Update();
  }
  canvasname = canvasname + ".png";
  if (i == 0)
    c1->Print(Form("%s %s %s", "results/", backupfile.c_str(), "/HIRES.pdf("), "pdf");
  else
    c1->Print(Form("%s %s %s", "results/", backupfile.c_str(), "/HIRES.pdf"), "pdf");
  if (voltage == maxvolt)
    c1->Print(Form("%s %s %s", "results/", backupfile.c_str(), "/HIRES.pdf)"), "pdf");

  cout << "Vamos a los low res " << endl;

  // lets_pause();
  // c1->Print(canvasname.c_str());
  c1->Close();

  canvasname = "SPE LOW-RES V" + voltage_str;
  TCanvas *c = new TCanvas(canvasname.c_str(), canvasname.c_str());
  c->SetTitle(canvasname.c_str());
  c->Divide(3, 2);

  for (unsigned int ch = 1; ch <= 5; ch++)
  {
    string histn = "pmt" + std::to_string(ch + 5 * tanda);
    histn = histn + "_low" + voltage_str;
    cout << "Abro histograma " << histn << endl;
    hist = (TH1F *)data1.Get(histn.c_str());
    cout << "minimobin " << autorangemin(hist) << endl;
    cout << "maxbin " << autorangemax(hist) << endl;
    hist->SetTitle(histn.c_str());
    hist->GetXaxis()->SetRange(autorangemin(hist) * 0.9, autorangemax(hist) * 1.1);
    c->cd(ch);
    hist->Draw("HIST");
    gPad->SetLogy(); // pm2->Draw();
    c->cd(ch)->Modified();
    c->cd(ch)->Update();
  }
  canvasname = canvasname + ".png";
  if (i == 0)
    c->Print(Form("%s %s %s", "results/", backupfile.c_str(), "/LOWRES.pdf("), "pdf");
  else
    c->Print(Form("%s %s %s", "results/", backupfile.c_str(), "/LOWRES.pdf"), "pdf");
  if (voltage == maxvolt)
    c->Print(Form("%s %s %s", "results/", backupfile.c_str(), "/LOWRES.pdf)"), "pdf");

  // c->Print(canvasname.c_str());
  c->Close();
  data1.Close();
}

void print_fitter_pervoltage(int i, int tanda, string voltage_str, float voltage, float maxvolt, std::vector<string> plots[5], string fecha, int pausa, string backupfile, string mode)
{
  cout << "-------------START-print_fitter_pervoltage---------------" << endl;
  string canvasname = "SPE FIT - Double gaussian method - " + voltage_str + "V - " + fecha;
  TCanvas *c1 = new TCanvas("c1", canvasname.c_str(), 1200, 800);
  c1->SetTitle(canvasname.c_str());
  c1->Divide(3, 2);
  TImage *_img[5];

  for (unsigned int ch = 0; ch < 5; ch++)
  {

    cout << "VEAMOS MI VECTR " << i << " " << plots[ch][i] << endl;

    _img[ch] = TImage::Open(plots[ch][i].c_str());
    c1->cd(ch + 1);
    _img[ch]->Draw("xxx");
    _img[ch]->SetEditable(kTRUE);
    c1->cd(ch + 1)->Modified();
    c1->cd(ch + 1)->Update();
  }

  canvasname = "results/" + backupfile + "/" + mode + voltage_str + ".png";
  c1->Print(canvasname.c_str(), "png");
  // if(i==0) c1->Print("results/SPEperV.pdf(");
  // else c1->Print("results/SPEperV.pdf");
  // if (voltage==maxvolt){ c1->Print("results/SPEperV.pdf)");cout  << "ARCHIVO CERRADO"  << endl;}
  if (pausa == 1)
    lets_pause();
  c1->Delete();
  // lets_pause();
  cout << "--------------print_fitter_pervoltage-END--------------" << endl;
}

TH1F *quickrebin(TH1F *hist, string histn)
{
  TH1F *hnew;
  cout << autorangemax(hist) << " " << autorangemin(hist) << endl;

  if (autorangemax(hist) - autorangemin(hist) > 3000)
  {
    cout << histn << " rebinned 1:32" << endl;
    hnew = dynamic_cast<TH1F *>(hist->Rebin(32, histn.c_str()));
  }
  else
  {
    if (autorangemax(hist) - autorangemin(hist) > 2000)
    {
      cout << histn << " rebinned 1:16" << endl;
      hnew = dynamic_cast<TH1F *>(hist->Rebin(16, histn.c_str()));
    }
    else
    {
      if (autorangemax(hist) - autorangemin(hist) > 1000)
      {
        cout << histn << " rebinned 1:8" << endl;
        hnew = dynamic_cast<TH1F *>(hist->Rebin(8, histn.c_str()));
      }
      else
      {
        if (autorangemax(hist) - autorangemin(hist) > 500)
        {
          cout << histn << " rebinned 1:4" << endl;
          hnew = dynamic_cast<TH1F *>(hist->Rebin(4, histn.c_str()));
        }
        else
        {

          if (autorangemax(hist) - autorangemin(hist) > 250)
          {
            cout << histn << " rebinned 1:2" << endl;
            hnew = dynamic_cast<TH1F *>(hist->Rebin(2, histn.c_str()));
          }
          else
          {
            hnew = hist;
          }
        }
      }
    }
  }

  return hnew;
}

/*
void fitandprintgains(string round, int tanda,std::vector<Double_t>  xx[5],std::vector<Double_t>  yy[5],std::vector<Double_t>  xxe[5],std::vector<Double_t>  yye[5],string fecha_str, ofstream &ofile2)
{


  TCanvas *c3 = new TCanvas("c3","c3",1200, 800);
     //c3->Divide(3,2);

  TMultiGraph * mg1 = new TMultiGraph("mg1","");


  std::vector<TGraphErrors*> _graph;
  std::vector<TF1*> _myfit;
  std::vector<bool> _fit_successful;


  std::vector<Double_t> A, errorA, B, errorB, vv, vv9 ,chi2, ndf, covarianzaAB, fA, fB, fB9, errorV, errorV9, errorVtilde, errorVtilde9;

    Double_t gainn=1.e7;
    Double_t gain9=1.e9; //<-- nominal gain



  for(unsigned int i=0;i<5;i++)
  {
     cout << "Creo TGraph " << i << endl;

       _graph.push_back( new TGraphErrors(xx[i].size(),&(xx[i][0]),&(yy[i][0]),&(xxe[i][0]),&(yye[i][0])));
     _fit_successful.push_back(false);

     _graph[i]->SetName("gr1");
    int auxiliar=i+1+tanda*5;
     _graph[i]->SetTitle(Form("%s%i","PMT",auxiliar));
     _graph[i]->SetMarkerStyle(20+i);
     _graph[i]->SetMarkerColor(2+i);
     if (i==0) _graph[i]->SetDrawOption("AP");
     else _graph[i]->SetDrawOption("P");
     _graph[i]->Draw();
  //lets_pause();

     _graph[i]->SetFillStyle(0);

    cout << "Creo myfit " << i << endl;

    _myfit.push_back(new TF1(Form("%s %i","myfit_",i),"pow(10,[0])*pow(x,[1])", 1100, 1800));  //<-- RANGE!
    _myfit[i]->SetParName(0,"A");
    _myfit[i]->SetParName(1,"B");
    _myfit[i]->SetLineStyle(4);
    _myfit[i]->SetLineWidth(3);
    _myfit[i]->SetLineColor(i+2);


     cout << "Ajusto! " << endl;


    TFitResultPtr r = _graph[i]->Fit(Form("%s %i","myfit_",i),"QESV");
    TString estado=gMinuit->fCstatu;
    //float covarianzaAB;
    cout << "resultado ajuste: " << int(r) << " " << estado << endl;
     //lets_pause();
    if (int(r)==0  && estado=="SUCCESSFUL"){
      cout << _fit_successful[i] << " "<< i << endl;
      _fit_successful[i]=true;
      cout << _fit_successful[i] << " "<< i << endl << endl;
      //TMatrixDSym cov = r->GetCovarianceMatrix();
      //cov.Print();
      //cout << cov.Print();
      //cout << "numero columnas" << cov.GetNrows() << endl;//" " <<   cov.GetSub(1,2)<< endl;
      covarianzaAB.push_back(r->GetCovarianceMatrix()[0][1]);
      cout << " eerror ferist par " << covarianzaAB[i] <<endl;
    }
    else
    {
      covarianzaAB.push_back(0.0);
    }
    //lets_pause();

    gStyle->SetOptFit(0000);

    A.push_back(_myfit[i]->GetParameter(0));
    errorA.push_back(_myfit[i]->GetParError(0));
    B.push_back(_myfit[i]->GetParameter(1));
    errorB.push_back(_myfit[i]->GetParError(1));
    chi2.push_back(_myfit[i]->GetChisquare());
    ndf.push_back(_myfit[i]->GetNDF());

    fA.push_back(-1/B[i]);
    fB.push_back(-(TMath::Log10(gainn)-A[i])/(B[i]*B[i]));
    fB9.push_back(-(TMath::Log10(gain9)-A[i])/(B[i]*B[i]));

     errorVtilde.push_back(TMath::Sqrt(fA[i]*fA[i]*errorA[i]*errorA[i]+fB[i]*fB[i]*errorB[i]*errorB[i]+2*fA[i]*fB[i]*covarianzaAB[i]));
     errorVtilde9.push_back(TMath::Sqrt(fA[i]*fA[i]*errorA[i]*errorA[i]+fB9[i]*fB9[i]*errorB[i]*errorB[i]+2*fA[i]*fB9[i]*covarianzaAB[i]));

    //TMatrixDSym cov = _myfit[i]->GetCovarianceMatrix();
    //->GetChisquare();
    //--------------------------------------------------
    //             HV for nominal gain
    //--------------------------------------------------
    //G=10^A*V^B
    vv.push_back(pow(gainn/(pow(10,A[i])),1/B[i]));
    vv9.push_back(pow(gain9/(pow(10,A[i])),1/B[i]));
    errorV.push_back(errorVtilde[i]*vv[i]/TMath::Log10(TMath::E()));
    errorV9.push_back(errorVtilde9[i]*vv9[i]/TMath::Log10(TMath::E()));

    cout<<" "<<endl;
    cout<<"HV(gain="<<gainn<<")"<<"="<<vv[i]<<" "<<"V"<<endl;

    mg1-> Add(_graph[i]);
    cout << "Acabo fit " <<endl;
    //ofile2 << "Pos,Date,Analysis,log_10A,err,14k,Err,V_G=1e7,V_G=1e9" << endl;
    ofile2 << auxiliar<< "," << fecha_str<< ",jose_autogaussian,"<<Form("%.6lf",A[i])<< ","<<Form("%.6lf",errorA[i])<<","<<Form("%.6lf",B[i])<<","<<Form("%.6lf",errorB[i])<<","<<covarianzaAB[i]<<","<<Form("%.0lf",vv[i])<<","<<Form("%.6lf",errorV[i])<<","<<Form("%.0lf",vv9[i])<<","<<Form("%.6lf",errorV9[i]) <<","<<fA[i]<<","<<fB[i]<<","<<fB9[i]<<","<<errorVtilde[i]<<","<<errorVtilde9[i]<<endl;

  }

  mg1->Draw("AP");

  float x1=0.02;
  float y1=0.5;
  float x2=x1+0.4;//anchura
  float y2=0.97; //altura

     TLegend *leg = new TLegend(x1,y1,x2,y2);
     //leg->SetHeader(Form("%s %i %s %i %s","PMTs from ",1+tanda*5," to ",5+tanda*5," @ Room Temperature (G=10^{A}V^{B})"));   //<-- change
  leg->SetHeader(Form("%s %i %s %i %s %s %s","PMTs from ",1+tanda*5," to ",5+tanda*5," @ RoomT ",fecha_str.c_str()," (G=10^{A}V^{B})"));   //<-- change

  for(unsigned int i=0;i<_graph.size();i++)
  {
    int auxx3=i+1+tanda*5;

       //leg2->AddEntry(_graph[i],Form("%s %i %s %s  %s%.2lf %s %.2lf","PMT ",PMT_pos[i]," - ",serialnum[i].c_str(), " #chi^{2}/ndf = ",chi2[i],"/",ndf[i]));
    cout << "Fit successful?? " << _fit_successful[i] << endl;
     if (_fit_successful[i])
     {
       leg->AddEntry(_graph[i],Form("%s %i %s %s %s %.0lf\n %s %.0lf\n","PMT",auxx3,"-",Serialnum(round,i+5*tanda+1).c_str()," - #chi^{2}/ndf = ",chi2[i],"/",ndf[i]));       //<-- change
       leg->AddEntry(_myfit[i],Form("%s %.2lf\n %s %.2lf\n %s %.2lf\n %s %.2lf\n","Fit : A=",A[i],"#pm",errorA[i],"; B=",B[i],"#pm",errorB[i]),"L");
       leg->AddEntry(_graph[i],Form("%s %.1lf\n %s %.1lf\n %s %.1lf\n %s%.1lf\n %s ","HV (Gain=10^{7}) = ",vv[i],"#pm",errorV[i]," V - HV (Gain=10^{9}) = ",vv9[i],"#pm",errorV9[i],"V"),"");
     }
  }

  //  leg->AddEntry(myfit,"Fit: G=10^{A}V^{B}: A= -25.5231 #pm 0.0046; B= 10.6243 #pm 0.0014","L");
    leg->SetTextSize(0.034);
    leg->SetBorderSize(0.0);
    leg->Draw();

  setlogscale(c3,mg1,1000,1900);
lets_pause();
  cout << "lets log scale! " << endl;
  c3->Print("results/Superfit_exp.png");



}*/

void FitTGraphs(std::vector<TGraphErrors *> _graph, std::vector<string> label, string picname)
{

  TCanvas *c3 = new TCanvas("c3", "c3", 1200, 800);
  TMultiGraph *mg1 = new TMultiGraph("mg1", "");

  std::vector<TF1 *> _myfit;
  std::vector<bool> _fit_successful;
  std::vector<Double_t> A, errorA, B, errorB, vv, vv9, chi2, ndf, covarianzaAB, fA, fB, fB9, errorV, errorV9, errorVtilde, errorVtilde9;

  Double_t gainn = 1.e7;
  Double_t gain9 = 1.e9;

  for (unsigned int i = 0; i < _graph.size(); i++)
  {
    cout << "Creo TGraph " << i << endl;

    _fit_successful.push_back(false);

    _graph[i]->SetMarkerStyle(20 + i);
    _graph[i]->SetMarkerColor(2 + i);
    if (i == 0)
      _graph[i]->SetDrawOption("AP");
    else
      _graph[i]->SetDrawOption("P");
    _graph[i]->Draw();
    _graph[i]->SetFillStyle(0);

    cout << "Creo myfit " << i << endl;
    Double_t rangemin, rangemax;
    Double_t *x = _graph[i]->GetX();
    rangemin = x[0];
    rangemax = x[0];
    for (int j = 0; j < _graph[i]->GetN(); j++)
    {
      if (rangemin > x[j])
        rangemin = x[j];
      if (rangemax < x[j])
        rangemax = x[j];
    }
    cout << "Range fit " << rangemin << " " << rangemax << endl;
    _myfit.push_back(new TF1(Form("%s %i", "myfit_", i), "pow(10,[0])*pow(x,[1])", 1200, rangemax)); //<-- RANGE!
    _myfit[i]->SetParName(0, "A");
    _myfit[i]->SetParName(1, "B");
    _myfit[i]->SetLineStyle(4);
    _myfit[i]->SetLineWidth(3);
    _myfit[i]->SetLineColor(i + 2);

    cout << "Ajusto! " << endl;

    TFitResultPtr r = _graph[i]->Fit(Form("%s %i", "myfit_", i), "MSRE");
    TString estado = gMinuit->fCstatu;
    // float covarianzaAB;
    cout << "resultado ajuste: " << int(r) << " " << estado << endl;
    // lets_pause();
    if (int(r) == 0 && estado == "SUCCESSFUL")
    {
      cout << _fit_successful[i] << " " << i << endl;
      _fit_successful[i] = true;
      cout << _fit_successful[i] << " " << i << endl
           << endl;
      // TMatrixDSym cov = r->GetCovarianceMatrix();
      // cov.Print();
      // cout << cov.Print();
      // cout << "numero columnas" << cov.GetNrows() << endl;//" " <<   cov.GetSub(1,2)<< endl;
      covarianzaAB.push_back(r->GetCovarianceMatrix()[0][1]);
      cout << " eerror ferist par " << covarianzaAB[i] << endl;
    }
    else
    {
      _fit_successful[i] = true;
      covarianzaAB.push_back(0.0);
      cout << "ERROR!!!!" << endl;
      lets_pause();
    }
    // lets_pause();

    gStyle->SetOptFit(0000);

    A.push_back(_myfit[i]->GetParameter(0));
    errorA.push_back(_myfit[i]->GetParError(0));
    B.push_back(_myfit[i]->GetParameter(1));
    errorB.push_back(_myfit[i]->GetParError(1));
    chi2.push_back(_myfit[i]->GetChisquare());
    ndf.push_back(_myfit[i]->GetNDF());

    fA.push_back(-1 / B[i]);
    fB.push_back(-(TMath::Log10(gainn) - A[i]) / (B[i] * B[i]));
    fB9.push_back(-(TMath::Log10(gain9) - A[i]) / (B[i] * B[i]));

    errorVtilde.push_back(TMath::Sqrt(fA[i] * fA[i] * errorA[i] * errorA[i] + fB[i] * fB[i] * errorB[i] * errorB[i] + 2 * fA[i] * fB[i] * covarianzaAB[i]));
    errorVtilde9.push_back(TMath::Sqrt(fA[i] * fA[i] * errorA[i] * errorA[i] + fB9[i] * fB9[i] * errorB[i] * errorB[i] + 2 * fA[i] * fB9[i] * covarianzaAB[i]));

    // TMatrixDSym cov = _myfit[i]->GetCovarianceMatrix();
    //->GetChisquare();
    //--------------------------------------------------
    //              HV for nominal gain
    //--------------------------------------------------
    // G=10^A*V^B
    vv.push_back(pow(gainn / (pow(10, A[i])), 1 / B[i]));
    vv9.push_back(pow(gain9 / (pow(10, A[i])), 1 / B[i]));
    errorV.push_back(errorVtilde[i] * vv[i] / TMath::Log10(TMath::E()));
    errorV9.push_back(errorVtilde9[i] * vv9[i] / TMath::Log10(TMath::E()));

    cout << " " << endl;
    cout << "HV(gain=" << gainn << ")"
         << "=" << vv[i] << " "
         << "V" << endl;

    mg1->Add(_graph[i]);
    cout << "Acabo fit " << endl;
    // ofile2 << "Pos,Date,Analysis,log_10A,err,14k,Err,V_G=1e7,V_G=1e9" << endl;
    cout << Form("%.6lf", A[i]) << "," << Form("%.6lf", errorA[i]) << "," << Form("%.6lf", B[i]) << "," << Form("%.6lf", errorB[i]) << "," << covarianzaAB[i] << "," << Form("%.0lf", vv[i]) << "," << Form("%.6lf", errorV[i]) << "," << Form("%.0lf", vv9[i]) << "," << Form("%.6lf", errorV9[i]) << "," << fA[i] << "," << fB[i] << "," << fB9[i] << "," << errorVtilde[i] << "," << errorVtilde9[i] << endl;
  }

  mg1->Draw("AP");

  float x1 = 0.02;
  float y1 = 0.5;
  float x2 = x1 + 0.4; // anchura
  float y2 = 0.97;     // altura

  TLegend *leg = new TLegend(x1, y1, x2, y2);
  // leg->SetHeader(Form("%s %i %s %i %s","PMTs from ",1+tanda*5," to ",5+tanda*5," @ Room Temperature (G=10^{A}V^{B})"));   //<-- change
  leg->SetHeader(Form("%s (G=10^{A}V^{B})", picname.c_str())); //<-- change

  for (unsigned int i = 0; i < _graph.size(); i++)
  {
    // leg2->AddEntry(_graph[i],Form("%s %i %s %s  %s%.2lf %s %.2lf","PMT ",PMT_pos[i]," - ",serialnum[i].c_str(), " #chi^{2}/ndf = ",chi2[i],"/",ndf[i]));
    cout << "Fit successful?? " << _fit_successful[i] << endl;
    if (_fit_successful[i])
    {
      leg->AddEntry(_graph[i], Form("%s %s %.0lf\n %s %.0lf\n", label[i].c_str(), " - #chi^{2}/ndf = ", chi2[i], "/", ndf[i])); //<-- change
      leg->AddEntry(_myfit[i], Form("%s %.2lf\n %s %.2lf\n %s %.2lf\n %s %.2lf\n", "Fit : A=", A[i], "#pm", errorA[i], "; B=", B[i], "#pm", errorB[i]), "L");
      leg->AddEntry(_graph[i], Form("%s %.1lf\n %s %.1lf\n %s %.1lf\n %s%.1lf\n %s ", "HV (Gain=10^{7}) = ", vv[i], "#pm", errorV[i], " V - HV (Gain=10^{9}) = ", vv9[i], "#pm", errorV9[i], "V"), "");
    }
  }

  //  leg->AddEntry(myfit,"Fit: G=10^{A}V^{B}: A= -25.5231 #pm 0.0046; B= 10.6243 #pm 0.0014","L");
  leg->SetTextSize(0.034);
  leg->SetBorderSize(0.0);
  leg->Draw();

  setlogscale(c3, mg1, 1200, 1900);
  lets_pause();
  cout << "lets log scale! " << endl;
  c3->Print(Form("%s_superfit.png", picname.c_str()));
}

TH1D *FindTH1D(TFile *f1, string name)
{
  cout << "Holi" << endl;
  TIter next(f1->GetListOfKeys());
  TKey *key;
  while ((key = (TKey *)next()))
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1D"))
      continue;
    TH1D *h = (TH1D *)key->ReadObj();
    h->Print();
    string histname = h->GetName();
    if (histname.find(name) != std::string::npos)
    {
      cout << histname << " FOUND " << name << endl;
      return h;
    }
  }
  return NULL;
}
TH1F *FindTH1F(TFile *f1, string name)
{
  cout << "Holi" << endl;
  TIter next(f1->GetListOfKeys());
  TKey *key;
  while ((key = (TKey *)next()))
  {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1F"))
      continue;
    TH1F *h = (TH1F *)key->ReadObj();
    h->Print();
    string histname = h->GetName();
    if (histname.find(name) != std::string::npos)
    {
      cout << histname << " FOUND " << name << endl;
      return h;
    }
  }
  return NULL;
}
