#ifndef _THREEEXPOGAUSFITSCINTUTILSCH_
#define _THREEEXPOGAUSFITSCINTUTILSCH_

#include <TFitResultPtr.h>
#include <TF1.h>
#include <TProfile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TStyle.h>
#include <TPad.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TMatrixDSym.h>

using namespace std;

double func(double *x, double *par, int index = 3);
double scint_time_1(double *x, double *par);
double scint_time_2(double *x, double *par);
double scint_time_3(double *x, double *par);

class ThreeExpoGausFitScintUtilsCh : public TObject
{

public:
  ThreeExpoGausFitScintUtilsCh();
  ~ThreeExpoGausFitScintUtilsCh(){};

  void Init();
  void Reset();
  void SetLimits();

  void Fit(TH1 &h);

public:
  int _fit_n_comp;

  TF1 *_func;

  TFitResultPtr _fit_res;

  double _fit_x_min;
  double _fit_x_max;

  double _init_ped;
  double _init_t0;
  double _init_sigma;
  double _init_A_fast;
  double _init_t_fast;
  double _init_A_int;
  double _init_t_int;
  double _init_A_slow;
  double _init_t_slow;

  double _ped;
  double _t0;
  double _sigma;
  double _A_fast;
  double _t_fast;
  double _A_int;
  double _t_int;
  double _A_slow;
  double _t_slow;
  double slowrange = 4e-6;
  double sampling;
  void SetMaxRange(double m) { slowrange = m; }
};
double ExcludeMin, ExcludeMax;
Bool_t reject;

inline ThreeExpoGausFitScintUtilsCh::ThreeExpoGausFitScintUtilsCh()
{

  _fit_n_comp = 3;

  _func = 0;
  _fit_res = 0;

  _fit_x_min = 0.3e-6;
  _fit_x_max = 5.e-6;

  _init_t0 = -1;
  _init_sigma = -1;
  _init_A_fast = -1;
  _init_t_fast = -1;
  _init_A_int = -1;
  _init_t_int = -1;
  _init_A_slow = -1;
  _init_t_slow = -1;

  _t0 = -1;
  _sigma = -1;
  _A_fast = -1;
  _t_fast = -1;
  _A_int = -1;
  _t_int = -1;
  _A_slow = -1;
  _t_slow = -1;
}
inline void ThreeExpoGausFitScintUtilsCh::SetLimits()
{

  double factor = sampling / 16.e-9;

  _func->SetParLimits(0, 0.1e-6, 1.e-5);
  _func->SetParLimits(1, _init_t0 - 100e-9, _init_t0 + 100e-9);
  _func->SetParLimits(2, factor * 3e-9, 20e-9);
  _func->SetParLimits(3, factor * 0.1e-9, factor * 20e-9); // A_Fast
  _func->SetParLimits(5, factor * 10e-9, factor * 60e-9);   // A_int
  _func->SetParLimits(6, 10e-9, 50e-9);
  _func->SetParLimits(7, factor * 0.9e-9, 2e-9); // A_slow
  _func->SetParLimits(8, 0.1e-6, .9e-6);
}
inline void ThreeExpoGausFitScintUtilsCh::Init()
{

  if (_func)
  {
    delete _func;
  }

  if (_fit_n_comp == 3)
  {

    _func = new TF1("scint_time", scint_time_3, _fit_x_min, _fit_x_max, 9);

    _func->SetParName(0, "Ped");
    _func->SetParName(1, "t_{0}");
    _func->SetParName(2, "#sigma");
    _func->SetParName(3, "A_{fast}");
    _func->SetParName(4, "#tau_{fast}");
    _func->SetParName(5, "A_{int}");
    _func->SetParName(6, "#tau_{int}");
    _func->SetParName(7, "A_{slow}");
    _func->SetParName(8, "#tau_{slow}");

    _func->SetParameter(0, _init_ped);
    _func->SetParameter(1, _init_t0);
    _func->SetParameter(2, _init_sigma);
    _func->SetParameter(3, _init_A_fast);
    _func->SetParameter(4, _init_t_fast);
    _func->SetParameter(5, _init_A_int);
    _func->SetParameter(6, _init_t_int);
    _func->SetParameter(7, _init_A_slow);
    _func->SetParameter(8, _init_t_slow);

    // my setting parameters
    _func->FixParameter(4, _init_t_fast);
    SetLimits();
  }
  else
  {

    Fatal("ThreeExpoGausFitScintUtilsCh::Init", "_fit_n_comp = %i not implememnted", _fit_n_comp);
  }

  _func->SetNpx(1000);
}

inline void ThreeExpoGausFitScintUtilsCh::Reset()
{

  delete _func;
}

inline void ThreeExpoGausFitScintUtilsCh::Fit(TH1 &h)
{

  TString status_fit;

  // Pedestal subtraction

  sampling = h.GetBinWidth(1);
  _fit_x_min = h.GetXaxis()->GetBinCenter(h.GetMaximumBin()) - 500e-9;
  _fit_x_max = h.GetXaxis()->GetBinCenter(h.GetMaximumBin()) - 100e-9;
  TF1 *_func_pedestal = new TF1("pedestal", "pol0", _fit_x_min, _fit_x_max);
  _func_pedestal->SetParameter(0, h.GetBinContent(h.GetXaxis()->FindBin(h.GetXaxis()->GetBinCenter(h.GetMaximumBin()) - 250e-9)));
  TFitResultPtr r = h.Fit(_func_pedestal, "LNQ", "", _fit_x_min, _fit_x_max);
  _init_ped = _func_pedestal->GetParameter(0);
  std::cout << "Pedestal fit range is " << _fit_x_min << " " << _fit_x_max << std::endl;
  std::cout << "Pedestal ini is " << h.GetBinContent(3) << std::endl;
  std::cout << "Pedestal is " << _init_ped << std::endl;
  // if(isnan(_init_ped)) throw std::exception();
  int i_max = h.GetMaximumBin();
  double x_max = h.GetBinCenter(i_max);
  double max = h.GetMaximum();
  if (_init_ped < 0)
    _init_ped = -_init_ped;

  // T0 estimation
  _init_t0 = x_max;
  std::cout << "_init_t0 " << _init_t0 << std::endl;

  // Amplitude estimation
  int i_l_max = i_max;
  while (h.GetBinContent(i_l_max) > 10 * _init_ped)
  {
    i_l_max--; // std::cout << "i_l_max++ " << i_l_max++ << std::endl;
  }

  int i_r_max = i_max;
  while (h.GetBinContent(i_r_max) > max / 2)
  {
    i_r_max++; // std::cout << "i_r_max++ " << i_r_max++ << std::endl;
  }

  TF1 *_func_fast = new TF1("around_fast", "gaus", h.GetBinCenter(i_l_max), h.GetBinCenter(i_r_max));
  _func_fast->SetParameter(1, x_max);
  _func_fast->SetParameter(2, 1e-8);
  _func_fast->SetParLimits(1, x_max - 10e-9, x_max + 10e-9);
  _func_fast->SetParLimits(2, 0.2e-08, 7e-08);
  h.Fit(_func_fast, "LNQ", "", h.GetBinCenter(i_r_max) - 0.1e-6, h.GetBinCenter(i_r_max));
  _init_A_fast = 0.1 * _func_fast->Integral(_func_fast->GetParameter(1) - 5 * _func_fast->GetParameter(2), _func_fast->GetParameter(1) + 5 * _func_fast->GetParameter(2));
  _init_sigma = _func_fast->GetParameter(2);

  std::cout << "func limits is " << h.GetBinCenter(i_l_max) << " " << h.GetBinCenter(i_r_max) << std::endl;
  std::cout << "_init_sigma is " << _init_sigma << std::endl;

  // Int component estimation
  double min_expo = _init_t0 + 15e-9; // 7*_init_sigma;
  double max_expo = _init_t0 + 40e-9; // 12*_init_sigma;
  TF1 *_func_int = new TF1("int", "expo", min_expo, max_expo);
  _func_int->SetParameter(1, -1 / 50e-9);

  h.Fit(_func_int, "LNQ", "", min_expo, max_expo);

  std::cout << "min_expo, max_expo is " << min_expo << " " << max_expo << std::endl;
  _init_A_int = 6 * _init_A_fast;
  _init_t_int = -1 / _func_int->GetParameter(1);
  std::cout << "_init_t_int is " << _init_t_int << std::endl;

  ExcludeMin = x_max + 30e-9;
  ExcludeMax = x_max + 90e-9;
  // Slow component estimation
  TF1 *_func_slow = new TF1("slow", "expo", x_max + 0.6e-6, x_max + slowrange); // 3300);//500,1000 ,4000);

  h.Fit(_func_slow, "LNQ", "", x_max + 0.6e-6, x_max + slowrange); // lets_pause();lets_pause();
  h.Draw("HIST");
  _func_slow->Draw("same");
  _func_pedestal->Draw("same");
  _func_fast->Draw("same");
  _func_int->Draw("same");

  _init_t_fast = 6e-9;

  _init_A_slow = 0.3 * _init_A_fast;
  _init_t_slow = -1 / _func_slow->GetParameter(1);
  std::cout << "_init_t_slow is " << _init_t_slow << std::endl;

  cout << "           ---------   " << endl;
  cout << "Init Ped      : " << _init_ped << endl;
  cout << "Init t0       : " << _init_t0 << endl;
  cout << "Init sigma    : " << _init_sigma << endl;
  cout << "Init A_{Fast} : " << _init_A_fast << endl;
  cout << "Init t_{Fast} : " << _init_t_fast << endl;
  cout << "Init A_{Int}  : " << _init_A_int << endl;
  cout << "Init t_{Int}  : " << _init_t_int << endl;
  cout << "Init A_{Slow} : " << _init_A_slow << endl;
  cout << "Init t_{Slow} : " << _init_t_slow << endl;

  Init();
  cout << reject << " reject value " << endl;
  gStyle->SetOptFit(1112);
  cout << "Fit Range: " << x_max - 0.5e-6 << ", " << x_max + slowrange << endl;
  // TF1 *f = new TF1("scint_time", scint_time_3, _fit_x_min, _fit_x_max, 9);

  _func->SetRange(x_max - 0.5e-6, x_max + slowrange);             // h.Sumw2(1);
  _fit_res = h.Fit(_func, "LQESN", "", x_max - 0.5e-6, x_max + slowrange); // h.Draw("HIST"); lets_pause();
  _fit_res = h.Fit(_func, "LQE", "", x_max - 0.5e-6, x_max + slowrange);   // lets_pause(); lets_pause();
  _fit_res = h.Fit(_func, "LQESN", "", x_max - 0.5e-6, x_max + slowrange);
  _fit_res = h.Fit(_func, "LMESN", "", x_max - 0.5e-6, x_max + slowrange);
  _func->FixParameter(0, _func->GetParameter(0));
  _func->FixParameter(1, _func->GetParameter(1));
  _func->FixParameter(2, _func->GetParameter(2));
  _func->FixParameter(3, _func->GetParameter(3));
  _func->FixParameter(4, _func->GetParameter(4));
  _func->FixParameter(5, _func->GetParameter(5));
  _func->FixParameter(6, _func->GetParameter(6));
  //        _func->SetRange(x_max+1.e-6, x_max+4.e-6); //h.Sumw2(1);
  _fit_res = h.Fit(_func, "LMESN", "", x_max + 1.e-6, x_max + slowrange);
  _func->SetRange(x_max - 0.5e-6, x_max + slowrange); // h.Sumw2(1);

  _func->ReleaseParameter(0);
  _func->ReleaseParameter(1);
  _func->ReleaseParameter(2);
  _func->ReleaseParameter(3);
  _func->ReleaseParameter(5);
  _func->ReleaseParameter(6);
  SetLimits();
  _func->FixParameter(7, _func->GetParameter(7));
  _func->FixParameter(8, _func->GetParameter(8));
  _fit_res = h.Fit(_func, "LMESN", "", x_max - 0.5e-6, x_max + 0.2e-6);

  //        _func->ReleaseParameter(7);
  //        _func->ReleaseParameter(8);
  SetLimits();
  //        _func->FixParameter(0,_func->GetParameter(0));
  //        _func->FixParameter(1,_func->GetParameter(1));
  //        _func->FixParameter(2,_func->GetParameter(2));
  //        _func->FixParameter(3,_func->GetParameter(3));
  //        _func->FixParameter(5,_func->GetParameter(5));
  //        _func->FixParameter(6,_func->GetParameter(6));

  _fit_res = h.Fit(_func, "LMESN", "", x_max + 1.5e-6, x_max + slowrange);
  cout << reject << " reject value changed " << endl; // lets_pause();
  SetLimits();
  reject = kTRUE;
  cout << reject << " reject value changed " << endl; // lets_pause();
  _fit_res = h.Fit(_func, "LMSN", "", x_max - .5e-6, x_max + slowrange);
  _fit_res = h.Fit(_func, "LMSN", "", x_max - .5e-6, x_max + slowrange);
  _fit_res = h.Fit(_func, "LM", "", x_max - .5e-6, x_max + slowrange);

  reject = false; // delete rejected range to draw the full function

  h.Draw("HIST");
  h.GetYaxis()->SetRangeUser(1e-5, 2);
  _func->Draw("SAME"); // lets_pause();lets_pause();lets_pause();lets_pause();lets_pause();
  h.GetListOfFunctions()->Add(_func);
  //  ((TF1*)h.GetListOfFunctions()->FindObject("scint_time"))->SetRange(x_max-0.5e-6, x_max+4.e-6);
  cout << "           ---------   " << endl;
  cout << "     ----> I'm here!   " << endl;
  cout << "           ---------   " << endl;
}

// Fit functions
double func(double *x, double *par, int index)
{

  double P = par[0] / 3.0; // Pedestal

  double TT = par[1];   // PMT Transit Time
  double TTS = par[2]; // PMT Transit Time Spread (Jitter)

  double A = par[index];     // Normalisation
  double tau = par[index + 1]; // Time constant

  double t = x[0]; // Variable

  double c = A * 4.; //*1000;// if scale the histo
  c /= (2. * tau);

  double t1 = (TTS * TTS) / (2 * tau * tau);

  double t2 = t - TT;
  t2 /= tau;

  double t3 = (TTS * TTS) - (tau * (t - TT));
  t3 /= (TMath::Sqrt(2) * TTS * tau);

  return P + c * TMath::Exp(t1 - t2) * (1 - TMath::Erf(t3));
}

double scint_time_1(double *x, double *par)
{
  return func(x, par, 3);
}

double scint_time_2(double *x, double *par)
{
  return func(x, par, 3) + func(x, par, 5);
}

double scint_time_3(double *x, double *par)
{
  if (reject && x[0] > ExcludeMin && x[0] < ExcludeMax)
  {
    TF1::RejectPoint(); // cout << "point excluded " << x[0] << endl; lets_pause();
    return 0;
  }
  return func(x, par, 3) + func(x, par, 5) + func(x, par, 7);
}

#endif
