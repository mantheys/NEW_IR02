// ---------------------------------- ChargeSpectrum ---------------------------------- //

/* Macro para obtener el espectro de carga de por canal y run.
   Usa como input las ntuplas generadas con la macro SimpleAnalysis.C
   source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.18.04/x86_64-centos7-gcc48-opt/bin/thisroot.sh */

#include"lib/headers.h"

/* Draw_vector() -> Funcion con parametros vectorizados para todos los canales de manera que podemos representar todas las 
cargas/amplitudes de todos los dias en el mismo histograma. Ademas se puede hacer un fit a todos y guardar los valores en un txt :)
Se podria introducir una ganancia diferente para cada run aunque es más coherente usar un valor medio para todos los dias. */
void Draw_vector(string path, std::vector<int> mych, std::vector<int> runs, std::vector<string> name, bool rate, std::vector<string> var, std::vector<double> conv_fact,
double aef_sipm, double aef_pmt, double aef_sc, bool normalize=false, double Rebin=1, TVirtualPad *pad=NULL, float CustomMaxRange=0)
{
  std::vector<double> SPEAmp={38.6,24.8,25.5};//Amplitud del SPE en cuentas de ADC //de Rodrigo --> cambiar por nuevas
  int n = runs.size(); 
  std::vector<variable> myvar(n);
  myvar = MyVar(myvar, n, var, mych, conv_fact);
  if(CustomMaxRange!=0) for(int i=0;i<n;i++) myvar[i].limitup=CustomMaxRange;

  TFile *ifile[n];
  for(int i=0;i<n;i++) ifile[i]= new TFile(Form("%srun%02i_NTuple.root",path.c_str(),runs[i])); 

  TNtuple *nt[n];
  TNtuple *nt_ref[n];
  Float_t ch; 
  Float_t time; 
  Float_t _var; 

  TH1D *h[n];
  TH1D *aux[n];
  for(int i=0;i<n;i++)
  {
    nt_ref[i]= (TNtuple*)ifile[i]->Get("ntuple");
    nt[i]= (TNtuple*)ifile[i]->Get("charge");
    // lets_pause();
    // cout<<"hecho1"<<endl;
    nt_ref[i]->SetBranchAddress("ch",&ch);
    // lets_pause();
    // cout<<"hecho2"<<endl;
    nt_ref[i]->SetBranchAddress("time",&time);
    // lets_pause();
    // cout<<"hecho3"<<endl;
    nt[i]->SetBranchAddress(myvar[i].var.c_str(),&_var);
    // lets_pause();
    // cout<<"hecho4"<<endl;

    if(rate) h[i]=new TH1D(Form("h%i",i),Form("%s;%s; Events/s",name[i].c_str(),myvar[i].title.c_str()),400,myvar[i].limitdown,myvar[i].limitup);
    else h[i]=new TH1D(Form("h%i",i),Form("%s;%s; NEvents",name[i].c_str(),myvar[i].title.c_str()),400,myvar[i].limitdown,myvar[i].limitup);
    h[i]->SetLineColor(i+1);if (i==4) h[i]->SetLineColor(i+3); if (i==9) h[i]->SetLineColor(30); // avoid bad colors
    h[i]->SetLineWidth(2);
    if(rate) aux[i]=new TH1D(Form("h%i",i),Form("%s;%s; Events/s",name[i].c_str(),myvar[i].title.c_str()),400,myvar[i].limitdown,myvar[i].limitup);
    else aux[i]=new TH1D(Form("h%i",i),Form("%s;%s; NEvents",name[i].c_str(),myvar[i].title.c_str()),400,myvar[i].limitdown,myvar[i].limitup);
    aux[i]->SetLineColor(i+1);if (i==4) aux[i]->SetLineColor(i+3); if (i==9) aux[i]->SetLineColor(30); // avoid bad colors
    aux[i]->SetLineWidth(2);
  }

  cout << "Drawing " << n<< endl;
  double max_charge;
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<nt[i]->GetEntries();j++)
    { 
      nt[i]->GetEntry(j);
      nt_ref[i]->GetEntry(j);
      if(ch!=mych[i]) continue;
      //if (runs[i]==19 && ch==1) {cout << ch << " " << _var << " " << time << " " << mych[i] << " - " << _var*myvar[i].UnitConversion << endl; lets_pause();}
      h[i]->Fill(_var*myvar[i].UnitConversion);
      //cout<<_var*myvar[i].UnitConversion<<endl;//debug
    }
    aux[i]= (TH1D*)h[i]->Clone();
    // cout << h[i]->GetMaximum();
    h[i]->Rebin(Rebin);
    if(normalize)
    {
      h[i]->Scale(1./ h[i]->GetMaximum());
    }
  
    nt_ref[i]->GetEntry(0);
    double firstime=time;
    nt_ref[i]->GetEntry(nt_ref[i]->GetEntries()-1);
    double lasttime=time; 
    Double_t duration = (lasttime - firstime)*8.e-9;
    //cout << duration << endl; lets_pause();
    if(rate) h[i]->Scale(1./duration);
    // else ;//h[i]->Scale(1./h[i]->GetEntries());
  }

  TCanvas *c;
  if(!pad){ c  = new TCanvas("c"); c->cd();}
  else {pad->cd();}
  gStyle->SetOptTitle(0); gStyle->SetOptStat(0);
  // gPad->SetLogy();

  for(int i=0;i<n;i++)
  { 
    //gPad->SetGrid(1,1);
    h[i]->Draw("HIST SAME"); 
    // h[i]->Rebin(4);
    //h[i]->Smooth();

    gPad->Update();
    //h[i]->GetXaxis()->SetRange(h[i]->FindBin(400),h[i]->FindBin(700));
    //gPad->Update();
    //cout<<h[i]->GetMean()<<endl; 
  }

  bool fit=true;
  // bool fit=false;
  // draw the legend
  TLegend *legend=new TLegend(0.7, 0.6, 0.9, 0.9);
  //legend->SetTextFont(72);
  //legend->SetTextSize(0.04);
  if(fit!=true) for(int i=0;i<n;++i) legend->AddEntry(h[i],name[i].c_str());
  
  if(fit==true)
  {
    lets_pause();
    double min_fit, max_fit;
	  cout << "Min value for fit: "; cin >> min_fit;
	  cout << "Max value for fit: "; cin >> max_fit;
    string filename;
    cout << "File name?:"; cin >> filename;
    
    for(int i=0;i<n;i++)
    {
      int npeaks = 1;
	    double_t par[100];
	    h[i]->GetXaxis()->SetRangeUser(min_fit, max_fit);
	    TSpectrum *s = new TSpectrum(2*npeaks,1);
	    Int_t nfound = s->Search(h[i],2,"goff",0.01);
	    printf("Found %d candidate peaks to fit\n",nfound);
	    h[i]->GetXaxis()->SetRange(0,0);
	    Double_t *xpeaks;
	    xpeaks = s->GetPositionX();
	
   	  for (int p=0;p<npeaks;p++)
      {
        Double_t xp = xpeaks[p];
        Int_t bin = h[i]->GetXaxis()->FindBin(xp);
        Double_t yp = h[i]->GetBinContent(bin);
        //if (yp-TMath::Sqrt(yp) < fline->Eval(xp)) continue;
        par[3*p] = yp; // "height"
        par[3*p+1] = xp; // "mean"
        par[3*p+2] = 1; // "sigma"
        // cout << Form("sigma_%i: ", p); cin >> par[3*npeaks+2];
        //~ p++;
	    }
	   
     //int_t max_bin[npeaks], j=0;
	   string init = Form("gaus(0)");
	   string added;
	   for (int i=1; i<npeaks; i++){
   		added = Form("+gaus(%i)",i*3);
	   	init = init + added;
		}

    //int Nbins = h->GetNbinsX();
    TF1 *fb = new TF1("fb",init.c_str(),min_fit,max_fit);
    fb->SetParameters(par);
    fb->SetLineColor(i+1); if (i==4) fb->SetLineColor(i+3); if (i==9) fb->SetLineColor(30);
    fb->Draw("SAME");
    legend->AddEntry(aux[i],name[i].c_str(),"l");
    h[i]->Fit("fb","MER+","",min_fit, max_fit); gPad->Update();
    //  if (i>0 && i<n) h[i]->Delete();
    h[i]->SetLineWidth(0);
    cout << h[i]->GetBinWidth(h[i]->GetMaximumBin()) << endl;
    // TH1 *hd = (TH1*)gPad->GetPrimitive(Form("h%i",i)); gPad->GetListOfPrimitives()->Remove(hd);
    //  lets_pause();
    { // ------------ Volcamos a txt --------------- 
    double mu = fb->GetParameter(3*0+1), Dmu = fb->GetParError(3*0+1);
	  double sigma = fb->GetParameter(3*0+2), Dsigma = fb->GetParError(3*0+2);
	  double amp = fb->GetParameter(0), Damp = fb->GetParError(0);
    
	  // FILE *f = fopen(Form("prueba.txt"), "a");
	  FILE *f = fopen(Form("fit_data/%s.txt",filename.c_str()), "a");
	  if (f == NULL)
	  {
		printf("Error opening file!\n");
		exit(1);
	  }
	  fprintf(f, "%f\t%f\t%f\t%f\t%f\t%f\n", mu, Dmu, sigma, Dsigma, amp, Damp);
	  fclose(f);
	  } //*/
	  //  lets_pause();
    }
  }

  legend->Draw();
  //  gStyle->SetOptFit();
   
TPaveText *t = new TPaveText(0.3, 0.85, 0.7, 0.95, "brNDC"); // middle
t->AddText("SiPM #290 - OV = 4.0 V");
// t->AddText("SuperCell Charge Spectrum - OV = 7.0 V");
// t->AddText("SuperCell PE - OV = 3.5");
// t->AddText("SuperCell PE/mm^2 - OV = 3.5");
// t->AddText("PMT Charge (1400V)");
t->Draw();
lets_pause();
lets_pause();
}

void ChargeSpectrum(string input = "cs_config_file.txt")
{
  int ch;
  ch = IntInput(input, "CHANNEL");

  string path; string range_type;
  path = StringInput(input, "PATH"); range_type = StringInput(input, "RANGE_TYPE");

  bool debug;
  debug = BoolInput(input, "DEBUG");

  /* >> Empleamos el json para seleccionar las runes */
  const Run_map map_of_runs = Run_map("mapa_feb_2.json");
  std::vector<int> runs;
  std::vector<string> label;
  std::vector<double> gains_sipm289;
  std::vector<double> gains_sipm290;
  std::vector<double> gains_pmt;
  std::vector<double> gains_sc;
  std::vector<double> conv_fact;

  for (auto i:map_of_runs.json_map) 
  {
    if (i["tipo"] == "coincidencia" && i["threshold"] == 200 &&  i["ov_sipm"] == 4 &&  i["ov_pmt"] == 1200)
    {
      runs.push_back(i["run"]);
      label.push_back(to_string(i["dia"])+" Feb");
    }
    if (i["run"] == 1)
    {
      conv_fact.push_back(i["factor_conv_sipm"]);
      conv_fact.push_back(i["factor_conv_pmt"]);
      conv_fact.push_back(i["factor_conv_sc"]);
    }
  }

  int n = runs.size();
  std::vector<int> channels(n,ch);

  vector<string> var(n,range_type);

  /* >> Ganancias si no estan leidas del json */
  std::vector<double> gains_sipm1(n,6.77E+06); //Mean + 4OV + ped-1PE


  if (debug)
  {
    for (auto w:channels) cout << "Channels: " << w << endl;
    for (auto test:runs) cout << "Runs: " << test << endl;
    for (auto l:label) cout << l << endl;
    for (auto i:var) cout << "Charge Range: " << i << endl;
    for (auto test5:conv_fact) cout << test5 << endl;
  }

  /* En esta funcion las ganancias de todos los canales estan definidos como vectores */
  Draw_vector(path, channels, runs, label, true, var, conv_fact, 236, 415.47, 415.47, false, 1);
}