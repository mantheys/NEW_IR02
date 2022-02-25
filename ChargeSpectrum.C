/*
Macro para obtener el espectro de carga de por canal y run.
Usa como input las ntuplas generadas con la macro SimpleAnalysis.C
*/

#include"lib/headers.h"

void Draw(std::vector<int> mych, std::vector<int> runs, std::vector<string> name, bool rate, std::vector<string> var, TVirtualPad *pad=NULL, float CustomMaxRange=0, double Rebin=2 )
{
  int n = runs.size(); 
  std::vector<double> Gains={7.90191,3.61,3.83}; //ganancias en pC
  std::vector<double> SPEAmp={38.6,24.8,25.5};//Amplitud del SPE en cuentas de ADC
  std::vector<variable> myvar(n);

  myvar = MyVar(myvar, n, var, mych, Gains, SPEAmp, Rebin);

  if(CustomMaxRange!=0) for(int i=0;i<n;i++) myvar[i].limitup=CustomMaxRange;

  TFile *ifile[n]; TNtuple *nt[n]; TH1D *h[n];
  for(int i=0;i<n;i++) ifile[i]= new TFile(Form("/pnfs/ciemat.es/data/neutrinos/Super-cells_LAr/Feb22/AnalysisROOT/run%02i_NTuple.root",runs[i])); 
  
  Float_t ch; Float_t time; Float_t _var; 
  
  for(int i=0;i<n;i++)
  {
    nt[i]= (TNtuple*)ifile[i]->Get("ntuple");
    nt[i]->SetBranchAddress("ch",&ch);
    nt[i]->SetBranchAddress("time",&time);
    nt[i]->SetBranchAddress(myvar[i].var.c_str(),&_var);
    if(rate) h[i]=new TH1D(Form("h%i",i),Form("%s;%s; Events/s",name[i].c_str(),myvar[i].title.c_str()),400,myvar[i].limitdown,myvar[i].limitup);
    else h[i]=new TH1D(Form("h%i",i),Form("%s;%s; NEvents",name[i].c_str(),myvar[i].title.c_str()),400,myvar[i].limitdown,myvar[i].limitup);
    h[i]->SetLineColor(i+1);if (i==4) h[i]->SetLineColor(i+3);//yellow sucks hard
    h[i]->SetLineWidth(2);
  }
  
  cout << "Drawing " << n<< endl;

  for(int i=0;i<n;i++)
  {
    for(int j=0;j<nt[i]->GetEntries();j++)
    { 
      nt[i]->GetEntry(j);
      if(ch!=mych[i]) continue;
      //if (runs[i]==19 && ch==1) {cout << ch << " " << _var << " " << time << " " << mych[i] << " - " << _var*myvar[i].UnitConversion << endl; lets_pause();}
      h[i]->Fill(_var*myvar[i].UnitConversion);
      //~cout<<_var*myvar[i].UnitConversion<<endl;//debug
    }

    nt[i]->GetEntry(0);
    double firstime=time;
    nt[i]->GetEntry(nt[i]->GetEntries()-1);
    double lasttime=time; 
    Double_t duration = (lasttime - firstime)*8.e-9;
    //~cout << duration << endl; lets_pause();
    if(rate) h[i]->Scale(1./duration);
    else ;//h[i]->Scale(1./h[i]->GetEntries());
  }

  TCanvas *c;
  if(!pad){ c  = new TCanvas("c"); c->cd();}
  else {pad->cd();}
  gStyle->SetOptTitle(0); gStyle->SetOptStat(0);
  gPad->SetLogy();
    
  for(int i=0;i<n;i++)
  { 
  
    //~gPad->SetGrid(1,1);
  
    h[i]->Draw("HIST SAME"); 
    h[i]->Rebin(4);
    //~h[i]->Smooth();

    gPad->Update();
    //~h[i]->GetXaxis()->SetRange(h[i]->FindBin(400),h[i]->FindBin(700));
    //~gPad->Update();
    //~cout<<h[i]->GetMean()<<endl; 
  }

  //~bool fit=false;
  bool fit=false;

  TF1 *f[n] ;  //gaussianas para el ajuste
  if (fit==true)
  {
    //Patricia code
    ofstream ofs("Gaus_Results.txt",std::ofstream::out | std::ofstream::app);
    ofs <<"----------------------------------------------"<<endl;

    for(int i=0;i<n;++i)
    {
      h[i]->GetYaxis()->SetRangeUser(1e-3,30);
      //~double min=380,max=680;
      double min=25,max=60;
      //~f[i] = new TF1("f","[p0]*exp(-0.5*((x-[p1])/[p2])*((x-[p1])/[p2]))+[p3]*exp(-0.5*((x-[p4])/[p5])*((x-[p4])/[p5]))  ",min,max);//PMT Doble Gaussiana
      f[i] = new TF1("f","gaus",min,max);
      float par[6]={1,450,10,1,550,10};
      for (int j=0;j<6;j++) f[i]->SetParameter(j,par[j]);
      h[i]->Fit(f[i],"R","SAME",min,max);
      f[0]->SetRange(myvar[i].limitdown,myvar[i].limitup);f[0]->Draw("SAME");f[i]->Print();gPad->Update();
      //~lets_pause(); 
      cout << runs[i] << "\t"<< mych[i]<<"\t" << f[i]->GetParameter(1) << "\t" << f[i]->GetParameter(2)<<"\t"<<var[i]<< endl;
      
    } 

    float totalrate[n];
    for (int j=0; j<n;j++)for (int i=-1000;i<1000;i++)totalrate[j]+=f[j]->Eval(i);

    cout<<"{";	
    //~for(int i=0;i<n-1;++i) cout<<f[i]->Integral(-1000,1000)/h[i]->GetBinWidth(1)<<",";
    //~cout<<f[n-1]->Integral(-1000,1000)/h[n-1]->GetBinWidth(1);
    for(int i=0;i<n-1;++i) cout<<totalrate[i]/h[i]->GetBinWidth(1)<<",";
    cout<<totalrate[n-1]/h[n-1]->GetBinWidth(1);
    //~cout<<f[0]->Integral(-1000,1000)/h[n-1]->GetBinWidth(1);
    cout<<"}"<<endl;	
  }

  // draw the legend
  TLegend *legend=new TLegend(0.7, 0.6, 0.9, 0.9);
  //~legend->SetTextFont(72);
  //~legend->SetTextSize(0.04);
  for(int i=0;i<n;++i)legend->AddEntry(h[i],name[i].c_str());
  if (fit==true)legend->AddEntry(f[0],"Gauss Fit","l");
  legend->Draw();
  gStyle->SetOptFit();
   
  TPaveText *t = new TPaveText(0.3, 0.85, 0.7, 0.95, "brNDC"); // middle
  t->AddText("SiPM Charge Spectrum - Triggers");
  t->Draw();
  lets_pause();
}

void ChargeSpectrum()
{
  Draw({0},{14},{"QTotal"},false,{"QTotal"});
}
