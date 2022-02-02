/*
Macro para obtener el espectro de carga de por canal y run.
Usa como input las ntuplas generadas con la macro SimpleAnalysis.C
*/
#include"lib/headers.h"

struct variable
{
  string var;
  string title;
  double limitdown;
  double limitup;
  double UnitConversion;
};
void Draw(std::vector<int> mych, std::vector<int> runs, std::vector<string> name, bool rate, std::vector<string> var, TVirtualPad *pad=NULL, float CustomMaxRange=0, double Rebin=2 )
{
  //~std::vector<double> Gains={0.59,3.61,3.83}; //ganancias en pC
  std::vector<double> Gains={7.90191,3.61,3.83}; //ganancias en pC
  
  std::vector<double> SPEAmp={38.6,24.8,25.5};//Amplitud del SPE en cuentas de ADC
 int n=runs.size(); 
 std::vector<variable> myvar(n);
 for(int i=0;i<n;i++)
 {
  if (var[i]=="TStartQpeak" ) myvar[i]={"TStartQpeak","T (s)",1e-6,3e-6,1.0};
  if (var[i]=="TEndQPeak" ) myvar[i]={"TEndQPeak","Time (s)",1e-6,3e-6,1.0};
  if (var[i]=="TDiff" ) myvar[i]={"TDiff","Time (s)",0e-6,2e-6,1.0};
  
  if (var[i]=="Q1" && mych[i]==0) myvar[i]={"Q1","Charge collected (pC)",0,3000,1.0};
  if (var[i]=="Q1" && mych[i]!=0) myvar[i]={"Q1","Charge collected (pC)",0,500,1.0};
  
  if (var[i]=="Q2" && mych[i]==0) myvar[i]={"Q2","Charge collected (pC)",0,3000,1.0};
  if (var[i]=="Q2" && mych[i]!=0) myvar[i]={"Q2","Charge collected (pC)",0,500,1.0};
  
  if (var[i]=="Q3" && mych[i]==0) myvar[i]={"Q3","Charge collected (pC)",0,3000,1.0};
  if (var[i]=="Q3" && mych[i]!=0) myvar[i]={"Q3","Charge collected (pC)",0,500,1.0};
  
  if (var[i]=="QTotal" && mych[i]==0) myvar[i]={"QTotal","Charge collected (pC)",0,3000,1.0};
  if (var[i]=="QTotal" && mych[i]!=0) myvar[i]={"QTotal","Charge collected (pC)",0,500,1.0};
  
  
  if (var[i]=="Amp" && mych[i]==0) myvar[i]={"Amp","Amplitude (ADC)",0,2000,1.0};
  if (var[i]=="Amp" && mych[i]!=0) myvar[i]={"Amp","Amplitude (ADC)",0,2000,1.0};

  if (var[i]=="Q1PE"&& mych[i]!=0) myvar[i]={"Q1","Charge collected (PE)",0,130,1.0/Gains[mych[i]]};
  if (var[i]=="Q1PE"&& mych[i]==0) myvar[i]={"Q1","Charge collected (PE)",0,1000,1.0/Gains[mych[i]]};
  
  if (var[i]=="Q2PE"&& mych[i]!=0) myvar[i]={"Q2","Charge collected (PE)",0,130,1.0/Gains[mych[i]]};
  if (var[i]=="Q2PE"&& mych[i]==0) myvar[i]={"Q2","Charge collected (PE)",0,1000,1.0/Gains[mych[i]]};
  
  if (var[i]=="Q3PE"&& mych[i]!=0) myvar[i]={"Q3","Charge collected (PE)",0,130,1.0/Gains[mych[i]]};
  if (var[i]=="Q3PE"&& mych[i]==0) myvar[i]={"Q3","Charge collected (PE)",0,1000,1.0/Gains[mych[i]]};

  if (var[i]=="QPeak"&& mych[i]==0) myvar[i]={"QPeak","Charge Peak",0,2000,1.0};
  if (var[i]=="QPeak"&& mych[i]!=0) myvar[i]={"QPeak","Charge Peak",0,500,1.0};

  if (var[i]=="QPeakPE"&& mych[i]==0) myvar[i]={"QPeakPE","Charge Peak (PE)",0,6000*Rebin,1.0/Gains[mych[i]]};
  if (var[i]=="QPeakPE"&& mych[i]!=0) myvar[i]={"QPeakPE","Charge Peak (PE)",0,110*Rebin,1.0/Gains[mych[i]]};

  if (var[i]=="AmpPE"&& mych[i]!=0) myvar[i]={"Amp","Amplitude (PE)",0,200,1.0/SPEAmp[mych[i]]};
  if (var[i]=="AmpPE"&& mych[i]==0) myvar[i]={"Amp","Amplitude (PE)",0,1000,1.0/SPEAmp[mych[i]]};
 }
 if(CustomMaxRange!=0) for(int i=0;i<n;i++) myvar[i].limitup=CustomMaxRange;

 TFile *ifile[n];
 for(int i=0;i<n;i++) ifile[i]= new TFile(Form("AnalysisROOT/Run%i_NTuple.root",runs[i])); 
 TNtuple *nt[n];
 Float_t ch; 
 Float_t time; 
 Float_t _var; 
  TH1D *h[n];
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
 for(int i=0;i<n;i++){ 
	 
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
if (fit==true){
 //Patricia code
 ofstream ofs("Gaus_Results.txt",std::ofstream::out | std::ofstream::app);
 ofs <<"----------------------------------------------"<<endl;

	for(int i=0;i<n;++i){
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

void Draw_April(std::vector<int> mych, std::vector<int> runs, std::vector<string> name, bool rate, std::vector<string> var, TVirtualPad *pad=NULL, float CustomMaxRange=0, double Rebin=2 )
{
  std::vector<double> Gains={0.59,3.61,3.83}; //ganancias en pC Abril
  std::vector<double> SPEAmp={38.6,24.8,25.5};//Amplitud del SPE en cuentas de ADC
 int n=runs.size(); 
 std::vector<variable> myvar(n);
 for(int i=0;i<n;i++)
 {
  if (var[i]=="TStartQpeak" ) myvar[i]={"TStartQpeak","T (s)",1e-6,3e-6,1.0};
  if (var[i]=="TEndQPeak" ) myvar[i]={"TEndQPeak","Time (s)",1e-6,3e-6,1.0};
  if (var[i]=="TDiff" ) myvar[i]={"TDiff","Time (s)",0e-6,2e-6,1.0};
  
  if (var[i]=="Q1" && mych[i]==0) myvar[i]={"Q1","Charge collected (pC)",0,1000,1.0};
  if (var[i]=="Q1" && mych[i]!=0) myvar[i]={"Q1","Charge collected (pC)",0,500,1.0};
  
  if (var[i]=="Q2" && mych[i]==0) myvar[i]={"Q2","Charge collected (pC)",0,1000,1.0};
  if (var[i]=="Q2" && mych[i]!=0) myvar[i]={"Q2","Charge collected (pC)",0,500,1.0};
  
  if (var[i]=="Q3" && mych[i]==0) myvar[i]={"Q3","Charge collected (pC)",0,1000,1.0};
  if (var[i]=="Q3" && mych[i]!=0) myvar[i]={"Q3","Charge collected (pC)",0,500,1.0};
  
  if (var[i]=="QTotal" && mych[i]==0) myvar[i]={"QTotal","Charge collected (pC)",0,600,1.0};
  if (var[i]=="QTotal" && mych[i]!=0) myvar[i]={"QTotal","Charge collected (pC)",0,500,1.0};
  
  
  if (var[i]=="Amp" && mych[i]==0) myvar[i]={"Amp","Amplitude (ADC)",0,2000,1.0};
  if (var[i]=="Amp" && mych[i]!=0) myvar[i]={"Amp","Amplitude (ADC)",0,2000,1.0};

  if (var[i]=="Q1PE"&& mych[i]!=0) myvar[i]={"Q1","Charge collected (PE)",0,130,1.0/Gains[mych[i]]};
  if (var[i]=="Q1PE"&& mych[i]==0) myvar[i]={"Q1","Charge collected (PE)",0,1000,1.0/Gains[mych[i]]};
  
  if (var[i]=="Q2PE"&& mych[i]!=0) myvar[i]={"Q2","Charge collected (PE)",0,130,1.0/Gains[mych[i]]};
  if (var[i]=="Q2PE"&& mych[i]==0) myvar[i]={"Q2","Charge collected (PE)",0,1000,1.0/Gains[mych[i]]};
  
  if (var[i]=="Q3PE"&& mych[i]!=0) myvar[i]={"Q3","Charge collected (PE)",0,130,1.0/Gains[mych[i]]};
  if (var[i]=="Q3PE"&& mych[i]==0) myvar[i]={"Q3","Charge collected (PE)",0,1000,1.0/Gains[mych[i]]};

  if (var[i]=="QPeak"&& mych[i]==0) myvar[i]={"QPeak","Charge Peak",0,2000,1.0};
  if (var[i]=="QPeak"&& mych[i]!=0) myvar[i]={"QPeak","Charge Peak",0,500,1.0};

  if (var[i]=="QPeakPE"&& mych[i]==0) myvar[i]={"QPeakPE","Charge Peak (PE)",0,6000*Rebin,1.0/Gains[mych[i]]};
  if (var[i]=="QPeakPE"&& mych[i]!=0) myvar[i]={"QPeakPE","Charge Peak (PE)",0,110*Rebin,1.0/Gains[mych[i]]};

  if (var[i]=="AmpPE"&& mych[i]!=0) myvar[i]={"Amp","Amplitude (PE)",0,200,1.0/SPEAmp[mych[i]]};
  if (var[i]=="AmpPE"&& mych[i]==0) myvar[i]={"Amp","Amplitude (PE)",0,1000,1.0/SPEAmp[mych[i]]};
 }
 if(CustomMaxRange!=0) for(int i=0;i<n;i++) myvar[i].limitup=CustomMaxRange;

 TFile *ifile[n];
 for(int i=0;i<n;i++) ifile[i]= new TFile(Form("AnalysisROOT_April/Run%i_NTuple.root",runs[i])); 
 TNtuple *nt[n];
 Float_t ch; 
 Float_t time; 
 Float_t _var; 
  TH1D *h[n];
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
 for(int i=0;i<n;i++){ 
	 
	 //~gPad->SetGrid(1,1);
	 
	 h[i]->Draw("HIST SAME"); 
	 h[i]->Rebin(4);
	 //~h[i]->Smooth();

	gPad->Update();
	 //~h[i]->GetXaxis()->SetRange(h[i]->FindBin(400),h[i]->FindBin(700));
	 //~gPad->Update();
	 //~cout<<h[i]->GetMean()<<endl; 
	 }


  bool fit=false;
  //~bool fit=true;

TF1 *f[n] ;  //gaussianas para el ajuste
if (fit==true){
 //Patricia code
 ofstream ofs("Gaus_Results.txt",std::ofstream::out | std::ofstream::app);
 ofs <<"----------------------------------------------"<<endl;

	for(int i=0;i<n;++i){
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


void DrawSum(std::vector<int> mych, int n, int runs[], string name[], bool rate=true, string var="Q3PE")
{
 TFile *ifile[n];
 for(int i=0;i<n;i++) ifile[i]= new TFile(Form("AnalysisROOT/Run%i_NTuple.root",runs[i])); 
 TNtuple *nt[n];
 Double_t duration[n];
 Float_t evt;
 Float_t ch; 
 Float_t time; 
 Float_t Q3; 
 std::vector<std::map<int,double>> EvtToCharge(n); 
  TH1D *h[n];
 for(int i=0;i<n;i++)
 {
   nt[i]= (TNtuple*)ifile[i]->Get("ntuple");
   nt[i]->SetBranchAddress("time",&time);
   nt[i]->SetBranchAddress("evt",&evt);
   nt[i]->SetBranchAddress("ch",&ch);
   nt[i]->SetBranchAddress(var.c_str(),&Q3);
   if(rate) h[i]=new TH1D(Form("h%i",i),Form("%s;Charge collected on ch0 and 1 (PE); Events/s ",name[i].c_str()),500,0,100);
   else h[i]=new TH1D(Form("h%i",i),Form("%s;Charge collected on ch0 and 1 (PE); Events",name[i].c_str()),500,0,100);
   h[i]->SetLineColor(i+1);if (i==4) h[i]->SetLineColor(i+3);//yellow sucks hard
   h[i]->SetLineWidth(2);
   //EvtToCharge[i]=std::map<int,double>();
 }
 cout << "Drawing " << n<< endl;
 for(int i=0;i<n;i++)
 {
// cout << "Drawing " << i<< endl;
  for(int j=0;j<nt[i]->GetEntries();j++)
  { 
//   cout << "Drawing " << j<< endl;
   nt[i]->GetEntry(j);
   if(find(mych.begin(),mych.end(),(int)ch)==mych.end()) continue;
   if(EvtToCharge[i].find((int)evt) == EvtToCharge[i].end() )
   {
     EvtToCharge[i].emplace((int)evt,Q3);
   }
   else
   {
     EvtToCharge[i][(int)evt]+=(Q3);
   }
  }

  float firstime=time;
  nt[i]->GetEntry(nt[i]->GetEntries()-1);
  float lasttime=time;
  duration[i] = (lasttime - firstime)*8.e-9;
 }
 for(int i=0;i<n;i++)
 {
   for (auto e : EvtToCharge[i]){h[i]->Fill(e.second);}
   if(rate) h[i]->Scale(1./duration[i]);
   //else h[i]->Scale(1./h[i]->GetEntries());
 }

 cout << "Drawing4 " << n<< endl;

 TCanvas *c = new TCanvas("c");
 c->cd();
 for(int i=0;i<n;i++){h[i]->Rebin();h[i]->Draw("HIST SAME");}
 gPad->BuildLegend();
 lets_pause();
 
 //Patricia code
 ofstream ofs("Gaus_Results.txt");
	for(int i=0;i<n;++i){
			h[i]->GetYaxis()->SetRangeUser(1e-5,1);
			double min=30,max=90;
		   TF1 *f = new TF1("f","gaus",min,max);  //ajuste a una gaussiana
		   h[i]->Fit(f,"R","SAME",min,max);
		   f->Draw("SAME");f->Print();
		   lets_pause(); 
		   ofs << runs[i] << "\t" << f->GetParameter(1) << "\t" << f->GetParError(1)<< endl;
       
    } 
}

void ChargeSpectrum()
{
//(std::vector<int> mych, std::vector<int> runs, std::vector<string> name, bool rate=true, string var="Q3", TVirtualPad *pad=NULL, float CustomMaxRange=0 )

//~Draw({0,0,0}, {15,27,39},{"26/05","27/05","28/05"},false,{"Q1PE","Q1PE","Q1PE"});lets_pause();
//~Draw({1,1,1}, {15,27,39},{"26/05","27/05","28/05"},false,{"Q1PE","Q1PE","Q1PE"});lets_pause();
//~Draw_April({1,1,1,1,1}, {12,26,38,50,62},{"30/04","01/05","02/05","03/05","04/05"},false,{"Q1PE","Q1PE","Q1PE","Q1PE","Q1PE"});lets_pause();
//~Draw_April({0,0,0,0,0}, {12,26,38,50,62},{"30/04","01/05","02/05","03/05","04/05"},false,{"Q1PE","Q1PE","Q1PE","Q1PE","Q1PE"});lets_pause();

//~Draw_April({0}, {12},{"QTotal"},false,{"QTotal"});lets_pause();
Draw({0}, {15},{"QTotal"},false,{"QTotal"});lets_pause();

//~//SiPM coin vs auto trigger, Amp
//~Draw_April({0,0,0,0,0}, {12,12,12,12,12}, {"5#mus Tras pico","2.5#mus Tras pico","0.5 #mus Tras pico","Pico hasta Pedestal","Total Charge in event"}, true,{"Q1","Q2","Q3","QPeak","QTotal"}); lets_pause();

//~Draw_April({1,1}, {12,12},{"TStartQpeak","TEndQpeak"},false,{"TStartQpeak","TEndQPeak"});lets_pause();
//~Draw_April({0,0}, {12,12},{"TStartQpeak","TEndQpeak"},false,{"TStartQpeak","TEndQPeak"});lets_pause();
//~Draw_April({0,0}, {12,12},{"TStartQpeak","TEndQpeak"},false,{"TStartQpeak","TEndQPeak"});lets_pause();
//~Draw({1}, {15},{"TDiff"},false,{"TDiff"});lets_pause();
//~Draw({0}, {15},{"TDiff"},false,{"TDiff"});lets_pause();
//~Draw({0,0}, {15,15},{"TStartQpeak","TEndQpeak"},false,{"TStartQpeak","TEndQPeak"});lets_pause();
//~Draw({1,1}, {15,15},{"TStartQpeak","TEndQpeak"},false,{"TStartQpeak","TEndQPeak"});lets_pause();
//~Draw({0}, {15},{"TDiff"},false,{"TDiff"});lets_pause();
//--------------------------------------------
//~//SiPM coin vs auto trigger, Amp
//~Draw({0,0,0,0,0}, {15,15,15,15,15}, {"5#mus Tras pico","2.5#mus Tras pico","0.5 #mus Tras pico","Pico hasta Pedestal","Total Charge in event"}, true,{"Q1","Q2","Q3","QPeak","QTotal"}); lets_pause();
//--------------------------------------------
//~//SiPM coin vs auto trigger, Amp
//~Draw({1,1,1,1,1}, {15,15,15,15,15}, {"5#mus Tras pico","2.5#mus Tras pico","0.5 #mus Tras pico","Pico hasta Pedestal","Total Charge in event"}, true,{"Q1","Q2","Q3","QPeak","QTotal"}); lets_pause();
//--------------------------------------------
//~//SiPM coin vs auto trigger, Amp
//~Draw({0,0}, {14,15}, {"Coincidencias 25","Coincidencias 200"}, true,{"Q3","Q3"}); lets_pause();
//--------------------------------------------
//~//SiPM coin vs auto trigger, Amp
//~Draw({1,1,1}, {5,14,15}, {"AutoTrigger","Coincidencias 25","Coincidencias 200"}, true,{"Amp","Amp","Amp"}); lets_pause();
//--------------------------------------------
//~//SiPM coin vs auto trigger, Amp
//~Draw({1,1,1}, {5,14,15}, {"AutoTrigger","Coincidencias 25","Coincidencias 200"}, true,{"Amp","Amp","Amp"}); lets_pause();
//--------------------------------------------

//SiPM Solo trigger vs coincidencias 25 y 200, ambos canales
//~Draw({1,2,1,2,1,2}, {5,5,11,11,12,12}, {"Th25-Ch4 SingleTrigger","Th25-Ch5 SingleTrigger","Th25-Ch4 Coincidencias ","Th25-Ch5 Coincidencias","Th200-Ch4 Coincidencias","Th200-Ch5 Coincidencias"}, true,{"Q3PE","Q3PE","Q3PE","Q3PE","Q3PE","Q3PE"}); lets_pause();

//SiPM autotrigger 4-25, charge
//~Draw({1,1,1,1,1}, {5,19,31,43,55}, {"04/30","05/01","05/02","05/03","05/04"}, true,{"Q3PE","Q3PE","Q3PE","Q3PE","Q3PE"}); lets_pause();
//~//--------------------------------------------
//~//PMT coin 25, charge
//~Draw({0,0,0,0,0}, {11,25,37,49,61}, {"04/30","05/01","05/02","05/03","05/04"}, true,{"Q3PE","Q3PE","Q3PE","Q3PE","Q3PE"}); lets_pause();
//~//--------------------------------------------
//~//SiPM coin 25, charge
//~Draw({1,1,1,1,1}, {11,25,37,49,61}, {"04/30","05/01","05/02","05/03","05/04"}, true,{"Q3PE","Q3PE","Q3PE","Q3PE","Q3PE"}); lets_pause();
//--------------------------------------------
//~//PMT coin vs auto trigger, charge
//~Draw({0,0}, {11,5}, {"Coincidence trigger","Single channel trigger"}, true,{"Q3PE","Q3PE"}); lets_pause();
//~//--------------------------------------------
//SiPM coin vs auto trigger, charge
//~Draw({1,1}, {11,5}, {"Coincidence trigger","Single channel trigger"}, true,{"Q3PE","Q3PE"}); lets_pause();
//~//--------------------------------------------
//~//SiPM coin vs auto trigger, charge
//~Draw({0,0}, {5,11}, {"Auto-Trigger","Coincidence Trigger"}, true,{"Q3PE","Q3PE"}); lets_pause();
//---------------------------------------------
//SiPM coin 200, charge
//~Draw({1,1,1,1,1,2,2,2,2,2}, {12,26,38,50,62,12,26,38,50,62}, {"","","","","","","","","",""}, false,{"Q3PE","Q3PE","Q3PE","Q3PE","Q3PE","Q3PE","Q3PE","Q3PE","Q3PE","Q3PE"}); lets_pause();
//~Draw({1,1,1,1,1,2,2,2,2,2}, {6,20,32,44,56,9,23,35,47,59}, {"","","","","","","","","",""}, false,{"Q3PE","Q3PE","Q3PE","Q3PE","Q3PE","Q3PE","Q3PE","Q3PE","Q3PE","Q3PE"}); lets_pause();
//----------------------------------------------
//~//PMT coin 200, 
//~Draw({1}, {12}, {"SIPM Coincidencias"}, false,{"Amp"}); gPad->SetGrid(1,1);lets_pause();
//~Draw({1}, {12}, {"SIPM Coincidencias"}, false,{"Q3"}); gPad->SetGrid(1,1);lets_pause();
//~Draw({1}, {12}, {"SIPM Coincidencias"}, false,{"AmpPE"});gPad->SetGrid(1,1); lets_pause();
//~Draw({1}, {12}, {"SIPM Coincidencias"}, false,{"Q3PE"}); gPad->SetGrid(1,1);lets_pause();

//~Draw({0,0}, {12,12}, {"PMT Coincidencias",""}, true,{"Q3","QPeak"}); gPad->SetGrid(1,1);lets_pause();
//~Draw({0,0,0,0,0}, {12,26,38,50,62}, {"","","","",""}, false,{"Amp","Amp","Amp","Amp","Amp"}); lets_pause();
//~Draw({0,0,0,0,0}, {12,26,38,50,62}, {"","","","",""}, false,{"AmpPE","AmpPE","AmpPE","AmpPE","AmpPE"}); lets_pause();
//~Draw({0,0,0,0,0}, {12,26,38,50,62}, {"","","","",""}, false,{"Q3","Q3","Q3","Q3","Q3"}); lets_pause();
//~Draw({0,0,0,0,0}, {2,3,4,11,12}, {"PMT-25Thr","PMT-200 Thr","PMT-600 Thr","Coincidencias-25 Thr","Coincidencias-200 Thr"}, true,{"Q3PE","Q3PE","Q3PE","Q3PE","Q3PE"}); lets_pause();

//SiPM coin 200, Amp
//~Draw({1,1,1,1,1,2,2,2,2,2}, {12,26,38,50,62,12,26,38,50,62}, {"","","","","","","","","",""}, false,{"AmpPE","AmpPE","AmpPE","AmpPE","AmpPE","AmpPE","AmpPE","AmpPE","AmpPE","AmpPE"}); lets_pause();
//~Draw({1,1,1,1,1,2,2,2,2,2}, {6,20,32,44,56,9,23,35,47,59}, {"","","","","","","","","",""}, false,{"AmpPE","AmpPE","AmpPE","AmpPE","AmpPE","AmpPE","AmpPE","AmpPE","AmpPE","AmpPE"}); lets_pause();

//PMT 900
//~Draw({0,0,0,0}, {4,18,30,42}, {"PMT 30Abr","PMT 1 Mayo","PMT 2Mayo","PMT 3Mayo"}, true,{"Amp", "Amp", "Amp", "Amp"}); lets_pause();

}
