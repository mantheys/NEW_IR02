/*
Macro para dibujar el perfil de centelleo de distintos runes superpuesto uno sobre otro.
Toma los perfiles ya generados a partir de la macro AverageWaveform.C, y los plotea juntos en el mismo canvas.
*/
#include"lib/headers.h"
#include"lib/ThreeExpoGausFitScintUtilsCh_MultiRanges.h"
TH1D* CenterHist(TH1D* h, int Rebin=1, double bincenter=1.3e-6)
{
  h->Rebin(Rebin);
  h->Smooth();
  cout << "Entro " << endl;
  int binmax = h->GetMaximumBin();
  double xmax = h->GetXaxis()->GetBinCenter(binmax);
  double ymax = h->GetBinContent(binmax);
  cout << xmax << " " << ymax << endl;

  //binmax = h->GetMaximumBin();
  ymax = h->GetBinContent(binmax);
  binmax = h->FindFirstBinAbove(0.8*ymax);
  //bincenter=h->FindBin(1.5e-6);
  //binmax=bincenter;

  //~ymax=1;
	
  int binshift = h->FindBin(bincenter);
  TH1D *aux=(TH1D*)h->Clone("haux"); aux->Reset();
  for(int i=0; i<aux->GetSize();i++) aux->SetBinContent(i-binmax+binshift,(h->GetBinContent(i)/ymax));
  //~aux->Scale(1/aux->GetBinContent(aux->FindBin(0.5e-6)));
  //~cout<<aux->GetBinContent(aux->FindBin(0.5e-6))<<endl;
  cout<<aux->Integral(aux->FindFixBin(0.25e-6), aux->FindFixBin(0.4e-6), "")<<endl;
  cout << "Salgo <" << endl;
  return aux;
  
 
}
void Draw(int n, int *runs, string* name, int pm)
{
ana::HistogramCollection_t mycol[n];
 for(int i=0;i<n;i++) mycol[i].GetFromFile(Form("AnalysisROOT/Run%i_ScintProfFirstSignalBin.root",runs[i]));
 cout << "Hist loaded " << endl;
  TH1 *h[n];
 for(int i=0;i<n;i++)
 {
   h[i]=(TH1D*)(mycol[i].GetByOpChannel(pm))[0];
//   h[i]=CenterHist((TH1D*)(mycol[i].GetByOpChannel(pm))[0]);
   h[i]->SetLineColor(i+1);
   if (i==5) h[i]->SetLineColor(i+3);//yellow sucks hard
   h[i]->SetLineWidth(2);
   h[i]->SetTitle(name[i].c_str());
 }
 cout << "Hist loaded " << endl;
 TCanvas *c = new TCanvas("c");
 c->cd();
 for(int i=0;i<n;i++) h[i]->Draw("HIST SAME");
 gPad->BuildLegend();
 lets_pause();

  
}

void Clean(TH1D* h)
{
  cout << "Cleaining" << endl;
  for(int i=1; i<=h->GetSize();i++) if(h->GetBinCenter(i)<3.5e-6) h->SetBinContent(i,0.0);
}
void Draw2(int n, string *runs, string* name, int* pm,int *Rebin)
{
ana::HistogramCollection_t mycol[n];
 for(int i=0;i<n;i++) mycol[i].GetFromFile(runs[i].c_str());

 cout << "Hist loaded " << endl;
  TH1 *h[n];
 for(int i=0;i<n;i++)
 {
   TH1D* h1=(TH1D*)(mycol[i].GetByOpChannel(pm[i]))[0];
   if(runs[i]==13) Clean(h1);// h1->Draw("HIST"); lets_pause();
   h[i]=CenterHist(h1);


         ThreeExpoGausFitScintUtilsCh ChiaraFunction; 
         ChiaraFunction.Fit(*h[i]); TF1 *newf;
         h[i]->Draw("HIST");gStyle->SetOptFit(1112);gStyle->SetOptStat(0); double xmin, xmax;
         if(h[i]->GetListOfFunctions()->FindObject("scint_time"))
         {
           h[i]->GetListOfFunctions()->FindObject("scint_time")->Draw("SAME");
           ((TF1*)h[i]->GetListOfFunctions()->FindObject("scint_time"))->GetRange(xmin,xmax); 
           h[i]->GetXaxis()->SetRangeUser(-0.1e-6+xmin,xmax+0.6e-6); h[i]->GetYaxis()->SetRangeUser(0.00005,3);
//         TF1* f=(TF1*)h->GetListOfFunctions()->FindObject("scint_time"); cout << "scint_time FOUND!" << endl;
//         newf=f->Clone("scint_time2");
           TPaveStats *st = (TPaveStats*)h[i]->FindObject("stats");
           st->SetY1NDC(0.55); //new x start position
           st->SetY2NDC(0.9); //new x end position
           lets_pause();
         }
         



   h[i]->SetLineColor(i+1);
   if (i==5) h[i]->SetLineColor(i+4);//yellow sucks hard
   h[i]->SetLineWidth(2);
   h[i]->SetTitle(Form("%s - #tau = %.2f #mus",name[i].c_str(),1.e6*((TF1*)h[i]->GetListOfFunctions()->FindObject("scint_time"))->GetParameter(8)));
 }
 cout << "Hist loaded " << endl;
 TCanvas *c = new TCanvas("c");gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
 c->cd();
 for(int i=0;i<n;i++){ h[i]->Draw("HIST SAME");}
 gPad->BuildLegend();
 for(int i=0;i<n;i++){ h[i]->GetListOfFunctions()->FindObject("scint_time")->Draw("SAME");}
 lets_pause();  
}

void Draw3(int n, string *runs, string* name, int* pm,int *Rebin)
{
ana::HistogramCollection_t mycol[n];
 for(int i=0;i<n;i++) mycol[i].GetFromFile(runs[i].c_str());

 cout << "Hist loaded " << endl;
  TH1 *h[n];
 for(int i=0;i<n;i++)
 {
   TH1D* h1=(TH1D*)(mycol[i].GetByOpChannel(pm[i]))[0];
   //Clean(h1);
   h[i]=CenterHist(h1,Rebin[i]);   h[i]->SetLineColor(i+1);   if (i==4) h[i]->SetLineColor(i+3);//yellow sucks hard
   h[i]->SetLineWidth(2);
   h[i]->SetTitle(Form("%s",name[i].c_str()));
   //~h[i]->Smooth();
 }
 cout << "Hist loaded " << endl;
 TCanvas *c = new TCanvas("c");gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
 c->cd();
 for(int i=0;i<n;i++) {
	 //~h[i]->Smooth();
	 h[i]->Draw("HIST SAME C");
	 	 }


bool fit=true;
//~bool fit=true;
TF1 *f[n] ;  //exp para el ajuste
if (fit==true){
 //Patricia code
 ofstream ofs("Gaus_Results.txt",std::ofstream::out | std::ofstream::app);
 ofs <<"----------------------------------------------"<<endl;

	for(int i=0;i<n;++i){
			h[i]->GetYaxis()->SetRangeUser(1e-5,1.2);
			h[i]->GetXaxis()->SetRangeUser(0.7e-6,4e-6);
			double min=2.0e-6,max=3e-6;
		   f[i] = new TF1("f","exp([0]-(x)*1000000/[1])",min,max);
		   f[i]->SetParameter(0,-2.95);
		   f[i]->SetParameter(1,1/1.47402);
		   f[i]->SetParName(1,"#tau_{slow}");

		   //~f[i] = new TF1("f","expo",min,max);
		   h[i]->Fit(f[i],"L R","SAME",min,max);
		   f[0]->Draw("SAME");f[i]->Print();gPad->Update();
		   //~lets_pause(); 
		   cout << runs[i] << "\t" << f[i]->GetParameter(1) << "\t" << f[i]->GetParameter(2)<<"\t"<< endl;
       
    } 
	//~cout<<"{";	
    //~for(int i=0;i<n-1;++i) cout<<f[i]->Integral(-1000,1000)/h[i]->GetBinWidth(1)<<",";
	//~cout<<f[n-1]->Integral(-1000,1000)/h[n-1]->GetBinWidth(1);
	//~cout<<"}"<<endl;	
   gStyle->SetOptFit();

}
    // draw the legend
   TLegend *legend=new TLegend(0.7, 0.6, 0.9, 0.9);
   //~legend->SetTextFont(72);
   //~legend->SetTextSize(0.04);
   for(int i=0;i<n;++i){
		h[i]->GetXaxis()->SetRangeUser(0.5e-6,4e-6);
		legend->AddEntry(h[i],Form("%s - #tau_{slow}= %.2f #mus",name[i].c_str(),f[i]->GetParameter(1)));}//pmts ajuste exponencial
		//~legend->AddEntry(h[i],name[i].c_str());}/ leyenda sin ajuste exponencial

   if (fit==true)legend->AddEntry(f[0],"Exponential Fit","l");
   legend->Draw();
   gPad->SetLogy();
	//~h[i]->GetYaxis()->SetRangeUser(1e-5,1.2);

TPaveText *t = new TPaveText(0.3, 0.85, 0.7, 0.95, "brNDC"); // middle
t->AddText("PMT Scintillation Profile");
t->Draw();

 lets_pause();  
}

void Scintillation()
{
/*Argumentos:
1. número de runes a plotear
2. vector de enteros con los runes a plotear
3. vector de strings, con los títulos de los runes para poner en la leyenda.
4. Canal a dibujar.
*/

// SiPM 200 Coincidencias
	const int n = 2;
	
   string name[n]={
		"Abril",
		"Mayo"};
                  
   string runs[n]={
		"AnalysisROOT_PMT_APRIL/Run12_ScintProfFirstSignalBin.root",
		"AnalysisROOT/Run15_ScintProfFirstSignalBin.root"};
   
   int pms[n]={0,0};
   int rebin [n]={1,1};
   Draw3(n,runs,name,pms,rebin);

//---------------------------------------------------------------------------
//~// SiPM 200 Coincidencias
	//~const int n = 2;
	
   //~string name[n]={
		//~"Primer pico",
		//~"Sgundo pico"};
                  
   //~string runs[n]={
		//~"AnalysisROOT/Run15_ScintProf_PMT_primera_gaussiana.root",
		//~"AnalysisROOT/Run15_ScintProf_PMT_segunda_gaussiana.root"};
   
   //~int pms[n]={0,0};
   //~int rebin [n]={1,1};
   //~Draw3(n,runs,name,pms,rebin);

//~//---------------------------------------------------------------------------
// SiPM 200 Coincidencias
	//~const int n = 4;
	
   //~string name[n]={
		//~"05/26",
		//~"05/27",
		//~"05/28",
		//~"Laser"};
                  
   //~string runs[n]={
		//~"AnalysisROOT/Run15_ScintProfFirstSignalBin.root",
		//~"AnalysisROOT/Run27_ScintProfFirstSignalBin.root",
		//~"AnalysisROOT/Run39_ScintProfFirstSignalBin.root",
		//~"AnalysisROOT/Run1002_ScintProfFirstSignalBin.root"};
   
   //~int pms[n]={0,0,0,0};
   //~int rebin [n]={1,1,1,1};
   //~Draw3(n,runs,name,pms,rebin);

//---------------------------------------------------------------------------
//~// SiPM 200 (PMT channel)
	//~const int n = 5;
	
   //~string name[n]={
		//~"04/30",
        //~"05/01",
		//~"05/02",
		//~"05/03",
		//~"05/04"};
		//~//"Laser"};
                  
   //~string runs[n]={
		//~"AnalysisROOT_PMT/Run12_ScintProfFirstSignalBin.root",
		//~"AnalysisROOT_PMT/Run26_ScintProfFirstSignalBin.root",
		//~"AnalysisROOT_PMT/Run38_ScintProfFirstSignalBin.root",
		//~"AnalysisROOT_PMT/Run50_ScintProfFirstSignalBin.root",
		//~"AnalysisROOT_PMT/Run62_ScintProfFirstSignalBin.root"};
		//~//"AnalysisROOT/Run1001_ScintProfFirstSignalBin.root"};
   
   //~int pms[n]={0,0,0,0,0};
   //~int rebin [n]={1,1,1,1,1};
   //~Draw3(n,runs,name,pms,rebin);

//---------------------------------------------------------------------------
// PMT 200
/*
	const int n = 5;
	
   string name[n]={
		"Run 12 - 30/04",
        "Run 26 - 01/05",
		"Run 38 - 02/05",
		"Run 50 - 03/05",
		"Run 62 - 04/05"};
                  
   string runs[n]={
		"AnalysisROOT_PMT/Run12_ScintProfFirstSignalBin.root",
		"AnalysisROOT_PMT/Run26_ScintProfFirstSignalBin.root",
		"AnalysisROOT_PMT/Run38_ScintProfFirstSignalBin.root",
		"AnalysisROOT_PMT/Run50_ScintProfFirstSignalBin.root",
		"AnalysisROOT_PMT/Run62_ScintProfFirstSignalBin.root"};
   
   int pms[n]={0,0,0,0,0};
   int rebin [n]={1,1,1,1,1};
   Draw3(n,runs,name,pms,rebin);
//*/
//---------------------------------------------------------------------------
// PMT Feb
/*
	const int n = 5;
	
   string name[n]={
		"Run 12 - 11/02",
        "Run 23 - 12/02",
		"Run 34 - 13/02",
		"Run 45 - 14/02",
		"Run 56 - 15/02"};
                  
   string runs[n]={
		"AnalysisROOT_PMT_FEB/Run12_ScintProfFirstSignalBin.root",
		"AnalysisROOT_PMT_FEB/Run23_ScintProfFirstSignalBin.root",
		"AnalysisROOT_PMT_FEB/Run34_ScintProfFirstSignalBin.root",
		"AnalysisROOT_PMT_FEB/Run45_ScintProfFirstSignalBin.root",
		"AnalysisROOT_PMT_FEB/Run56_ScintProfFirstSignalBin.root"};
   
   int pms[n]={2,2,2,2,2};
   int rebin [n]={1,1,1,1,1};
   Draw3(n,runs,name,pms,rebin);
//*/
//-----------------------------------------------------------------------------------
//~// PMT 200 Feb Vs April

	//~const int n = 4;
	
   //~string name[n]={
		//~"Abril - 30/04",
		//~"Abril - 04/05",
		//~"Feb - 11/02",
		//~"Feb - 15/02"};
                  
   //~string runs[n]={
		//~"AnalysisROOT_PMT/Run12_ScintProfFirstSignalBin.root",
		//~"AnalysisROOT_PMT/Run62_ScintProfFirstSignalBin.root",
		//~"AnalysisROOT_PMT_FEB/Run12_ScintProfFirstSignalBin.root",
		//~"AnalysisROOT_PMT_FEB/Run56_ScintProfFirstSignalBin.root"};
   
   //~int pms[n]={0,0,2,2};
   //~int rebin [n]={1,1,1,1};
   //~Draw3(n,runs,name,pms,rebin);
//~//*/
//-----------------------------------------------------------------------------------
// SiPM 200 Feb Vs April
/*
	const int n = 4;
	
   string name[n]={
		"Abril - 30/04",
		"Abril - 04/05",
		"Feb - 11/02",
		"Feb - 15/02"};
                  
   string runs[n]={
		"AnalysisROOT/Run12_ScintProfFirstSignalBin.root",
		"AnalysisROOT/Run62_ScintProfFirstSignalBin.root",
		"AnalysisROOT_JASOTO_FEB/Run12_ScintProfFirstSignalBin_alpha.root",
		"AnalysisROOT_JASOTO_FEB/Run56_ScintProfFirstSignalBin.root"};
   
   int pms[n]={1,1,1,1};
   int rebin [n]={1,1,1,1};
   Draw3(n,runs,name,pms,rebin);
//*/
//--------------------------------------------------------------------------------------
// SiPM SPE Feb Abril
/*
	const int n = 4;
	
   string name[n]={
		"Abril - SPE Alfas",
		"Abril - SPE Calibracion",
		"Feb - Alfas SPE",
		"Feb - AvWf_SPE"};
                  
   string runs[n]={
		"AnalysisROOT/Run19_ScintProfFirstSignalBin_SPE.root",
		"AnalysisROOT/Run2_ScintProfFirstSignalBin_SPE_cal.root",
		"AnalysisROOT_JASOTO_FEB/Run3_SPE_ScintProfFirstSignalBin.root",
		"AnalysisROOT_JASOTO_FEB/AvWf_SPE.root"};
   
   int pms[n]={1,0,0,0};
   int rebin [n]={1,1,1,1};
   Draw3(n,runs,name,pms,rebin);
//*/







//   string name[5]={"Run 12 - 11/02","Run 23 - 12/02","Run 34 - 13/02", "Run 45 - 14/02", "Run 56 - 15/02"};
//   int runs[5]={12,23,34,45,56};
//   int pms[5]={1,1,1,1,1};
//   Draw2( 5,runs,name,pms,10,"Deconvoluted");

/*
   string name[5]={"Run 12 - 11/02","Run 23 - 12/02","Run 34 - 13/02", "Run 45 - 14/02", "Run 56 - 15/02"};
   int runs[5]={12,23,34,45,56};
   int pms[5]={1,1,1,1,1};
   Draw2( 5,runs,name,pms,1);
*/
//   string name[3]={"Run 12 - (0,1) 35ADC 11/02","Run 13 - (0,1) 750ADC - 11/02","LASER - 11/02"};
//   int runs[3]={12,13,1001};
//   int pms[3]={1,1,1};
//   Draw2( 3,runs,name,pms,1);
/*
   string name[2]={"Run 10 - (2) ALFA","(2) LED - 11/02"};
   string runs[2]={"AnalysisROOT/Run11_ScintProfFirstSignalBin.root",
                   "AnalysisROOT/Run1001_ScintProfFirstSignalBin.root"};
   int pms[2]={2,2};
   int rebin[2]={1,1};
   Draw3( 2,runs,name,pms,rebin);
//*/
/*
   string name[2]={"Run 13 - SPE ALPHA - 11/02","SPE Laser - 11/02"};
   int runs[2]={13,500};
   int pms[2]={1,1};
   int rebin[2]={1,1};
   Draw2( 2,runs,name,pms,1,"SPEPeak");
*/
/*
   string name[2]={"Run 13 - SPE ALPHA - 11/02 - FirstBin","Run 13 - SPE ALPHA - 11/02 - Peak"};
   string runs[2]={"AnalysisROOT/Run13_SPE_ScintProfBeginPeak.root",
                "AnalysisROOT/Run13_SPE_ScintProf.root"};
   int pms[2]={0,0};
   Draw3( 2,runs,name,pms);
*/

/*
   string name[2]={"SPE Laser Osciloscopio","SPE Laser ADC"};
   string runs[2]={"AnalysisROOT/AvWf_SPE_Osc.root",
                   "AnalysisROOT/Run500_SPEPeak.root"};
   int pms[2]={1,1};
   int rebin[2]={2,1};
   Draw3( 2,runs,name,pms,rebin);
//*/
/*
   string name[2]={"Run 3 - SPE ALPHA - 11/02","SPE Laser - 11/02"};
   string runs[2]={"AnalysisROOT/Run3_SPE_ScintProfFirstSignalBin.root",
                   "AnalysisROOT/Run500_SPEPeak.root"};
   int pms[2]={0,0};
   int rebin[2]={1,1};
   Draw3( 2,runs,name,pms, rebin);
//*/
/*
   string name[2]={"Run 13 - MUONS - Amp>1100ADC - 11/02","Run 12 - ALPHAS - 250<Amp<1100ADC - 11/02"};
   string runs[2]={"AnalysisROOT/Run13_ScintProfFirstSignalBin_muon.root",
                   "AnalysisROOT/Run12_ScintProfFirstSignalBin_alpha.root"};
   int pms[2]={0,0};
   int rebin[2]={1,1};
   Draw3( 2,runs,name,pms,rebin);
//*/


/*
   string name[3]={"Run 12 ADC ALPHA - 11/02","Run 56 - ADC ALPHA 15/02", "Run 2002 - OSC 16/02"};
   int runs[3]={12,56,2002};
   int pms[3]={1,1,1};
   int rebin[3]={1,1,1};
   Draw2( 3,runs,name,pms,1);
//*/
/*
   string name[3]={"Run 10 ADC ALPHA - 11/02","Run 54 - ADC ALPHA 15/02", "Run 2002 - OSC 16/02"};
   int runs[3]={10,54,2002};
   int pms[3]={2,2,2};
   Draw2( 3,runs,name,pms,1);
*/

/*
   string name[2]={"More light ALPHA (ADC)","More light Laser (ADC)"};
   string runs[2]={"AnalysisROOT/Run2_ScintProfFirstSignalBin.root",
                   "AnalysisROOT/Run1001_ScintProfFirstSignalBin.root"};
   int pms[2]={0,0};
   int rebin[2]={1,1};
   Draw3( 2,runs,name,pms, rebin);
*/

/*
   string name[2]={"More light Alfa (Osciloscope)","More light Laser (Osciloscope)"};
   string runs[2]={"AnalysisROOT/Run2002_ScintProfFirstSignalBin.root",
                   "AnalysisROOT/Run2004_ScintProfFirstSignalBin.root"};
   int pms[2]={0,0};
   int rebin[2]={1,1};
   Draw3( 2,runs,name,pms, rebin);
*/
/*
   string name[1]={"SPE Laser Osciloscopio"};
   string runs[1]={"AnalysisROOT/AvWf_SPE_Osc.root"};
   int pms[1]={1};
   int rebin[1]={1};
   Draw3( 1,runs,name,pms,rebin);
*/
/*
   string name[2]={"ALPHA","MUON"};
   string runs[2]={"AnalysisROOT/Run12_ScintProfFirstSignalBin_alpha_PMT.root",
                   "AnalysisROOT/Run13_ScintProfFirstSignalBin_muon_PMT.root"};
   int pms[2]={2,2};
   int rebin[2]={1,1};
   Draw2( 2,runs,name,pms, rebin);
*/
/*
   string name[2]={"1st Gaussian","2nd Gaussian"};
   string runs[2]={"AnalysisROOT/Run11_ScintProfFirstSignalBin_1stGaussian_PMT2.root",
                   "AnalysisROOT/Run11_ScintProfFirstSignalBin_2ndGaussian_PMT2.root"};
   int pms[2]={2,2};
   int rebin[2]={1,1};
   Draw2( 2,runs,name,pms, rebin);


//*/
/*

   string name[2]={"SPE alfa Osciloscopio","SPE Laser Osciloscopio"};
   string runs[2]={"AnalysisROOT/Run2001_AvWf_SPE_Osc.root",
                   "AnalysisROOT/AvWf_SPE_Osc.root"};
   int pms[2]={1,1};
   int rebin[2]={1,1};
   Draw3( 2,runs,name,pms,rebin);


   string name[1]={"PMT LED Osc"};
   string runs[1]={"AnalysisROOT/Run2003_ScintProfFirstSignalBin_PMT.root"};
   int pms[1]={0};
   int rebin[1]={1};
   Draw3( 1,runs,name,pms,rebin);

   string name[1]={"PMT ALFA Osc"};
   string runs[1]={"AnalysisROOT/Run2002_ScintProfFirstSignalBin_PMT.root"};
   int pms[1]={0};
   int rebin[1]={1};
   Draw3( 1,runs,name,pms,rebin);
*/
/*
   string name[1]={"PMT LED SPE Osc"};
   string runs[1]={"AnalysisROOT/RunCaliPMTOsci_ScintProf_PMT.root"};
   int pms[1]={0};
   int rebin[1]={1};
   Draw3( 1,runs,name,pms,rebin);

   string name[1]={"PMT ALFA SPE ADC"};
   string runs[1]={"AnalysisROOT/Run9_ScintProf_SPE_PMT.root"};
   int pms[1]={2};
   int rebin[1]={1};
   Draw3( 1,runs,name,pms,rebin);


   string name[2]={"PMT LED Osc", "SiPM Laser"};
   string runs[2]={"AnalysisROOT/Run2003_ScintProfFirstSignalBin_PMT.root",
                   "AnalysisROOT/Run2004_ScintProfFirstSignalBin.root"};
   int pms[2]={0,0};
   int rebin[2]={1,1};
   Draw3( 2,runs,name,pms,rebin);

*/
/*
   string name[2]={"SiPM ALFA Osc","PMT ALFA Osc"};
   string runs[2]={"AnalysisROOT/Run2002_ScintProfFirstSignalBin.root"
                  ,"AnalysisROOT/Run2002_ScintProfFirstSignalBin_PMT.root"};
   int pms[2]={0,0};
   int rebin[2]={1,1};
   Draw3( 2,runs,name,pms,rebin);
*/
/*
{10,{2},-1},
{21,{2},-1},
{32,{2},-1},
{43,{2},-1},
{54,{2},-1},
*/

}
