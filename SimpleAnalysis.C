/*
Esta macro vuelca analiza de forma rápida un run.
1. Plotea los pedestales, la pedestal STD, el peaktime
2. Aplica unos cortes básicos para quitar eventos ruidosos.
3. Imprime por consola y a un fichero SimpleAnalysis.log la carga y amplitud media detectada por canal a priori, y pidiendo coincidencias.
4. Vuelca el espectro de carga del run a un fichero AnalysisROOT/Run*_Q3MaxPeak500ns.root
5. Vuelca a una ntupla las variables de carga y amplitud por evento y por canal en el fichero AnalysisROOT/Run%i_NTuple.root
*/
#include"lib/headers.h"


//ofstream ofile("Run.log");
void Analyse(int r, double range1, double range2, std::vector<int> triggerchannels, ofstream &ofs, ofstream &ofs2,  int Nevts=-1, string adc="DT5725" )
{
  cout << "run " << r << " Trigger channels "; for(auto t : triggerchannels) cout << t << " "; cout << endl;
  
  //TODO: Ganancias de las calibraciones de Febrero: corregir cuando se pueda
  std::vector<double> Gains={0.56,
                             3.61,
                             3.61,
                             3.83}; //ganancias en pC
  std::vector<double> SPEAmp={38.6,
                              24.8,
                              24.8,
                              25.5};//Amplitud del SPE en cuentas de ADC
  for(int i=0;i<4;i++) Gains[i]=Gains[i]/1.602e-7; //(1e-19+12)
  ana::Run_t myrun(r,{{Form("/pc/choozdsk01/palomare/SiPM/SC_Fuente_Alpha_Oct/ROOT/run%i_ch2.root",r),"ADC0"},
                      {Form("/pc/choozdsk01/palomare/SiPM/SC_Fuente_Alpha_Oct/ROOT/run%i_ch4.root",r),"ADC1"},
                      {Form("/pc/choozdsk01/palomare/SiPM/SC_Fuente_Alpha_Oct/ROOT/run%i_ch5.root",r),"ADC2"},
                      {Form("/pc/choozdsk01/palomare/SiPM/SC_Fuente_Alpha_Oct/ROOT/run%i_ch6.root",r),"ADC3"}},
                      adc,range1, range2,250,Nevts);
  myrun.SetGains(Gains);
  myrun.SetSPEAmps(SPEAmp);
  //~myrun.SelectChannels({0,1,2});
  myrun.SelectChannels({0,1,2,3});
  myrun.ParSet->t1=5000e-9;
  myrun.ParSet->t2=2500e-9;
  myrun.ParSet->t3=500e-9;
  
  myrun.ParSet->Peak_ranges={-0.1e-6,0.1e-6,3e-6,4e-6,5e-6};
  myrun.ParSet->Fixed_ranges={1.2e-6,1.4e-6,6.3e-6,7.3e-6,8.3e-6};

  
  
  myrun.Process();
  //  myrun.PlotPeakTimes();
    //~myrun.PlotPedSTDs();
  //  myrun.Plot36("PreTriggerSTD","",0,1);
  //  myrun.LoopWaveforms();
    //~myrun.SetCutPedSTD();
    //~myrun.SetCutVariableVector("PreTriggerSTD",std::map<int,std::pair<double,double>>({{0,{0,myrun.PedestalSTDCUT[0]}},{1,{0,myrun.PedestalSTDCUT[1]}},{2,{0,myrun.PedestalSTDCUT[2]}}}) );
    //~myrun.SetCutTriggerWaveformCuts(triggerchannels);
    //~myrun.SetCutPeakTimeRange(range1,range2);
  //  myrun.Plot36("PreTriggerSTD","",2,1);
  //  myrun.Plot36("PreTriggerCharge","",2,1);

  /*Carga/amplitud media detectada por canal:*/
  double Q0, Q1, Q2, Q3, Q4, A0, A1, A2, A3, A4;
  /*Carga/amplitud media detectada por canal:*/
  Q0 = myrun.getAverage(0,"Q3MaxPeakRange")*1.e12; //en pC
  Q1 = myrun.getAverage(1,"Q3MaxPeakRange")*1.e12;
  Q2 = myrun.getAverage(2,"Q3MaxPeakRange")*1.e12;
  Q3 = myrun.getAverage(3,"Q3MaxPeakRange")*1.e12;
  Q4 = myrun.getAverage(4,"Q3MaxPeakRange")*1.e12;
  A0 = myrun.getAverage(0,"MaxAmplitudeRange"); //en ADC
  A1 = myrun.getAverage(1,"MaxAmplitudeRange"); 
  A2 = myrun.getAverage(1,"MaxAmplitudeRange"); 
  A3 = myrun.getAverage(3,"MaxAmplitudeRange"); 
  A4 = myrun.getAverage(4,"MaxAmplitudeRange");  
  cout << r << " - Duration: " << Form("%.2f",1.0*myrun.getRunDuration()/60) << "min" << " - Evts after cuts" << myrun.getEventsAfterCuts() <<  endl;
  cout << "Charge integrated on 500ns (before signal cuts):" << endl;
  cout << r << " Q - Ch0 " << Q0 << " Ch1 " << Q1 << " Ch2 " << Q2 << " Ch3 " << Q3 << " Ch4 " << Q4 << endl;
  cout << "MaxAmplitudeRange (before signal cuts):" << endl;
  cout << r << " Amp - Ch0 " << A0 << " Ch1 " << A1 << " Ch2 " << A2 << " Ch3 " << A3 << " Ch4 " << A4 << endl;

  /*Signal cuts: Quitamos las waveforms sin señal. Solo afecta a los canales que no son de trigger.*/
  //~myrun.SetCutMaxAmplitudeRange(7,300000);

  /*Carga/amplitud media detectada por canal:*/
  Q0 = myrun.getAverage(0,"Q3MaxPeakRange")*1.e12; //en pC
  Q1 = myrun.getAverage(1,"Q3MaxPeakRange")*1.e12;
  Q2 = myrun.getAverage(2,"Q3MaxPeakRange")*1.e12;
  Q3 = myrun.getAverage(3,"Q3MaxPeakRange")*1.e12;
  Q4 = myrun.getAverage(4,"Q3MaxPeakRange")*1.e12;
  A0 = myrun.getAverage(0,"MaxAmplitudeRange"); //en ADC
  A1 = myrun.getAverage(1,"MaxAmplitudeRange"); 
  A2 = myrun.getAverage(1,"MaxAmplitudeRange"); 
  A3 = myrun.getAverage(3,"MaxAmplitudeRange"); 
  A4 = myrun.getAverage(4,"MaxAmplitudeRange");


  cout << "Charge integrated on 500ns (after signal cuts):" << endl;
  cout << r << " Q - Ch0 " << Q0 << " Ch1 " << Q1 << " Ch2 " << Q2 << " Ch3 " << Q3 << " Ch4 " << Q4 << endl;
  cout << "MaxAmplitudeRange (after signal cuts):" << endl;
  cout << r << " Amp - Ch0 " << A0 << " Ch1 " << A1 << " Ch2 " << A2 << " Ch3 " << A3 << " Ch4 " << A4 << endl;
  ofs << r << "\t" << myrun.GetEventsToProcess()
           << "\t" << Form("%.2f",1.0*myrun.getRunDurationNEvents(myrun.GetEventsToProcess())/60)
           << "\t" << myrun.getEventsAfterCuts() 
           << "\t" << myrun.getWaveformsAfterCuts(0) 
           << "\t" << Q0
           << "\t" << myrun.getWaveformsAfterCuts(1) 
           << "\t" << Q1
           << "\t" << myrun.getWaveformsAfterCuts(2) 
           << "\t" << Q2
           <<  endl;

  ofs2 << r
       << "\t" << myrun.PedestalSTD[0]
       << "\t" << myrun.PedestalSTDCUT[0]
       << "\t" << myrun.PedestalSTD[1]
       << "\t" << myrun.PedestalSTDCUT[1]
       << "\t" << myrun.PedestalSTD[2]
       << "\t" << myrun.PedestalSTDCUT[2]
       << endl;
  //  ofs << "Charge integrated on 500ns (after signal cuts):" << endl;
  //  ofs << r << " Q - Ch0 " << Q0 << " Ch1 " << Q1 << " Ch2 " << Q2 << " Ch3 " << Q3 << " Ch4 " << Q4 << endl;
  //  ofs << "MaxAmplitudeRange (after signal cuts):" << endl;
  //  ofs << r << " Amp - Ch0 " << A0 << " Ch1 " << A1 << " Ch2 " << A2 << " Ch3 " << A3 << " Ch4 " << A4 << endl;

  //  myrun.LoopWaveforms(0,"pqr",NULL,1,0);
  // myrun.Plot36("Charge_Q3MaxPeakRange",Form("AnalysisROOT/Run%i_Q3MaxPeak500ns.root",r));
  myrun.DumpVariableToNtuple(Form("AnalysisROOT/Run%i_NTuple.root",r),0);
  myrun.Close();
}

void SimpleAnalysis()
{

  /*Argumentos:
  1. numero de Run
  2. tiempo inicial donde empezar a buscar la señal de trigger (s)
  3. tiempo final donde dejar de buscar la señal de trigger (s)
  4. Canales de trigger, para limpiar ruido.
  5. Versión del ADC usada (v1720 o DT5725)
  */
  struct runprop
  {
     int run;
     std::vector<int> triggerChannels;
     int Nevts;
  };

  runprop myvector_april[]={

  {37,{1,2},100}

  };
  
  ofstream ofs("AnalysisSimple.log",std::ofstream::out | std::ofstream::app);
  ofstream ofs2("AnalysisSimple_noise.log",std::ofstream::out | std::ofstream::app);
  ofs << "Run\tEvents\tDuration\tEventsAfterCuts\tEventsWithSignal_0\t<Q>_0\tEventsWithSignal_1\t<Q>_1\tEventsWithSignal_2\t<Q>_2" << endl;
  ofs2 << "Run\tPedestalSTD_0\tPedestalSTDCUT_0\tPedestalSTD_1\tPedestalSTDCUT_1\tPedestalSTD_2\tPedestalSTDCUT_2\t" << endl;
  //1.25-1.35
  
  //Normal Runs: 1-39 (5120 samples x 4ns)
  
  for(auto r : myvector_april) Analyse(r.run, 1.2e-6,1.6e-6,r.triggerChannels,ofs,ofs2, r.Nevts);
  
  //~for(auto r : myvector) Analyse(r.run, 1.2e-6,1.6e-6,r.triggerChannels,ofs,ofs2, r.Nevts);
  
  
  //LASER
  //~Analyse(1001, 1.6e-6,1.7e-6,{2},ofs,ofs2, 10000);
  
  ofs.close();
  
}
