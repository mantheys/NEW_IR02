#include "lib/headers.h"

void Analyse(string path, int r, int ch, int ped, double range1, double range2, std::vector<bool> conditions) 
{// Macro para visualizar eventos y ver cómo afectan los cortes que queremos establecer
  ana::Run_t myrun(r,{{path+Form("run%i_ch%i.root",r,ch),"ADC0"}},"DT5725",range1, range2, ped, -1);
  
  std::vector<double> SPEAmp={38.6,24.8,25.5};//Amplitud del SPE en cuentas de ADC
  myrun.SetSPEAmps(SPEAmp);
 
  myrun.SelectChannels({0});
  myrun.Process();
  myrun.ParSet->t3 = 500e-9; //Fijamos el rango de integración de Q3 como 500ns tras el pico.
	
  if (conditions[0] == true){myrun.PlotPedestals();}
  if (conditions[1] == true){myrun.PlotPeakTimes();}
  if (conditions[2] == true)
  {
    TH1F* h0 = myrun.TH1Charge(0,"Range","pC");h0->Draw();gPad->Update();lets_pause();
    if (conditions[3] == true){myrun.autofit(h0);lets_pause();}
  }
  if (conditions[4] == true){myrun.LoopWaveforms(0,"paqr",NULL);}

  myrun.Close();
}

void EventDisplay(string input = "config_file,txt")
{
  int run; int ch; int ped;
  run = IntInput(input, "RUN"); 
  ch = IntInput(input, "CHANNEL");
  ped = IntInput(input, "PEDESTALRANGE");

  double isignaltime; double fsignaltime;
  isignaltime = DoubleInput(input, "ISIGNALTIME");
  fsignaltime = DoubleInput(input, "FSIGNALTIME");

  string path;
  path = StringInput(input, "PATH");

  std::vector<string> keywords; std::vector<bool> conditions; conditions = {};
  keywords = {"PLOTPEDESTALS","PLOTPEAKTIMES","CHARGEHIST", "CHARGEHISTAUTOFIT", "EVENTDISPLAY"};

  for(vector<string>::const_iterator key = keywords.begin(); key != keywords.end(); ++key)
  {bool condition; condition = BoolInput(input, *key); conditions.push_back(condition);}
  
  Analyse(path, run, ch, ped, isignaltime, fsignaltime, conditions);
  /*  0.  Path de la carpeta que incluye los archivos .root
      1.  Numero de Run
      2.  Canal del ADC que figura en el nombre del .root
      3.  Rango de tiempo que caracteriza el pedestal (en unidades de 4us)
      3.  Tiempo inicial donde empezar a buscar la señal de trigger (s)
      4.  Tiempo final donde dejar de buscar la señal de trigger (s)
      5.  Path de la carpeta que incluye los archivos .root
      6.  Vector de condiciones (0/1 false/true) para activar las funciones en orden de pariencia.
  */
}
