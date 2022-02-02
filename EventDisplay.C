/*Macro para visualizar eventos y como afectan los cortes que queremos establecer a la selección de eventos*/

#include"lib/headers.h"

void Analyse(int r, int ch, double range1, double range2, std::vector<int> triggerchannels, string path = "ROOT/", string adc = "DT5725", string m = "Month", bool pedestal_plot = false, bool peaktimes_plot = false, bool chargehist_plot = false, bool event_plot = false)
{
  ana::Run_t myrun(r,{{path+Form("run%i_ch%i.root",r,ch),"ADC2"}},adc,range1, range2, 200, -1);
  
  std::vector<double> SPEAmp={38.6,24.8,25.5};//Amplitud del SPE en cuentas de ADC
  myrun.SetSPEAmps(SPEAmp);
 
  myrun.SelectChannels({0});
  myrun.Process();
  myrun.ParSet->t3 = 500e-9; //Fijamos el rango de integración de Q3 como 500ns tras el pico.
	
  if (pedestal_plot == true){myrun.PlotPedestals();}
  if (peaktimes_plot == true){myrun.PlotPeakTimes();}
  if (chargehist_plot == true){
    TH1F* h0 = myrun.TH1Charge(0,"Range","pC");
    h0->Draw();
    lets_pause();
  }
  // string BoxOption="", std::vector<double> *ranges=NULL, int Rebin=1, int CutOption=-1, //este loop nos muestra las waveforms evento a evento, marcando cuales nos estaríamos quitando con nuestros cortes.
  if (event_plot == true){myrun.LoopWaveforms(0,"paqr",NULL);}

  myrun.Close();
}

void EventDisplay(string input = "config_file,txt")
{
  int run; int ch; 
  run = IntInput(input, "RUN"); 
  ch = IntInput(input, "CHANNEL");

  double isignaltime; double fsignaltime;
  isignaltime = DoubleInput(input, "ISIGNALTIME");
  fsignaltime = DoubleInput(input, "FSIGNALTIME");

  string path;
  path = StringInput(input, "PATH");

  bool pedestal; bool peaktimes; bool chargehist; bool event;
  pedestal = BoolInput(input, "PLOTPEDESTALS");
  peaktimes = BoolInput(input, "PLOTPEAKTIMES");
  chargehist = BoolInput(input, "CHARGEHIST");
  event = BoolInput(input, "EVENTDISPLAY");
  

  Analyse(run, ch, isignaltime, fsignaltime, {0}, path, "DT5725", "Month", pedestal, peaktimes, chargehist, event);
  /*Argumentos:
  1.  Numero de Run
  2.  Tiempo inicial donde empezar a buscar la señal de trigger (s)
  3.  Tiempo final donde dejar de buscar la señal de trigger (s)
  4.  Vector con los canales de trigger, para limpiar ruido.
  5.  Path de la carpeta que incluye los archivos .root
  5.  Versión del ADC usada (v1720 o DT5725)
  6.  Nombre descriptivo
  7.  Condición para aplicar la función PlotPedestals
  8.  Condición para aplicar la función PlotPeakTimes
  9.  Condición para aplicar la función TH1Charge
  10. Condición para aplicar la función LoopWaveforms
  */
}
