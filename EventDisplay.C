#include "lib/headers.h"
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//___EVENTDISPLAY.C ES UNA MACRO PARA VISUALIZAR EVENTOS Y HACER UN ANALISIS PREELIMINAR E INDIVIDUALIZADO DE CORTES Y FACTORES DETERMINANTES___//
//___EXECUTE USING THE FOLLOWING COMMAND: root -l EventDisplay.C+(\"config_file.txt\") _________________________________________________________//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Analyse(string adc, string path, int r, int ch, int ped, double range1, double range2, std::vector<bool> conditions) 
{ /* Macro para visualizar eventos y ver cómo afectan los cortes que queremos establecer
  En Analyse se incluyen las variables que se pasan a la clase Run_t y las condiciones de activación del resto de funciones. */
  
  ana::Run_t myrun(r,{{path+Form("run%02i_ch%i.root",r,ch),"ADC0"}}, adc, range1, range2, ped, -1);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //___ESTE ES EL CUERPO PRINCIPAL DE LA MACRO DONDE SE CONFIGURA SU FUNCIONALIDAD Y FACTORES DE CONVERSIÓN___//
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Variables requeridas para todas las funciones que s emuestran a continuación
  std::vector<double> SPEAmp={38.6,24.8,25.5}; myrun.SetSPEAmps(SPEAmp);//Amplitud del SPE en cuentas de ADC
  myrun.SelectChannels({0}); myrun.Process();
  myrun.ParSet->t3 = 500e-9; //Fijamos el rango de integración de Q3 como 500ns tras el pico.
	
  // Funciones que se aplicna a la clase Run_t y tienen como finalidad visualizar eventos individuales o hacer un estudio preliminar
  if (conditions[0] == true){myrun.PlotPedestals();}
  if (conditions[1] == true){myrun.PlotPeakTimes();}
  if (conditions[2] == true)
  {
    TH1F *h0 = myrun.TH1Charge(0,"Range","pC");h0->Draw();gPad->Update();lets_pause();
    if (conditions[3] == true){myrun.autofit(h0, false, "test");}
  }
  if (conditions[4] == true){TH1F *h0 = myrun.TH1Amp(0);h0->Draw();gPad->Update();lets_pause();}
  if (conditions[5] == true){myrun.LoopWaveforms(0,"paqr",NULL);}
 
  myrun.Close();
}

void EventDisplay(string input = "config_file.txt")
{ /* Función principal de esta macro. Las direfentes funciones _Input() llaman al archivo de configuracion e importan las variables pertinentes */

  /////////////////////////////////////////////////////////////////////
  //___AQUI SE IMPORTAN LAS VARIABLES DEL ARCHIVO DE CONFIGURACIÓN___//
  /////////////////////////////////////////////////////////////////////  
  
  int irun; int frun; int ch; int ped;
  irun = IntInput(input, "I_RUN"); frun = IntInput(input, "F_RUN"); ch = IntInput(input, "CHANNEL"); ped = IntInput(input, "PEDESTAL_RANGE");

  double isignaltime; double fsignaltime;
  isignaltime = DoubleInput(input, "I_SIGNALTIME"); fsignaltime = DoubleInput(input, "F_SIGNALTIME");

  string adc; string path;
  adc = StringInput(input, "ADC"); path = StringInput(input, "PATH");

  std::vector<string> keywords; std::vector<bool> conditions; conditions = {};
  keywords = {"PLOT_PEDESTALS","PLOT_PEAKTIMES","CHARGE_HIST","CHARGE_HIST_AUTOFIT","MAX_AMP_HIST","EVENT_DISPLAY"};

  for(vector<string>::const_iterator key = keywords.begin(); key != keywords.end(); ++key)
  {bool condition; condition = BoolInput(input, *key); conditions.push_back(condition);}

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //___LAS VARIABLES QUE SE HAN IMPORTADO SE PASAN A LA FUNCIÓN ANALYSE QUE A SU VEZ LLAMA AL RUN_T PERTINENTE Y EJECUTA LAS FUNCIONES ESCOGIDAS___//  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for (int run=irun; run<=frun; run++) Analyse(adc, path, run, ch, ped, isignaltime, fsignaltime, conditions);
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
