#include "lib/headers.h"

void Analyse(string adc, string path, string output, int r, int ch, int ped, double range1, double range2, int min_amp, int max_amp, double conv_factor, std::vector<bool> conditions)
{ /* Macro para obtener la waveform promedio de un run.
    La macro lee los ficheros con los runes guardados en ROOT, y crea un fichero en AnalysisROOT con los perfiles de centelleo para todos los canales.
  */

  ana::Run_t myrun(r,{{path+Form("run%02i_ch%i.root",r,ch),"ADC0"}}, adc, range1, range2, ped, -1, conv_factor);

  // Variables requeridas para todas las funciones que s emuestran a continuación
  std::vector<double> SPEAmp = {38.6, 24.8, 25.5}; myrun.SetSPEAmps(SPEAmp);
  myrun.SelectChannels({0});

  int show = 0;

  if (conditions[0] == true){myrun.SetCutMaxAmplitudeRange(min_amp,max_amp);}
  if (conditions[1] == true){myrun.SetCutPedSTD();}
  if (conditions[2] == true){myrun.SetMaximumWaveformsToProcess(-1);}
  if (conditions[3] == true){myrun.ParSet->setADCAmplitudeThreshold(-1000);}
  if (conditions[4] == true){show = 1;}
  myrun.Plot36("ScintProfFirstSignalBin", output+Form("_run%02i_ch%i.root",r,ch), 0, show);
  
  myrun.Close();
}

void AverageWaveform(string input = "CONFIG/aw_config_file.txt")
{
  int irun; int frun; int ch; int ped; int min_amp; int max_amp;
  irun = IntInput(input, "I_RUN"); frun = IntInput(input, "F_RUN"); ch = IntInput(input, "CHANNEL"); ped = IntInput(input, "PEDESTAL_RANGE");
  min_amp = IntInput(input, "MIN_AMP"); max_amp = IntInput(input, "MAX_AMP");
  
  double isignaltime; double fsignaltime; double conv_factor;
  isignaltime = DoubleInput(input, "I_SIGNALTIME"); fsignaltime = DoubleInput(input, "F_SIGNALTIME"); conv_factor = DoubleInput(input, "CONVERSION_FACTOR");

  string adc; string path; string output; 
  adc = StringInput(input, "ADCMODE"); path = StringInput(input, "PATH"); output = StringInput(input, "OUTPUT_FILE");

  std::vector<string> keywords; std::vector<bool> conditions; conditions = {};
  keywords = {"CUT_MAX_AMPLITUDE","CUT_PED_STD","AMP_WAVEFORMS","ADC_AMP_THRESHOLD","SHOW_AVERAGE"};

  for(vector<string>::const_iterator key = keywords.begin(); key != keywords.end(); ++key)
  {bool condition; condition = BoolInput(input, *key); conditions.push_back(condition);}

  for (int run=irun; run<=frun; run++){Analyse(adc, path, output, run, ch, ped, isignaltime, fsignaltime, min_amp, max_amp, conv_factor, conditions);}
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
