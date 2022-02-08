#include "lib/headers.h"

void Analyse(int r, double range1, double range2, std::vector<int> triggerchannels, int NEvts = -1, string adc = "DT5725")
{ /* Macro para obtener la waveform promedio de un run.
    La macro lee los ficheros con los runes guardados en ROOT, y crea un fichero en AnalysisROOT con los perfiles de centelleo para todos los canales.
  */
  std::vector<double> Gains = {0.56, 3.61, 3.83};  // ganancias en pC
  std::vector<double> SPEAmp = {38.6, 24.8, 25.5}; // Amplitud del SPE en cuentas de ADC
  for (int i = 0; i < 3; i++)
  Gains[i] = Gains[i] / 1.602e-7; //(1e-19+12)
  
  ana::Run_t myrun(r, {{Form("/pc/choozdsk01/palomare/SiPM/SC_Fuente_Alpha_Jan/ROOT/run%i_ch6.root", r),"ADC2"}}, adc, range1, range2, 425, -1);
  myrun.SetGains(Gains);
  myrun.SetSPEAmps(SPEAmp);
  myrun.SelectChannels({0});

  myrun.SetCutMaxAmplitudeRange(0,20);
  myrun.SetMaximumWaveformsToProcess(NEvts);
  myrun.ParSet->setADCAmplitudeThreshold(-1000);

  myrun.Plot36("ScintProfFirstSignalBin", Form("/pc/choozdsk01/palomare/SiPM/SC_Fuente_Alpha_Jan/AnalysisROOT/run%i_ScintProfFirstSignalBin_SC_Noise.root", r), 0, 1);
  //myrun.Plot36("ScintProf", Form("AnalysisROOT/run%i_ScintProfNoise.root", r), 0, 1);

  myrun.Close();
}

void AverageWaveform()
{
  Analyse(81, 0.0, 4e-6, {0}, 1);
  /*Argumentos:
  1. numero de Run
  2. tiempo inicial donde empezar a buscar la señal de trigger (s)
  3. tiempo final donde dejar de buscar la señal de trigger (s)
  4. Canales de trigger, para limpiar ruido.
  5. Versión del ADC usada (v1720 o DT5725)
  */
}