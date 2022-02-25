struct variable {string var; string title; double limitdown; double limitup; double UnitConversion;};

std::vector<variable> MyVar(std::vector<variable> myvar, int n, std::vector<string> var, std::vector<int> mych, std::vector<double> Gains, std::vector<double> SPEAmp, double Rebin)
{   
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
    return myvar;
}    
    