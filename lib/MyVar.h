struct variable {string var; string title; double limitdown; double limitup; double UnitConversion;};

std::vector<variable> MyVar(std::vector<variable> myvar, int n, std::vector<string> var, std::vector<int> mych, std::vector<double> conv_fact)
{   
    for(int i=0;i<n;i++)
    {
        if (var[i]=="TStartQpeak") myvar[i]={"TStartQpeak","T (s)",1e-6,3e-6,1.0};
        if (var[i]=="TEndQPeak") myvar[i]={"TEndQPeak","Time (s)",1e-6,3e-6,1.0};
        if (var[i]=="TDiff") myvar[i]={"TDiff","Time (s)",0e-6,2e-6,1.0};

        if (var[i]=="Q1" && (mych[i]==0 || mych[i]==1)) myvar[i]={"Q1","Charge collected (pC)",0,100,1.0/conv_fact[0]}; //SiPM
        if (var[i]=="Q1" && mych[i]==2) myvar[i]={"Q1","Charge collected (pC)",0,5000,1.0/conv_fact[1]}; // PMT
        if (var[i]=="Q1" && mych[i]==3) myvar[i]={"Q1","Charge collected (pC)",0,200,1.0/conv_fact[2]}; //SC

        if (var[i]=="Q2" && (mych[i]==0 || mych[i]==1)) myvar[i]={"Q2","Charge collected (pC)",0,100,1.0/conv_fact[0]}; //SiPM
        if (var[i]=="Q2" && mych[i]==2) myvar[i]={"Q2","Charge collected (pC)",0,5000,1.0/conv_fact[1]}; // PMT
        if (var[i]=="Q2" && mych[i]==3) myvar[i]={"Q2","Charge collected (pC)",0,200,1.0/conv_fact[2]}; //SC

        if (var[i]=="Q3" && (mych[i]==0 || mych[i]==1)) myvar[i]={"Q3","Charge collected (pC)",0,100,1.0/conv_fact[0]}; //SiPM
        if (var[i]=="Q3" && mych[i]==2) myvar[i]={"Q3","Charge collected (pC)",0,5000,1.0/conv_fact[1]}; // PMT
        if (var[i]=="Q3" && mych[i]==3) myvar[i]={"Q3","Charge collected (pC)",0,200,1.0/conv_fact[2]}; //SC
        
        if (var[i]=="QTotal" && (mych[i]==0 || mych[i]==1)) myvar[i]={"QTotal","Charge collected (pC)",0,100,1.0/conv_fact[0]}; //SiPM
        if (var[i]=="QTotal" && mych[i]==2) myvar[i]={"QTotal","Charge collected (pC)",0,5000,1.0/conv_fact[1]}; // PMT
        if (var[i]=="QTotal" && mych[i]==3) myvar[i]={"QTotal","Charge collected (pC)",0,100,1.0/conv_fact[2]}; //SC
        
        if (var[i]=="Amp" && (mych[i]==0 || mych[i]==1)) myvar[i]={"Amp","Amplitude (ADC)",0,1500,1.0}; //SiPM
        if (var[i]=="Amp" && mych[i]==2) myvar[i]={"Amp","Amplitude (ADC)",0,15000,1.0}; // PMT
        if (var[i]=="Amp" && mych[i]==3) myvar[i]={"Amp","Amplitude (ADC)",0,12000,1.0}; //SC

        if (var[i]=="QPeak" && (mych[i]==0 || mych[i]==1)) myvar[i]={"QPeak","Charge Peak",0,200,1.0/(conv_fact[0])}; //SiPM
        if (var[i]=="QPeak" && mych[i]==2) myvar[i]={"QPeak","Charge Peak",0,200,1.0/(conv_fact[1])}; // PMT
        if (var[i]=="QPeak" && mych[i]==3) myvar[i]={"QPeak","Charge Peak",0,200,1.0/(conv_fact[2])}; //SC

        if (var[i]=="QFixRange1" && (mych[i]==0 || mych[i]==1)) myvar[i]={"QFixRange1","Charge collected (pC)",0,120,1.0/(conv_fact[0])}; //SiPM
        if (var[i]=="QFixRange1" && (mych[i]==2)) myvar[i]={"QFixRange1","Charge collected (pC)",0,8000,1.0/(conv_fact[1])}; //PMT
        if (var[i]=="QFixRange1" && (mych[i]==3)) myvar[i]={"QFixRange1","Charge collected (pC)",0,300,1.0/(conv_fact[2])}; //SuperCell

        if (var[i]=="QFixRange2" && (mych[i]==0 || mych[i]==1)) myvar[i]={"QFixRange2","Charge collected (pC)",0,120,1.0/(conv_fact[0])}; //SiPM
        if (var[i]=="QFixRange2" && (mych[i]==2)) myvar[i]={"QFixRange2","Charge collected (pC)",0,8000,1.0/(conv_fact[1])}; //PMT
        if (var[i]=="QFixRange2" && (mych[i]==3)) myvar[i]={"QFixRange2","Charge collected (pC)",0,300,1.0/(conv_fact[2])}; //SuperCell

        if (var[i]=="QFixRange3" && (mych[i]==0 || mych[i]==1)) myvar[i]={"QFixRange3","Charge collected (pC)",0,120,1.0/(conv_fact[0])}; //SiPM
        if (var[i]=="QFixRange3" && (mych[i]==2)) myvar[i]={"QFixRange3","Charge collected (pC)",0,8000,1.0/(conv_fact[1])}; //PMT
        if (var[i]=="QFixRange3" && (mych[i]==3)) myvar[i]={"QFixRange3","Charge collected (pC)",0,300,1.0/(conv_fact[2])}; //SuperCell

        if (var[i]=="QFixRange4" && (mych[i]==0 || mych[i]==1)) myvar[i]={"QFixRange4","Charge collected (pC)",0,120,1.0/(conv_fact[0])}; //SiPM
        if (var[i]=="QFixRange4" && (mych[i]==2)) myvar[i]={"QFixRange4","Charge collected (pC)",0,8000,1.0/(conv_fact[1])}; //PMT
        if (var[i]=="QFixRange4" && (mych[i]==3)) myvar[i]={"QFixRange4","Charge collected (pC)",0,300,1.0/(conv_fact[2])}; //SuperCell

        if (var[i]=="QPeakRange1" && (mych[i]==0 || mych[i]==1)) myvar[i]={"QPeakRange1","Charge collected (pC)",0,400,1.0/(conv_fact[0])}; //SiPM
        if (var[i]=="QPeakRange1" && (mych[i]==2)) myvar[i]={"QPeakRange1","Charge collected (pC)",0,5000,1.0/(conv_fact[1])}; //PMT
        if (var[i]=="QPeakRange1" && (mych[i]==3)) myvar[i]={"QPeakRange1","Charge collected (pC)",0,300,1.0/(conv_fact[2])}; //SuperCell

        if (var[i]=="QPeakRange2" && (mych[i]==0 || mych[i]==1)) myvar[i]={"QPeakRange2","Charge collected (pC)",0,400,1.0/(conv_fact[0])}; //SiPM
        if (var[i]=="QPeakRange2" && (mych[i]==2)) myvar[i]={"QPeakRange2","Charge collected (pC)",0,5000,1.0/(conv_fact[1])}; //PMT
        if (var[i]=="QPeakRange2" && (mych[i]==3)) myvar[i]={"QPeakRange2","Charge collected (pC)",0,300,1.0/(conv_fact[2])}; //SuperCell

        if (var[i]=="QPeakRange3" && (mych[i]==0 || mych[i]==1)) myvar[i]={"QPeakRange3","Charge collected (pC)",0,400,1.0/(conv_fact[0])}; //SiPM
        if (var[i]=="QPeakRange3" && (mych[i]==2)) myvar[i]={"QPeakRange3","Charge collected (pC)",0,5000,1.0/(conv_fact[1])}; //PMT
        if (var[i]=="QPeakRange3" && (mych[i]==3)) myvar[i]={"QPeakRange3","Charge collected (pC)",0,300,1.0/(conv_fact[2])}; //SuperCell

        if (var[i]=="QPeakRange4" && (mych[i]==0 || mych[i]==1)) myvar[i]={"QPeakRange4","Charge collected (pC)",0,400,1.0/(conv_fact[0])}; //SiPM
        if (var[i]=="QPeakRange4" && (mych[i]==2)) myvar[i]={"QPeakRange4","Charge collected (pC)",0,5000,1.0/(conv_fact[1])}; //PMT
        if (var[i]=="QPeakRange4" && (mych[i]==3)) myvar[i]={"QPeakRange4","Charge collected (pC)",0,300,1.0/(conv_fact[2])}; //SuperCell
    }
    return myvar;
}    
    