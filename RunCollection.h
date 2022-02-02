class AnaRun_t
{
public:
  int run;
  int trigger;
  string threshold;
  TDatime date;
  AnaRun_t(int r, int t, string th, TDatime d) : run(r), trigger(t), threshold(th), date(d) {}
};
class RunCollection_t
{
  std::vector<AnaRun_t> list;

public:
  RunCollection_t(string file = "Runs.csv")
  {
    ifstream ifs(file);
    int r, t;
    string th, b, d;
    getline(ifs, th);
    cout << th << endl;
    while (1)
    {
      ifs >> r >> t >> th >> b >> d; // cout << r << " " << t << " " << th << " " << b << " " << d << endl;
      int dia = stoi(d.substr(0, 2));
      int mes = stoi(d.substr(3, 2));
      int anyo = stoi(d.substr(6, 4));
      if (ifs.eof())
        break;
      AnaRun_t ar(r, ConvertTrigger(t), th, TDatime(anyo, mes, dia, 0, 0, 0));
      list.push_back(ar);
    }
  }

  int ConvertTrigger(int t)
  {
    if (t == 1)
      return 0;
    else if (t == 2)
      return 1;
    else if (t == 3)
      return 2;
    else if (t == 12)
      return 10;
    else
      throw std::exception();
  }
  string TriggerToString(int t)
  {
    if (t == 0)
      return "(0)";
    else if (t == 1)
      return "(1)";
    else if (t == 2)
      return "(2)";
    else if (t == 10)
      return "(0,1)";
    else
      throw std::exception();
  }

  std::vector<std::pair<int, TDatime>> GetByTriggerThreshold(int tr, string th)
  {
    std::vector<std::pair<int, TDatime>> v;
    for (auto l : list)
    {
      if (l.trigger == tr && l.threshold == th)
        v.push_back(std::pair<int, TDatime>({l.run, l.date}));
    }
    return v;
  }
  std::vector<std::pair<int, string>> GetByTriggerDate(int tr, TDatime d)
  {
    std::vector<std::pair<int, string>> v; // returns run,threshold
    for (auto l : list)
    {
      if (l.trigger == tr && l.date == d)
        v.push_back(std::pair<int, string>({l.run, l.threshold}));
    }
    return v;
  }
  std::vector<std::pair<int, int>> GetByThresholdDate(string th, TDatime d)
  {
    std::vector<std::pair<int, int>> v;
    for (auto l : list)
    {
      if (l.threshold == th && l.date == d)
        v.push_back(std::pair<int, int>({l.run, l.trigger}));
    }
    return v;
  }
};
