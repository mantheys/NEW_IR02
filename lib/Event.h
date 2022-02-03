namespace ana
{

  class Event_t
  {
  public:
    std::vector<waveana::Waveform_t> waveform;
    bool cut = false;
    double TriggerTime = 0;
    double TimeStamp = 0;
    std::map<int, int> ChannelToVectorPosition;

  public:
    Event_t() {}
    void resize(int nevents)
    {
      waveform.resize(nevents);
    }
    void AddWaveform(waveana::Waveform_t w, int channel)
    {
      waveform.push_back(w);
      ChannelToVectorPosition[channel] = waveform.size() - 1;
    }
    waveana::Waveform_t *getChannel(int pm) { return &(waveform[ChannelToVectorPosition[pm]]); }
    void setCut(bool c) { cut = c; }
    bool getCut() { return cut; }
    double getTriggerTime() { return TriggerTime; }
    void setTriggerTime(double c) { TriggerTime = c; }
    double getTimeStamp() { return TimeStamp; }
    void setTimeStamp(double c) { TimeStamp = c; }
    void Clear()
    {
      waveform.clear();
      ChannelToVectorPosition.clear();
    }
  };

}
