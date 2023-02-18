%module csound 
%{
#include "csound.hpp"
%}

%inline %{
typedef MYFLT cfloat;
%}

class Csound
{
public:
  virtual CSOUND *GetCsound();
  virtual void SetCsound(CSOUND *csound_);
  virtual int InitializeCscore(FILE *insco, FILE *outsco);
  virtual int LoadPlugins(const char *dir);
  virtual int GetVersion();
  virtual int GetAPIVersion();
  virtual void *GetHostData();
  virtual void SetHostData(void *hostData);
  virtual const char *GetEnv(const char *name);
  virtual int SetGlobalEnv(const char *name,const char *value);
  virtual int SetOption(const char *option);
  virtual void SetParams(CSOUND_PARAMS *p);
  virtual void GetParams(CSOUND_PARAMS *p);
  virtual void SetOutput(const char *name,const char *type,const char *format);
  virtual void SetInput(const char *name);
  virtual void SetMIDIInput(const char *name);
  virtual void SetMIDIFileInput(const char *name);
  virtual void SetMIDIOutput(const char *name);
  virtual void SetMIDIFileOutput(const char *name);
  virtual TREE *ParseOrc(const char *str);
  virtual int CompileTree(TREE *root);
  virtual void DeleteTree(TREE *root);
  virtual int CompileOrc(const char *str);
  virtual cfloat EvalCode(const char *str);
  virtual int ReadScore(const char *str);
  virtual int CompileArgs(int argc,const char **argv);
  virtual int Compile(int argc,const char **argv);
  virtual int Compile(const char *csdName);
  virtual int Compile(const char *orcName,const  char *scoName);
  virtual int Compile(const char *arg1,const  char *arg2,const  char *arg3);
  virtual int Compile(const char *arg1,const  char *arg2,const  char *arg3,const  char *arg4);
  virtual int Compile(const char *arg1,const  char *arg2,const char *arg3, const char *arg4,const  char *arg5);
  virtual int CompileCsd(const char *csd);
  virtual int CompileCsdText(const char *csd_text);
  virtual int Start();
  virtual int Perform();
  virtual int Perform(int argc, const char **argv);
  virtual int Perform(const char *csdName);
  virtual int Perform(const char *orcName, const char *scoName);
  virtual int Perform(const char *arg1, const char *arg2, const char *arg3);
  virtual int Perform(const char *arg1, const  char *arg2, const  char *arg3, const  char *arg4);
  virtual int Perform(const char *arg1, const  char *arg2, const  char *arg3,const char *arg4, const  char *arg5);
  virtual int PerformKsmps();
  virtual int PerformBuffer();
  virtual void Stop();
  virtual int Cleanup();
  virtual void Reset();
  virtual cfloat GetSr();
  virtual cfloat GetKr();
  virtual int GetKsmps();
  virtual int GetNchnls();
  virtual int GetNchnlsInput();
  virtual cfloat Get0dBFS();
  virtual long GetInputBufferSize();
  virtual long GetOutputBufferSize();
  virtual cfloat *GetInputBuffer();
  virtual cfloat *GetOutputBuffer();
  virtual cfloat *GetSpin();
  virtual cfloat *GetSpout();
  virtual const char *GetOutputName();
  virtual long GetCurrentTimeSamples();
  virtual void SetHostImplementedAudioIO(int state, int bufSize);
  virtual void SetHostImplementedMIDIIO(int state);
  virtual double GetScoreTime();
  virtual int IsScorePending();
  virtual void SetScorePending(int pending);
  virtual cfloat GetScoreOffsetSeconds();
  virtual void SetScoreOffsetSeconds(double time);
  virtual void RewindScore();
  virtual void SetCscoreCallback(void (*cscoreCallback_)(CSOUND *));
  virtual int ScoreSort(FILE *inFile, FILE *outFile);
  virtual int ScoreExtract(FILE *inFile, FILE *outFile, FILE *extractFile);
  
  //virtual void Message(const char *format, ...)
  //virtual void MessageS(int attr, const char *format, ...)
  //virtual void MessageV(int attr, const char *format, va_list args)

  /*
  virtual void SetMessageCallback(
      void (*csoundMessageCallback_)(CSOUND *, int attr,
                                     const char *format, va_list valist))
  */                                
  virtual int GetMessageLevel();
  virtual void SetMessageLevel(int messageLevel);
  virtual void InputMessage(const char *message);
  virtual void KeyPressed(char c);
  virtual int ScoreEvent(char type, const cfloat *pFields, long numFields);
  virtual int ScoreEventAbsolute(char type, const cfloat *pFields,long numFields, double time_ofs);
  
  //virtual void SetExternalMidiInOpenCallback(int (*func)(CSOUND *, void **, const char *));
  //virtual void SetExternalMidiReadCallback(int (*func)(CSOUND *, void *, unsigned char *, int));
  //virtual void SetExternalMidiInCloseCallback(int (*func)(CSOUND *, void *));
  //virtual void SetExternalMidiOutOpenCallback(int (*func)(CSOUND *, void **, const char *));
  //virtual void SetExternalMidiWriteCallback(int (*func)(CSOUND *, void *, const unsigned char *, int));
  //virtual void SetExternalMidiOutCloseCallback(int (*func)(CSOUND *, void *));
  //virtual void SetExternalMidiErrorStringCallback(const char *(*func)(int));
  
  virtual int SetIsGraphable(int isGraphable);
  
  //virtual void SetMakeGraphCallback(void (*makeGraphCallback_)(CSOUND *, WINDAT *windat, const char *name));
  //virtual void SetDrawGraphCallback(void (*drawGraphCallback_)(CSOUND *, WINDAT *windat));
  //virtual void SetKillGraphCallback(void (*killGraphCallback_)(CSOUND *, WINDAT *windat));
  //virtual void SetExitGraphCallback(int (*exitGraphCallback_)(CSOUND *));
  
  virtual int NewOpcodeList(opcodeListEntry* &opcodelist);
  virtual void DisposeOpcodeList(opcodeListEntry *opcodelist);
  /*
  virtual int AppendOpcode(const char *opname, int dsblksiz, int flags,
                           int thread,
                           const char *outypes, const char *intypes,
                           int (*iopadr)(CSOUND *, void *),
                           int (*kopadr)(CSOUND *, void *),
                           int (*aopadr)(CSOUND *, void *))
*/

  virtual void SetYieldCallback(int (*yieldCallback_)(CSOUND *));
  virtual void SetPlayopenCallback(int (*playopen__)(CSOUND *, const csRtAudioParams *parm));
  virtual void SetRtplayCallback(void (*rtplay__)(CSOUND *, const cfloat *outBuf, int nbytes));
  virtual void SetRecopenCallback(int (*recopen_)(CSOUND *, const csRtAudioParams *parm));
  virtual void SetRtrecordCallback(int (*rtrecord__)(CSOUND *, cfloat *inBuf, int nbytes));
  virtual void SetRtcloseCallback(void (*rtclose__)(CSOUND *));
  virtual int GetDebug();
  virtual void SetDebug(int debug);
  virtual int TableLength(int table);
  virtual cfloat TableGet(int table, int index);
  virtual void TableSet(int table, int index, double value);
  virtual int GetTable(cfloat* &tablePtr, int tableNum);
  virtual void TableCopyOut(int table, cfloat *dest);
  virtual void TableCopyIn(int table, cfloat *src);
  virtual int CreateGlobalVariable(const char *name, size_t nbytes);
  virtual void *QueryGlobalVariable(const char *name);
  virtual void *QueryGlobalVariableNoCheck(const char *name);
  virtual int DestroyGlobalVariable(const char *name);
  virtual void **GetRtRecordUserData();
  virtual void **GetRtPlayUserData();
  virtual int RegisterSenseEventCallback(void (*func)(CSOUND *, void *),void *userData);
  virtual int RunUtility(const char *name, int argc, char **argv);
  virtual char **ListUtilities();
  virtual void DeleteUtilityList(char **lst);
  virtual const char *GetUtilityDescription(const char *utilName);
  virtual int GetChannelPtr(cfloat* &p, const char *name, int type);
  virtual int ListChannels(controlChannelInfo_t* &lst);
  virtual void DeleteChannelList(controlChannelInfo_t *lst);
  virtual int SetControlChannelHints(const char *name, controlChannelHints_t hints);
  virtual int GetControlChannelHints(const char *name, controlChannelHints_t *hints);
  virtual void SetChannel(const char *name, double value);
  virtual void SetControlChannel(const char *name, double value);
  virtual void SetChannel(const char *name, char *string);
  virtual void SetStringChannel(const char *name, char *string);
  virtual void SetChannel(const char *name, cfloat *samples);
  virtual cfloat GetChannel(const char *name, int *err = NULL);
  virtual cfloat GetControlChannel(const char *name, int *err = NULL);
  virtual void GetStringChannel(const char *name, char *string);
  virtual void GetAudioChannel(const char *name, cfloat *samples);
  virtual int PvsinSet(const PVSDATEXT* value, const char *name);
  virtual int PvsoutGet(PVSDATEXT* value, const char *name);
  virtual void SetInputChannelCallback(channelCallback_t inputChannelCalback);
  virtual void SetOutputChannelCallback(channelCallback_t outputChannelCalback);
  
  virtual int CreateConfigurationVariable(const char *name, void *p,
                                          int type, int flags,
                                          void *min, void *max,
                                          const char *shortDesc,
                                          const char *longDesc);
  virtual int SetConfigurationVariable(const char *name, void *value);
  virtual int ParseConfigurationVariable(const char *name, const char *value);
  virtual csCfgVariable_t *QueryConfigurationVariable(const char *name);
  virtual csCfgVariable_t **ListConfigurationVariables();
  virtual int DeleteConfigurationVariable(const char *name);
  
  Csound();
  Csound(CSOUND *csound_);
  Csound(void *hostData);
  virtual ~Csound();
  virtual void CreateMessageBuffer(int toStdOut);
  virtual const char *GetFirstMessage();
  virtual int GetFirstMessageAttr();
  virtual void PopFirstMessage();
  virtual int GetMessageCnt();
  virtual void DestroyMessageBuffer();
  virtual void AddSpinSample(int frame, int channel, cfloat sample);
  virtual void SetSpinSample(int frame, int channel, cfloat sample);
  virtual cfloat GetSpoutSample(int frame, int channel) const;
  virtual const char *GetInputName();
  virtual void SetAudioChannel(const char *name, cfloat *samples);
  virtual cfloat SystemSr(cfloat value);
};

class CsoundThreadLock {
public:
  int Lock(size_t milliseconds);
  void Lock();
  int TryLock();
  void Unlock();

  CsoundThreadLock();
  CsoundThreadLock(int locked);
  ~CsoundThreadLock();
};

class CsoundMutex {
public:
  void Lock();
  int TryLock();
  void Unlock();

  CsoundMutex();

  explicit CsoundMutex(int isRecursive);
  ~CsoundMutex();
};

class CsoundRandMT 
{
public:
  uint32_t Random();
  void Seed(uint32_t seedVal);
  void Seed(const uint32_t *initKey, int keyLength);

  CsoundRandMT();

  explicit CsoundRandMT(uint32_t seedVal);
  CsoundRandMT(const uint32_t *initKey, int keyLength);
  ~CsoundRandMT();
};

// timer (csoundInitialize() should be called before using this)

class CsoundTimer {
public:
  double GetRealTime();
  double GetCPUTime();
  void Reset();

  CsoundTimer();
  ~CsoundTimer();
};
