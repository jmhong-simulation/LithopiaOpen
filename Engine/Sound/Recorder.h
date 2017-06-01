// Recorder.h: interface for the CRecorder class.
//
//////////////////////////////////////////////////////////////////////

#pragma once

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include <windows.h>
#include <mmsystem.h>
//#pragma comment(lib,"winmm.lib")
#include <stdio.h>

class BEAT_DETECTOR
{
public:
	int fft_length;
	int subband_length;
	int number_of_temporal_samples;
	int number_of_recent_samples;

	float one_over_number_of_recent_samples;
	float one_over_number_of_temporal_samples;

public:
	float **temporal_subbands;
	float *current_subbands;
	float *subband_beats;
	float *subbands;// copy the values current_subbands to communicate with outside variables
	float *normalized_subbands;//normalized values of subbands by mean values

public:
	BEAT_DETECTOR(){};//cannot be used before being initialized

	BEAT_DETECTOR(const int fft_length_input,const int subband_length_input,const int number_of_temporal_samples_input,const int number_of_recent_samples_input)
	{Initialize(fft_length_input,subband_length_input,number_of_temporal_samples_input,number_of_recent_samples_input);}

	~BEAT_DETECTOR(){
//		for(int i=0;i<number_of_temporal_samples;i++)delete [] temporal_subbands[i];
		delete [] temporal_subbands;
		delete [] subband_beats;
		delete [] subbands;
		delete [] normalized_subbands;
	}

public:
	void Initialize(const int fft_length_input,const int subband_length_input,const int number_of_temporal_samples_input,const int number_of_recent_samples_input)
	{
		fft_length=fft_length_input;
		subband_length=subband_length_input;
		number_of_temporal_samples=number_of_temporal_samples_input;
		number_of_recent_samples=number_of_recent_samples_input;
		
		one_over_number_of_recent_samples=1.0f/(float)number_of_recent_samples;
		one_over_number_of_temporal_samples=1.0f/(float)number_of_temporal_samples;

		temporal_subbands=new float*[number_of_temporal_samples];
		for(int i=0;i<number_of_temporal_samples;i++)temporal_subbands[i]=new float [subband_length];
		for(int i=0;i<number_of_temporal_samples;i++)for(int j=0;j<subband_length;j++)temporal_subbands[i][j]=0.0f;

		current_subbands=temporal_subbands[0];
		subband_beats=new float [subband_length];
		
		subbands=new float [subband_length];
		for (int i = 0; i < subband_length; i++) subbands[i] = 0.0f;

		normalized_subbands=new float [subband_length];
	}

	void Update(const float* spectrum)// length of spectrum SHOULD be fft_length
	{
	/*	static DWORD start_time=timeGetTime();
		static int time_count=0;
		time_count++;
		if(time_count>100){
			printf("%f ",100.0f/(float)(timeGetTime()-start_time));
			start_time=timeGetTime();
			time_count=0;}*/
		delete temporal_subbands[number_of_temporal_samples-1];
		for(int i=number_of_temporal_samples-1;i>=1;i--)temporal_subbands[i]=temporal_subbands[i-1];
		current_subbands=new float [subband_length];
		for(int i=0;i<subband_length;i++)current_subbands[i]=0.0f;
		temporal_subbands[0]=current_subbands;
		
		static int a=fft_length/(subband_length*(subband_length+1)/2);
		int count=1;	
		for(int i=0;i<subband_length;i++){//subband_length error
			float sum=0.0f;	
			for(int j=1;j<=i*a;j++){
				sum+=spectrum[count];
				count++;}
			current_subbands[i]=sum;}

		//Detect beats
		for(int j=0;j<subband_length;j++){
			float mean=0.0f,recent_mean=0.0f;
			for(int i=0;i<number_of_temporal_samples;i++){
				mean+=temporal_subbands[i][j];
				if(i<number_of_recent_samples)recent_mean+=temporal_subbands[i][j];}
			mean*=one_over_number_of_temporal_samples;
			recent_mean*=one_over_number_of_recent_samples;
			
			if(mean<=0.0f)normalized_subbands[j]=0.0f;//filter white noise
			else normalized_subbands[j]=current_subbands[j]/mean;
			
			if(recent_mean>mean*1.5f)subband_beats[j]=current_subbands[j];
			else subband_beats[j]=0.0f;}

		for(int i=0;i<subband_length;i++)subbands[i]=current_subbands[i];
	}
};

#define FFT_LEN 2048
#define MAXNUMOFBUFFER 16

typedef struct _PCMFORMAT{
	WORD    wBitsPerSample;//no.of bits per sample for Recording 
	WORD	wChannels;//no.of channels for Recording
	DWORD	dwSampleRate;//Sampling rate for Recording
}PCMFORMAT, *LPPCMFORMAT;

struct WAVECLASS{HWAVE hWave;void* lpData;};
struct DataHolder{void* pData;void* pData2;};

typedef BOOL (*ProcessBuffer)(void* lpData, LPWAVEHDR pwh);

class Recorder
{
public:
	Recorder(int nBufferLength=FFT_LEN);
	Recorder(PCMFORMAT pcm,int nBufferLength=FFT_LEN);
	Recorder(WORD wBitsPerSample,WORD wChannels,DWORD dwSampleRate,int nBufferLength=FFT_LEN);
	virtual ~Recorder();

public:
	void Open(DWORD dwCallBack=NULL,DWORD dwCallbackType=CALLBACK_FUNCTION,MCIDEVICEID wMCIDeviceID=WAVE_MAPPER);
	void Close();
	void Start();
	void Stop();
	void ProcessNextBuffer(LPWAVEHDR pwh);
	void SetFormat(LPPCMFORMAT lpPcmFormat);
	void SetFormat(WORD wBitsPerSample,WORD wChannels,DWORD dwSampleRate);
	BOOL IsRecording();
	BOOL IsDeviceOpen();

	void SetBufferFunction(){m_lpData=(void*)&data_holder_temp;}
	
	DWORD GetPosition();
	BOOL Pause();
	BOOL Continue();
	BOOL IsFormatSupported(WAVEFORMATEX wfEx, UINT nDev=WAVE_MAPPER);
	BOOL Spectral_Analysis(void* lpData,LPWAVEHDR pwh);
	HANDLE m_hEvent;
	LPWAVEHDR m_lpWaveHdr;

public:
	float spectrum[FFT_LEN/2];

public:
	BEAT_DETECTOR beat_detector;

protected:
	DataHolder data_holder_temp;
	void *m_lpData;
	HANDLE m_hThread;
	WAVECLASS m_waveClass;
	GLOBALHANDLE m_hWaveInHdr[MAXNUMOFBUFFER];
	LPWAVEHDR m_lpWaveInHdr[MAXNUMOFBUFFER];
	GLOBALHANDLE m_hInBuffer[MAXNUMOFBUFFER];
	LPBYTE m_lpInBuffer [MAXNUMOFBUFFER];
	PCMFORMAT m_PcmFormat;
	WAVEFORMATEX m_WaveFormat;
	HWAVEIN m_hWaveIn;
	BOOL m_bRecording;
	BOOL m_bDeviceOpen;
	DWORD m_dwBufferSize;
};
