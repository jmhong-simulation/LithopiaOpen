// Recorder.cpp: implementation of the CRecorder class.
//
//////////////////////////////////////////////////////////////////////

//#include "waveInFFT.h"
#include <iostream>
#include <windows.h>
#include "Recorder.h"
#include "Fourier.h"
#include <math.h>

#define mag_sqrd(re,im) (re*re+im*im)
#define Decibels(re,im) ((re==0&&im==0)?(0):10.0*log10(float(mag_sqrd(re,im))))
#define Amplitude(re,im,len) (GetFrequencyIntensity(re,im)/(len))
#define AmplitudeScaled(re,im,len,scale) ((int)Amplitude(re,im,len)%scale)
inline float GetFrequencyIntensity(float re, float im){return sqrt((re*re)+(im*im));}
int mode_0_3=4;// 0 - 8192 , 1 - 4096 , 2 - 2048, 3 - 1024, others - 512

//BOOL Recorder::Spectral_Analysis(void* lpData,LPWAVEHDR pwh)
//{
//	static float finleft[FFT_LEN/2],finright[FFT_LEN/2],fout[FFT_LEN],foutimg[FFT_LEN];
//	
//	DWORD nCount=0;
//	for (DWORD dw=0;dw<FFT_LEN;dw++){
//		finleft[nCount]=(float)((short*)pwh->lpData)[dw++];
//		finright[nCount++]=(float)((short*)pwh->lpData)[dw];
//	
//		finleft[nCount-1]+=finright[nCount-1];
//		finleft[nCount-1]*=0.5f;
//	}
//
//	FFT(FFT_LEN/2,0,finleft,NULL,fout,foutimg);
//
//	float re,im,fmax=-99999.9f,fmin=99999.9f;
//	for(int i=1;i<FFT_LEN/4;i++){
//		re=(float)fout[i];
//		im=(float)foutimg[i];
////		spectrum[i]=Decibels(re,im);
//		spectrum[i]=GetFrequencyIntensity(re,im);
//		if(spectrum[i]>fmax)fmax=(float)spectrum[i];
//		if(spectrum[i]<fmin)fmin=(float)spectrum[i];}
//
//	return TRUE;
//}


BOOL Recorder::Spectral_Analysis(void* lpData,LPWAVEHDR pwh)
{
	float finleft[FFT_LEN/2],finright[FFT_LEN/2],fout[FFT_LEN],foutimg[FFT_LEN];
	static float data_buffer[FFT_LEN];	
	static int in_point=0;
	static bool isReady=false;
	for (DWORD dw=0;dw<pwh->dwBufferLength/4;dw++){
		data_buffer[in_point++]=((float)((short*)pwh->lpData)[dw++]);
		data_buffer[in_point]=((float)((short*)pwh->lpData)[dw]);
		if(!isReady)
			if(in_point+1==(FFT_LEN))isReady=true;
		in_point=(++in_point)%FFT_LEN;}
	
	if(!isReady)return false;

	DWORD nCount=0;
	int out_point=in_point;
	for (DWORD dw=0;dw<FFT_LEN;dw+=2){
		finleft[nCount]=data_buffer[out_point++];
		finright[nCount++]=data_buffer[out_point];
		finleft[nCount-1]+=finright[nCount-1]; 
		finleft[nCount-1]*=0.5f;
		out_point=(++out_point)%FFT_LEN;
	}
	FFT(FFT_LEN/2,0,finleft,NULL,fout,foutimg);
	float re,im,fmax=-99999.9f,fmin=99999.9f;
	for(int i=1;i<FFT_LEN/4;i++){
		re=(float)fout[i];
		im=(float)foutimg[i];
//		spectrum[i]=Decibels(re,im);
		spectrum[i]=GetFrequencyIntensity(re,im);
		if(spectrum[i]>fmax)fmax=(float)spectrum[i];
		if(spectrum[i]<fmin)fmin=(float)spectrum[i];
	}
	//static DWORD start_time=timeGetTime();
	//static int time_count=0;
	//time_count++;
	//if(time_count>100){
	//	printf("%f ",100.0f*1000.0f/(float)(timeGetTime()-start_time));
	//	start_time=timeGetTime();
	//	time_count=0;}
	return TRUE;
	
}
// Recorder.cpp: implementation of the CRecorder class.

//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void CALLBACK waveInProc(HWAVEIN hwi, UINT uMsg, DWORD dwInstance, DWORD dwParam1, DWORD dwParam2)
{
	Recorder *pRecorder=NULL;
	switch(uMsg){
	case WIM_OPEN:
		break;
	case WIM_DATA:/*if((((LPWAVEHDR)dwParam1)->dwFlags)==WHDR_DONE)*/
		pRecorder=(Recorder*)(((LPWAVEHDR)dwParam1)->dwUser);
		pRecorder->m_lpWaveHdr=(LPWAVEHDR)dwParam1;
		SetEvent(pRecorder->m_hEvent);
		break;
	case WIM_CLOSE:
		break;
	default:
		break;}
}

DWORD CALLBACK RecorderThreadFunc(LPVOID lpThreadData)
{
	Recorder *pRecorder=NULL;
	pRecorder=(Recorder*)lpThreadData;
	while(pRecorder->IsRecording()){
		
	//static DWORD start_time=timeGetTime();
	//static int time_count=0;
	//time_count++;
	//if(time_count>100){
	//	printf("%f ",100.0f*1000.0f/(float)(timeGetTime()-start_time));
	//	start_time=timeGetTime();
	//	time_count=0;}

		WaitForSingleObject(pRecorder->m_hEvent,INFINITE);
		pRecorder->ProcessNextBuffer(pRecorder->m_lpWaveHdr);}
	return 0;
}
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Recorder::Recorder(int nBufferLength)
{
	m_bRecording=FALSE;
	m_bDeviceOpen=FALSE;
	m_PcmFormat.wBitsPerSample=16;
	m_PcmFormat.wChannels=2;
	m_PcmFormat.dwSampleRate=44100;	
	
	if(mode_0_3==0)
		m_dwBufferSize=(nBufferLength*m_PcmFormat.wChannels*m_PcmFormat.wBitsPerSample/8);//8192
	else if(mode_0_3==1)
		m_dwBufferSize=(nBufferLength*m_PcmFormat.wChannels*m_PcmFormat.wBitsPerSample/16);//4096
	else if(mode_0_3==2)
		m_dwBufferSize=(nBufferLength*m_PcmFormat.wChannels*m_PcmFormat.wBitsPerSample/32);//2048
	else if(mode_0_3==3)
		m_dwBufferSize=(nBufferLength*m_PcmFormat.wChannels*m_PcmFormat.wBitsPerSample/64);//1024
	else
		m_dwBufferSize=(nBufferLength*m_PcmFormat.wChannels*m_PcmFormat.wBitsPerSample/128);//512


//	m_PcmFormat.dwSampleRate=44100;
//	m_dwBufferSize=(nBufferLength*m_PcmFormat.wChannels*m_PcmFormat.wBitsPerSample/16);

	m_lpWaveHdr=NULL;
	m_hEvent=NULL;
	m_hThread=NULL;
	for(int i=0;i<MAXNUMOFBUFFER;i++){
		m_hWaveInHdr[i]=NULL;
		m_hInBuffer[i]=NULL;}
}

Recorder::Recorder(PCMFORMAT pcm,int nBufferLength)
{
	m_bRecording=FALSE;
	m_bDeviceOpen=FALSE;
	m_PcmFormat.wBitsPerSample=pcm.wBitsPerSample;
	m_PcmFormat.wChannels=pcm.wChannels;
	m_PcmFormat.dwSampleRate=pcm.dwSampleRate;
	m_dwBufferSize=(nBufferLength*m_PcmFormat.wChannels*m_PcmFormat.wBitsPerSample/8);
	m_lpWaveHdr=NULL;
	
	m_hEvent=NULL;
	m_hThread=NULL;
	for(int i=0;i<MAXNUMOFBUFFER;i++){
		m_hWaveInHdr[i]=NULL;
		m_hInBuffer[i]=NULL;}
}

Recorder::Recorder(WORD wBitsPerSample,WORD wChannels,DWORD dwSampleRate,int nBufferLength)
{
	m_bRecording=FALSE;
	m_bDeviceOpen=FALSE;
	m_PcmFormat.wBitsPerSample=wBitsPerSample;
	m_PcmFormat.wChannels=wChannels;
	m_PcmFormat.dwSampleRate=dwSampleRate;
	m_dwBufferSize=(nBufferLength*m_PcmFormat.wChannels*m_PcmFormat.wBitsPerSample/8);
	m_lpWaveHdr=NULL;
	m_hEvent=NULL;
	m_hThread=NULL;
	for(int i=0;i<MAXNUMOFBUFFER;i++){
		m_hWaveInHdr[i]=NULL;
		m_hInBuffer[i]=NULL;}
}

Recorder::~Recorder()
{
	if(m_bRecording)Stop();
	if(m_bDeviceOpen)Close();
}

void Recorder::Open(DWORD dwCallBack, DWORD dwCallbackType,MCIDEVICEID wMCIDeviceID)
{
    if(m_bDeviceOpen){
//		TRACE("Device Already Opened. Please Stop Recorder before attempting to Open\n");
		std::cout << "Device Already Opened. Please Stop Recorder before attempting to Open" << std::endl;
//		MessageBox(NULL,"Device Already Opend.","Error", MB_OK);
		return;}
	if (dwCallBack==NULL)dwCallBack=(DWORD)waveInProc;

	for(int i=0; i<MAXNUMOFBUFFER;i++){
		m_hWaveInHdr[i]=GlobalAlloc(GHND|GMEM_SHARE,sizeof(WAVEHDR));
		m_lpWaveInHdr[i]=(LPWAVEHDR)GlobalLock(m_hWaveInHdr[i]);
		m_hInBuffer[i]=GlobalAlloc(GHND|GMEM_SHARE,m_dwBufferSize);
		m_lpInBuffer[i]=(LPBYTE)GlobalLock(m_hInBuffer[i]);
		m_lpWaveInHdr[i]->lpData=(LPSTR)m_lpInBuffer[i];
		m_lpWaveInHdr[i]->dwBufferLength=m_dwBufferSize;
		m_lpWaveInHdr[i]->dwBytesRecorded=0L;
		m_lpWaveInHdr[i]->dwUser=(DWORD)(void*)this;
		m_lpWaveInHdr[i]->dwFlags=0L;
		m_lpWaveInHdr[i]->dwLoops=1L;
		m_lpWaveInHdr[i]->lpNext=NULL;
		m_lpWaveInHdr[i]->reserved=0L;}

	m_WaveFormat.wFormatTag = WAVE_FORMAT_PCM;
	m_WaveFormat.nChannels = m_PcmFormat.wChannels;
	m_WaveFormat.wBitsPerSample = m_PcmFormat.wBitsPerSample;
	m_WaveFormat.nSamplesPerSec = m_PcmFormat.dwSampleRate;
	m_WaveFormat.nBlockAlign = m_WaveFormat.nChannels * m_WaveFormat.wBitsPerSample/8;
	m_WaveFormat.nAvgBytesPerSec = m_WaveFormat.nBlockAlign * m_WaveFormat.nSamplesPerSec;
	m_waveClass.lpData = this;
//	if(waveInOpen((LPHWAVEIN)&m_waveClass,wMCIDeviceID,&m_WaveFormat,(DWORD)dwCallBack,(DWORD)0L, dwCallbackType)||m_waveClass.hWave==0)return;
	if(waveInOpen((LPHWAVEIN)&m_waveClass,wMCIDeviceID,&m_WaveFormat,(DWORD_PTR)waveInProc,(DWORD)0L,dwCallbackType)||m_waveClass.hWave==0)return;
	m_waveClass.lpData = this;
	m_hWaveIn = (HWAVEIN)m_waveClass.hWave;
	m_hEvent = CreateEvent(NULL,FALSE,FALSE,NULL);
    m_bDeviceOpen=TRUE;
}

void Recorder::Start()
{
	if(!m_bDeviceOpen){
//		TRACE("Device not Opened. Please open device before attempting to Start\n");
//		MessageBox(NULL,"Device Not Opened","Error",MB_OK);

		std::cout << "device not opened" << std::endl;
		return;}

	for(int i=0; i<MAXNUMOFBUFFER; i++){
		// Prepare wave in header
		if(waveInPrepareHeader(m_hWaveIn, m_lpWaveInHdr[i], sizeof(WAVEHDR)) != MMSYSERR_NOERROR)return;

		// Add buffer into recording queue
		if(waveInAddBuffer(m_hWaveIn, m_lpWaveInHdr[i], sizeof(WAVEHDR)) != MMSYSERR_NOERROR)return;}

	// Begin sampling
	m_bRecording = TRUE;
	m_hThread = CreateThread(NULL,NULL,RecorderThreadFunc,this,NULL,NULL);
	waveInStart(m_hWaveIn);
//	ASSERT(m_hThread!=NULL);
	SetPriorityClass(m_hThread,REALTIME_PRIORITY_CLASS);
	SetThreadPriority(m_hThread,THREAD_PRIORITY_HIGHEST);
}

void Recorder::Stop()
{
	if(!m_bDeviceOpen || !m_bRecording)return;
	if(waveInStop(m_hWaveIn) != MMSYSERR_NOERROR)return;
	else m_bRecording=FALSE;
}

void Recorder::Close()
{
	if(m_bRecording)Stop();
	if (m_hThread != NULL)CloseHandle(m_hThread);
	if(m_bDeviceOpen)waveInClose(m_hWaveIn);
	for(int i=0; i<MAXNUMOFBUFFER; i++){
		if (m_hWaveInHdr[i] != NULL){
			if (GlobalUnlock(m_hWaveInHdr[i]))GlobalFree(m_hWaveInHdr[i]);
			if (GlobalUnlock(m_hInBuffer[i]))GlobalFree(m_hInBuffer[i]);
			m_hWaveInHdr[i] = NULL;
			m_hInBuffer[i] = NULL;}}
	m_bDeviceOpen=FALSE;
	m_bRecording = FALSE;
	m_hThread = NULL;
}
void Recorder::SetFormat(LPPCMFORMAT lpPcmFormat)
{
//	ASSERT(m_bDeviceOpen==false);
	m_PcmFormat.wBitsPerSample=lpPcmFormat->wBitsPerSample;
	m_PcmFormat.wChannels=lpPcmFormat->wChannels;
	m_PcmFormat.dwSampleRate=lpPcmFormat->dwSampleRate;
}
void Recorder::SetFormat(WORD wBitsPerSample,WORD wChannels,DWORD dwSampleRate)
{
//	ASSERT(m_bDeviceOpen==false);
	m_PcmFormat.wBitsPerSample=wBitsPerSample;
	m_PcmFormat.wChannels=wChannels;
	m_PcmFormat.dwSampleRate=dwSampleRate;
}
BOOL Recorder::IsRecording()
{
	return m_bRecording;
}
BOOL Recorder::IsDeviceOpen()
{
	return m_bDeviceOpen;
}

void Recorder::ProcessNextBuffer(LPWAVEHDR pwh)
{
	if(Spectral_Analysis(m_lpData,pwh))beat_detector.Update(spectrum);
	waveInUnprepareHeader(m_hWaveIn,pwh,sizeof(WAVEHDR));
	waveInPrepareHeader (m_hWaveIn,pwh,sizeof(WAVEHDR));
	waveInAddBuffer(m_hWaveIn,pwh,sizeof(WAVEHDR));
}

DWORD Recorder::GetPosition()
{
	if(m_hWaveIn){
		MMTIME mmtime;
		mmtime.wType = TIME_SAMPLES;
		if (waveInGetPosition(m_hWaveIn, &mmtime, sizeof(MMTIME)) != MMSYSERR_NOERROR)return -1;
		else return mmtime.u.sample;}
	return -1;
}

BOOL Recorder::Pause()
{
	if(m_hWaveIn){
		if(waveInStop(m_hWaveIn)==MMSYSERR_NOERROR){
			m_bRecording=FALSE;
			return TRUE;}}
	return FALSE;
}

BOOL Recorder::Continue()
{
	if(m_hWaveIn){
		if(waveInStart(m_hWaveIn) == MMSYSERR_NOERROR){
			m_bRecording = FALSE;
			return TRUE;}}
	return FALSE;
}

BOOL Recorder::IsFormatSupported(WAVEFORMATEX wfEx, UINT nDev)
{
	MMRESULT mm = waveInOpen(0,nDev,&wfEx,0,0,WAVE_FORMAT_QUERY);
	return (BOOL)mm == MMSYSERR_NOERROR;		
}
