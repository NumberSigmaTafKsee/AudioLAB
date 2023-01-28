#pragma once

#include <wx/wxprec.h>
#ifndef WX_PRECOMP
#include <wx/wx.h>
#endif

#include <wx/stattext.h>
#include <wx/pen.h>
#include <wx/filedlg.h>

#include <iostream>
#include <future>
#include <thread>
#include <chrono>

#include "../audio/Player.h"
#include "Animator.h"

#define PLAY_BUTTON 0xb1
#define RESTART_BUTTON 0xb2
#define CONTINUE_BUTTON 0xb3
#define STOP_BUTTON 0xb4
#define LOAD_FILE_BUTTON 0xb5

#define PARENT NULL
#define ID wxID_ANY
#define TITLE "Sound Editor"
#define WINODW_POSITION wxDefaultPosition
#define WINDOW_SIZE wxSize(1000, 700) 
#define WINDOW_STYLE wxDEFAULT_FRAME_STYLE | wxMAXIMIZE

class Window : public wxFrame
{
public:
    Window();
    ~Window();

    void Play(wxCommandEvent& event );
    void Restart(wxCommandEvent& event );
    void Continue(wxCommandEvent& event );
    void Stop(wxCommandEvent& event );
    void LoadFile(wxCommandEvent& event );

    void DrawAnimation(long long value, double freq);
private:
    WAVReader* reader = nullptr;
	WAV_FILE* sound = nullptr;

    Player* player;
    Animator* animator;

    void DrawSoundData(WAV_FILE* sound);

    DECLARE_EVENT_TABLE();
};