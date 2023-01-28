#include "Window.h"

BEGIN_EVENT_TABLE (Window, wxFrame)
    EVT_BUTTON ( PLAY_BUTTON, Window::Play)
    EVT_BUTTON ( RESTART_BUTTON, Window::Restart)
    EVT_BUTTON ( CONTINUE_BUTTON, Window::Continue)
    EVT_BUTTON ( STOP_BUTTON, Window::Stop)
    EVT_BUTTON ( LOAD_FILE_BUTTON, Window::LoadFile)
END_EVENT_TABLE()


std::future<void> animation;

Window::Window()
    : wxFrame(PARENT, ID, TITLE, WINODW_POSITION, WINDOW_SIZE, WINDOW_STYLE)
{   
    animator = new Animator(this);

    reader = new WAVReader();
    player = new Player(animator);
    
    sound = nullptr;

    DrawSoundData(sound);
}

Window::~Window(){
    delete reader;
    delete player;

    animator->ActivateRenderLoop(false);
    delete animator;
}

std::string soundPath = "";

void Window::LoadFile(wxCommandEvent& event){
     wxFileDialog 
        openFileDialog(this, _("Load a WAV file"), "", "",
                       "Wav files (*.wav)|*.wav", wxFD_OPEN);
    if (openFileDialog.ShowModal() == wxID_CANCEL)
        return; 

    soundPath = openFileDialog.GetPath().ToStdString();

    sound = reader->ReadFileData(soundPath);

    DrawSoundData(sound);
}

void Window::DrawSoundData(WAV_FILE* sound){
    if(!sound){
        wxButton* LoadFileButton = new wxButton(this, LOAD_FILE_BUTTON, _T("Find file"),
        wxPoint(200, 200), wxSize(200, 50), 0);
        return;
    }

    this->DestroyChildren();

    wxStaticText* staticText = new wxStaticText(this, wxID_ANY, "RIFF_chunk_ID", wxPoint(30, 10));
    wxTextCtrl* RIFF_chunk_ID = new wxTextCtrl(this, wxID_ANY, "", wxPoint(30, 30), wxSize(300, 50));
	RIFF_chunk_ID->SetValue(std::to_string(sound->header->RIFF_chunk_ID));
	
    wxStaticText* staticText2 = new wxStaticText(this, wxID_ANY, "RIFF_chunk_size", wxPoint(30, 90));
    wxTextCtrl* RIFF_chunk_size = new wxTextCtrl(this, wxID_ANY, "", wxPoint(30, 110), wxSize(300, 50));
	RIFF_chunk_size->SetValue(std::to_string(sound->header->RIFF_chunk_size));

    wxStaticText* staticText3 = new wxStaticText(this, wxID_ANY, "RIFF_format", wxPoint(30, 170));
    wxTextCtrl* RIFF_format = new wxTextCtrl(this, wxID_ANY, "", wxPoint(30, 190), wxSize(300, 50));
	RIFF_format->SetValue(std::to_string(sound->header->RIFF_format));

    wxStaticText* staticText4 = new wxStaticText(this, wxID_ANY, "FMT_chunk_ID", wxPoint(350, 10));
    wxTextCtrl* FMT_chunk_ID = new wxTextCtrl(this, wxID_ANY, "", wxPoint(350, 30), wxSize(300, 50));
	FMT_chunk_ID->SetValue(std::to_string(sound->header->FMT_chunk_ID));

    wxStaticText* staticText5 = new wxStaticText(this, wxID_ANY, "FMT_chunk_ID", wxPoint(350, 90));
    wxTextCtrl* FMT_chunk_size = new wxTextCtrl(this, wxID_ANY, "", wxPoint(350, 110), wxSize(300, 50));
	FMT_chunk_size->SetValue(std::to_string(sound->header->FMT_chunk_size));

    wxStaticText* staticText6 = new wxStaticText(this, wxID_ANY, "FMT_chunk_ID", wxPoint(350, 170));
    wxTextCtrl* FMT_audio_format = new wxTextCtrl(this, wxID_ANY, "", wxPoint(350, 190), wxSize(300, 50));
	FMT_audio_format->SetValue(std::to_string(sound->header->FMT_audio_format));

    wxStaticText* staticText7 = new wxStaticText(this, wxID_ANY, "FMT_num_channels", wxPoint(350, 250));
    wxTextCtrl* FMT_num_channels = new wxTextCtrl(this, wxID_ANY, "", wxPoint(350, 270), wxSize(300, 50));
	FMT_num_channels->SetValue(std::to_string(sound->header->FMT_num_channels));

    wxStaticText* staticText8 = new wxStaticText(this, wxID_ANY, "FMT_sample_rate", wxPoint(350, 330));
    wxTextCtrl* FMT_sample_rate = new wxTextCtrl(this, wxID_ANY, "", wxPoint(350, 350), wxSize(300, 50));
	FMT_sample_rate->SetValue(std::to_string(sound->header->FMT_sample_rate));

    wxStaticText* staticText9 = new wxStaticText(this, wxID_ANY, "FMT_byte_rate", wxPoint(350, 410));
    wxTextCtrl* FMT_byte_rate = new wxTextCtrl(this, wxID_ANY, "", wxPoint(350, 430), wxSize(300, 50));
	FMT_byte_rate->SetValue(std::to_string(sound->header->FMT_byte_rate));

    wxStaticText* staticText10 = new wxStaticText(this, wxID_ANY, "FMT_block_align", wxPoint(350, 490));
    wxTextCtrl* FMT_block_align = new wxTextCtrl(this, wxID_ANY, "", wxPoint(350, 510), wxSize(300, 50));
	FMT_block_align->SetValue(std::to_string(sound->header->FMT_block_align));

    wxStaticText* staticText11 = new wxStaticText(this, wxID_ANY, "FMT_bits_per_sample", wxPoint(350, 570));
    wxTextCtrl* FMT_bits_per_sample = new wxTextCtrl(this, wxID_ANY, "", wxPoint(350, 590), wxSize(300, 50));
	FMT_bits_per_sample->SetValue(std::to_string(sound->header->FMT_bits_per_sample));

    wxStaticText* staticText12 = new wxStaticText(this, wxID_ANY, "DATA_chunk_size", wxPoint(670, 10));
    wxTextCtrl* DATA_chunk_size = new wxTextCtrl(this, wxID_ANY, "", wxPoint(670, 30), wxSize(300, 50));
	DATA_chunk_size->SetValue(std::to_string(sound->header->DATA_chunk_size));

    wxButton* PlayButton = new wxButton(this, PLAY_BUTTON, _T("Play sound"),
    wxPoint(720, 110), wxSize(200, 50), 0);

    wxButton* StopButton = new wxButton(this, STOP_BUTTON, _T("Stop sound"),
    wxPoint(720, 170), wxSize(200, 50), 0);

    wxButton* ContinueButton = new wxButton(this, CONTINUE_BUTTON, _T("Continue sound"),
    wxPoint(720, 230), wxSize(200, 50), 0);

    wxButton* RestartButton = new wxButton(this, RESTART_BUTTON, _T("Restart sound"),
    wxPoint(720, 290), wxSize(200, 50), 0);

    wxButton* LoadFileButton = new wxButton(this, LOAD_FILE_BUTTON, _T("Find file"),
    wxPoint(720, 350), wxSize(200, 50), 0);
}   

void Window::Play(wxCommandEvent& event ){
    animator->ActivateRenderLoop(true);

    player->playSound(soundPath.c_str());
}

void Window::Restart(wxCommandEvent& event ){
    player->restartSound();
}

void Window::Continue(wxCommandEvent& event ){
    player->continueSound();
}

void Window::Stop(wxCommandEvent& event ){
    player->stopSound();
}