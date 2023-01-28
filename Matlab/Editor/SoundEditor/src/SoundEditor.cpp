#include "SoundEditor.h"

#include "./window/Window.h"

wxIMPLEMENT_APP(SoundEditor);

bool SoundEditor::OnInit()
{
    this->window = new Window();
    this->window->Show(true);
    
    return true;
}

