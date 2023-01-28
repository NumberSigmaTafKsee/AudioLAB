// wxWidgets "Hello world" Program
// For compilers that support precompilation, includes "wx/wx.h".
#pragma once

#include "./window/Window.h"
class SoundEditor : public wxApp
{
public:
    Window* window;
    virtual bool OnInit();
};
