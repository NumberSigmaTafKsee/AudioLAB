#pragma once

#include <iostream>
#include <iterator>
#include <list>
#include <mutex>
#include <cmath>

#include "wx/wx.h"

#define val(s) (s).first // Macro for sample value
#define freq(s) (s).second // Macro for sample freq

typedef std::pair<double, double> sample; // Sample definition

class Animator : public wxTimer{
public:
    Animator(wxWindow* window);
    ~Animator();

    void ActivateRenderLoop(bool on);

    void Notify();
    void TimerStart();

    void Draw();

    struct TSList{
        std::mutex dataMutex;

        std::list<sample> samples;

        void push(sample& s);

        void remove_front();

        void clearAll();

        sample& back();

        void forEach(void (*func)(sample, int));

        size_t size();
    } tslist;

private:
    wxPanel* panel;
    wxWindow* window;
    bool rendering = false;
};