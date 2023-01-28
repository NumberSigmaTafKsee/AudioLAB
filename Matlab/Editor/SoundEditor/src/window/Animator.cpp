#include "Animator.h"

Animator::Animator(wxWindow* window){
    this->window = window;

    this->rendering = false;
}

Animator::~Animator(){

}

void Animator::TimerStart(){
    Start(50);
}

double maxTime = 0;

void Animator::Notify(){
    Draw();
}

void Animator::ActivateRenderLoop(bool on){
    if(on && !rendering)
    {   
        panel = new wxPanel(window, wxID_ANY, {20, 350}, {300, 350});

        TimerStart();
        rendering = true;
    }
    else if(!on && rendering)
    {   
        tslist.clearAll();
        Stop();
        rendering = false;
    }
}

double last_r=0;

void Animator::Draw(){
    if(rendering){
        static wxClientDC dc(panel);

        static int width = 2;
        
        wxColor bckg = window->GetBackgroundColour();
        bckg.Set(bckg.Red(), bckg.Green(), bckg.Blue(), wxALPHA_OPAQUE);
        dc.SetBrush(wxBrush(bckg));
        dc.SetPen(wxPen(bckg, 3));
        dc.DrawRectangle(50, 0, 200, 50);

        tslist.forEach([](sample s, int idx){
            int height = std::min(50, (int)(600*val(s)/__SHRT_MAX__));

            int green = 150+freq(s)*25;
            int blue = std::min(255, std::max(green-255, 0));
            green = std::min(green, 255);
            
            dc.SetBrush(wxBrush(wxColor(0, green, blue)));
            dc.SetPen(wxPen(wxColor(0, green, blue)));
            dc.DrawRectangle(50+idx*width, 50-height, width, height);
        });

        sample s = tslist.back();
        int r = std::min(150, (int) std::abs(val(s)/40));

        int green = 150+freq(s)*25;
        int blue = std::min(255, std::max(green-255, 0));
        green = std::min(green, 255);

        dc.SetBrush(wxBrush(window->GetBackgroundColour()));
        dc.SetPen(wxPen(window->GetBackgroundColour(), 3));
        dc.DrawCircle(150, 200, last_r);

        dc.SetBrush(wxBrush(wxColor(0, green, blue)));
        dc.SetPen(wxPen(wxColor(0, green, blue)));
        dc.DrawCircle(150, 200, r);

        last_r = r;
    }
}

void Animator::TSList::push(sample& s){
    std::lock_guard<std::mutex> dataGuard(dataMutex);

    if(samples.size() >= 100)
        samples.pop_front();

    samples.push_back(s);
}

void Animator::TSList::remove_front(){
    std::lock_guard<std::mutex> dataGuard(dataMutex);

    samples.pop_front();
}

void Animator::TSList::clearAll(){
    std::lock_guard<std::mutex> dataGuard(dataMutex);

    samples.clear();
}

sample& Animator::TSList::back(){
    std::lock_guard<std::mutex> dataGuard(dataMutex);

    return samples.back();
}

void Animator::TSList::forEach(void (*func)(sample, int)){
    std::lock_guard<std::mutex> dataGuard(dataMutex);

    std::list<sample>::iterator it;
    int idx = 0;

    for (it = samples.begin(); it != samples.end(); ++it)
        func(*it, idx++);
}

size_t Animator::TSList::size(){
    std::lock_guard<std::mutex> dataGuard(dataMutex);

    return samples.size();
}
