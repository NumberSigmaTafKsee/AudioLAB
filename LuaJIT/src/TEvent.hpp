class Event {
  
public:
    
    enum Type {
        GATE,
        NOTE,
        CONTROLLER,
        PITCH,
        CLOCK,
        CLOCK_START,
        CLOCK_STOP
    };
    
    Event();
    Event(Event* e);
    Event(juce::String name, Type type);
    ~Event();
    
  juce::String getName();
    Type getType();
    int getValue();
    void setValue(int value);
    int getNote();
    void setNote(int note);
    int getNumber();
    void setNumber(int num);
private:
    
    Type type;
    juce::String name;
    int number = 1;
    int value = 0;
    int note = 0;
};



Skip to content
Pull requests
Issues
Marketplace
Explore
@TomCoalfish
mpue /
Synthlab
Public

Code
Issues 11
Pull requests
Actions
Projects
Wiki
Security

    Insights

Synthlab/Source/AudioEngine/Event.cpp
@InSL
InSL use juce namespace
Latest commit 28ec617 on May 4, 2018
History
2 contributors
@mpue
@InSL
62 lines (46 sloc) 830 Bytes
//
//  Event.cpp
//  Synthlab - App
//
//  Created by Matthias Pueski on 05.04.18.
//

#include "Event.h"

using juce::String;

Event::Event() {
    
}

Event::Event(Event* e){
    this->name = e->getName();
    this->type = e->getType();
    this->value = e->getValue();
    this->note = e->getNote();
}

Event::Event(String name, Type type) {
    this->name = name;
    this->type = type;
}

Event::~Event() {
    
}

String Event::getName() {
    return name;
}

Event::Type Event::getType() {
    return type;
}

void Event::setValue(int value) {
    this->value = value;
}

int Event::getValue() {
    return value;
}

void Event::setNote(int note) {
    this->note = note;
}

int Event::getNote(){
    return note;
}

void Event::setNumber(int num) {
    this->number = num;
}

int Event::getNumber() {
    return number;
}
