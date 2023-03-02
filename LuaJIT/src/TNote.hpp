class Note {
    
public:
    
    Note();
    ~Note();
    
    void setVelocity(int velocity);
    int getVelocity() const;
    void setMidiNote(int note);
    int getMidiNote() const;
    void setPlaying(bool playing);
    bool isPlaying() const;

    
private:
    int velocity;
    int midiNote;
    bool playing;

    
};

Note::Note() {
    this->velocity = 0;
    this->midiNote = 0;
    this->playing = false;
}

Note::~Note() {
    
}

void Note::setVelocity(int velocity) {
    this->velocity = velocity;
}

int Note::getVelocity() const {
    return this->velocity;
}

void Note::setMidiNote(int note) {
    this->midiNote = note;
}

int Note::getMidiNote() const {
    return this->midiNote;
}

void Note::setPlaying(bool playing) {
    this->playing = playing;
}

bool Note::isPlaying() const {
    return this->playing;
}