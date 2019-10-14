#include "../inc/MyBeamTestDetector.h"

void read()
{
    TFile f("20190816184039063v2.root");
    TTree *T = (TTree *)f.Get("tree2");
    MyBeamTestData *event = 0;
    T->SetBranchAddress("event", &event);
    Long64_t nentries = T->GetEntries();

    for (Long64_t ev = 0; ev < 10; ev++)
    {
        T->GetEntry(ev);
        cout << " \n\nid : " << ev
             << " \nEvent: " << event->event
             << " \nBoard: " << event->board.size() << " (" << event->board[0] << ", " << event->board[1] << ")"
             << " \nChip: " << event->chip.size() << " (" << event->chip[0] << ", " << event->chip[1] << ")"
             << " \nChannel: " << event->channel.size() << " (" << event->channel[0] << ", " << event->channel[1] << ")"
             << " \nwave: " << event->wave.size() << " " << event->wave[0].size() << endl
             << endl;
    }
}