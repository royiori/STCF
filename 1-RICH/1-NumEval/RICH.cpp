#include "MyGuiActionClass.h"
#include "MyMainFrameGui.h"

int main(int argc, char **argv)
{

    TApplication *theApp;
    theApp = new TApplication("App", &argc, argv);

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    gStyle->SetStatFont(42);

    gMyGuiActionClass = new MyGuiActionClass();
    gMyMainFrameGui = new MyMainFrameGui(gClient->GetRoot(), 1250, 600);
    
    // run ROOT application
    theApp->Run();

    delete theApp;
    return 0;
}