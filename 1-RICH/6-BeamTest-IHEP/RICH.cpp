#include "MyGuiMainFrame.h"
#include "MyGuiMainAction.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#define MAXBUFSIZE 1024

int main(int argc, char **argv)
{
    char resolved_path[100];
    realpath(argv[0], resolved_path);

    TApplication *theApp;
    theApp = new TApplication("App", &argc, argv);

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    gStyle->SetStatFont(42);

    gMyGuiMainAction = new MyGuiMainAction(gSystem->DirName(resolved_path));
    gMyGuiMainFrame = new MyGuiMainFrame(gClient->GetRoot(), 1250, 600);

    // run ROOT application
    theApp->Run();

    delete theApp;

    return 0;
}
