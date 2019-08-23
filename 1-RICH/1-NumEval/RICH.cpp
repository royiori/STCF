#include "MyGuiActionClass.h"
#include "MyMainFrameGui.h"

int main(int argc, char **argv)
{
    if (argc == 1) //默认采用gui模式
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

    //采用batch模式
    //输入参数： 保存的文件路径及文件名
    gMyGuiActionClass = new MyGuiActionClass(BATCH);
    gMyGuiActionClass->DoReadBatchFile(argv[1]); //batch file

    return 0;
}
