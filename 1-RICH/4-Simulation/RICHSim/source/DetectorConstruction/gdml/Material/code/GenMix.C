// generate optical parameters for the database folder
// ..
//
#include "MyReadData.h"
#include "MyReadData.cpp"

typedef vector<vector<double>> DTLIST;

void GenLiF();
void GenC6F14();
void GenQuartz();
void GenArMixGas(double O2ppm,  double H2Oppm);
void WriteXML(DTLIST fAbsData, DTLIST fRefData, DTLIST fRayData, double density, std::map<int, TString> mat, TString matName);

//.. Main ..
void GenMix()
{
    //GenLiF();
    //GenC6F14();
    //GenQuartz();
    GenArMixGas(10., 10.);
}

//.. sub_functions ..
void WriteXML(DTLIST fAbsData, DTLIST fRefData, DTLIST fRayData, double density, std::map<int, TString> mat, TString matName)
{
    TString fileName(matName + ".xml");
    ofstream fout(fileName, ios::out);
    if (!fout.is_open())
    {
        cout << "####> Can't open " << fileName << " to write. Please check!" << endl;
        return;
    }

    //1.0 begin define xml
    fout << "<define>" << endl;

    //1.1 折射率
    fout << "\t<!-- ********** Refractive index begin ********** -->" << endl;
    fout << "\t<matrix name=\"" << matName << "_RINDEX\" coldim=\"2\" values=\"";
    for (int i = 0; i < fRefData.size(); i++)
        fout << "\n\t\t\t" << 1239.8 / fRefData[i][0] << "*eV\t" << fRefData[i][1];
    fout << "\" />" << endl;
    fout << "\t<!-- ********** Refractive index end ********** -->\n\n";

    //1.2 吸收长度
    fout << "\t<!-- ********** Absorption length begin ********** -->" << endl;
    fout << "\t<matrix name=\"" << matName << "_ABS\" coldim=\"2\" values=\"";
    for (int i = 0; i < fAbsData.size(); i++)
        fout << "\n\t\t\t" << 1239.8 / fAbsData[i][0] << "*eV\t" << fAbsData[i][1] << "*cm";
    fout << "\" />" << endl;
    fout << "\t<!-- ********** Absorption length end ********** -->\n\n";

    //1.3 瑞利散射长度
    fout << "\t<!-- ********** Rayleigh scattering length begin ********** -->" << endl;
    fout << "\t<matrix name=\"" << matName << "_RAY\" coldim=\"2\" values=\"";
    for (int i = 0; i < fRayData.size(); i++)
        fout << "\n\t\t\t" << 1239.8 / fRayData[i][0] << "*eV\t" << fRayData[i][1] << "*cm";
    fout << "\" />" << endl;
    fout << "\t<!-- ********** Rayleigh scattering length end ********** -->\n\n";

    //1.4 end define xml
    fout << "</define>\n\n";

    //2. 定义材料
    fout << "<materials>" << endl;

    fout << "\t<material formula=\"" << matName << "\" name=\"" << matName << "\">" << endl;
    fout << "\t\t<property name=\"RINDEX\" ref=\"" << matName << "_RINDEX\" />" << endl;
    fout << "\t\t<property name=\"ABSLENGTH\" ref=\"" << matName << "_ABS\" />" << endl;
    fout << "\t\t<property name=\"RAYLEIGH\" ref=\"" << matName << "_RAY\" />" << endl;
    fout << "\t\t<D value=\"" << density << "\" />" << endl;
    
    for(std::map<int, TString>::const_iterator iter = mat.begin(); iter != mat.end(); ++iter)
        fout <<"\t\t<composite n=\""<<iter->first<<"\" ref=\""<<iter->second<<"\" />"<<endl;
    fout << "\t</material>" << endl;

    fout << "</materials>" << endl;
}

void GenLiF()
{
    MyReadData *fAbsData = new MyReadData();
    fAbsData->ReadData("../database/LiF_abs.txt");

    MyReadData *fRefData = new MyReadData();
    fRefData->ReadData("../database/LiF_ref.txt");

    DTLIST fRayData = {{120, 1000},{220, 1000}};

    double density = 2.64; //kg/cm^3
    std::map<int, TString> mat;
    mat.insert(std::make_pair(1, "Lithium"));
    mat.insert(std::make_pair(1, "Fluorine"));

    WriteXML(fAbsData->GetData(), fRefData->GetData(), fRayData, density, mat, "LiF");
}

void GenC6F14()
{
    MyReadData *fAbsData = new MyReadData();
    fAbsData->ReadData("../database/C6F14_abs.txt");

    MyReadData *fRefData = new MyReadData();
    fRefData->ReadData("../database/C6F14_ref.txt");

    DTLIST fRayData = {{160, 1000},{260, 1000}};

    double density = 1.669; //kg/cm^3
    std::map<int, TString> mat;
    mat.insert(std::make_pair(6, "Carbon"));
    mat.insert(std::make_pair(14, "Fluorine"));

    WriteXML(fAbsData->GetData(), fRefData->GetData(), fRayData, density, mat, "C6F14");
}

void GenQuartz()
{
    MyReadData *fAbsData = new MyReadData();
    fAbsData->ReadData("../database/Quartz_abs.txt");

    MyReadData *fRefData = new MyReadData();
    fRefData->ReadData("../database/Quartz_ref.txt");

    DTLIST fRayData = {{160, 1000},{260, 1000}};

    double density = 2.65; //kg/cm^3
    std::map<int, TString> mat;
    mat.insert(std::make_pair(1, "Silicon"));
    mat.insert(std::make_pair(2, "Oxygen"));

    WriteXML(fAbsData->GetData(), fRefData->GetData(), fRayData, density, mat, "Quartz");
}

void GenArMixGas(double O2ppm, double H2Oppm)
{
    //H2O & O2 吸收截面
    MyReadData *fH2OXsecData = new MyReadData();
    fH2OXsecData->ReadData("../database/H2O_abs.txt");

    MyReadData *fO2XsecData = new MyReadData();
    fO2XsecData->ReadData("../database/O2_abs.txt");

    //Ar相关参数
    MyReadData *fRefData = new MyReadData();
    fRefData->ReadData("../database/Ar_ref.txt");

    DTLIST fAbsData(100, vector<double>(2));
    DTLIST fRayData = {{130, 1000},{330, 1000}};

    //Mix Ar混合气需要根据 H2O 和 O2 计算Abs吸收长度
    for(int i=0; i<100; i++)
    {
        double lambda = 130 + i;
    
        double Xsec1 = fO2XsecData->GetValue(lambda)  * 1E6 * 1E-28; //截面为Mb，换算为m^2
        double Xsec2 = fH2OXsecData->GetValue(lambda)  * 1E6 * 1E-28; 

        double density1 = (O2ppm * 1E-6) / (22.4E-3) * 6.02214129 * 1E23; //ppm定义为单位体积的分子数mol占比
        double density2 = (H2Oppm * 1E-6) / (22.4E-3) * 6.02214129 * 1E23; //1mol理想气体所占的体积都约为22.4升=22.4E-3 m^3

        double absleng1 = 100. / Xsec1 / density1; //cm
        double absleng2 = 100. / Xsec2 / density2; //cm

        absleng1 = (Xsec1 == 0 || density1 == 0) ? 10000 : absleng1;
        absleng2 = (Xsec2 == 0 || density2 == 0) ? 10000 : absleng2;

        double len = absleng1 * absleng2 / (absleng1 + absleng2);

        fAbsData[i][0] = lambda;
        fAbsData[i][1] = len;
    }

    //写入数据
    double density = 1.7841E-3; //kg/cm^3 @ 1atm
    std::map<int, TString> mat;
    mat.insert(std::make_pair(1, "Argon"));

    WriteXML(fAbsData, fRefData->GetData(), fRayData, density, mat, "ArGasMix");
}
