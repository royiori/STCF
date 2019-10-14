
#include "MyReadData.h"

using namespace std;

//______________________________________________
// 1. read data
void MyReadData::ReadData(TString fn)
{
    // load file
    ifstream fp(fn);
    if (!fp.is_open())
    {
        status = -1;
        cout << "###Error: Can't read " << fn << "!" << endl;
        return;
    }

    data.clear();
    label.clear();

    TString line;
    while (!fp.eof())
    {
        line.ReadToDelim(fp);
        if (line.BeginsWith("#"))
            saveLabels(line);
        else
            saveRecord(line);
    }

    // add protection
    nrec = (int)data.size();

    ncol = 999;
    for (int i = 0; i < nrec; i++)
        ncol = (ncol > (int)data[i].size()) ? (int)data[i].size() : ncol;

    if ((int)label.size() < ncol)
        for (int i = (int)label.size(); i < ncol; i++)
            label.push_back("NONAME");

    status = 1;
    fileName = fn;
}

//...
void MyReadData::trimLine(TString &line)
{
    int trimflag = 1;
    while (trimflag > 0)
    {
        line.Remove(TString::kTrailing, '\n');
        line.Remove(TString::kBoth, ' ');
        line.Remove(TString::kBoth, '\t');
        trimflag = 0;
        if (line.BeginsWith(TString('\n')))
            trimflag++;
        if (line.EndsWith(TString('\n')))
            trimflag++;
        if (line.BeginsWith(" "))
            trimflag++;
        if (line.EndsWith(" "))
            trimflag++;
        if (line.BeginsWith("\t"))
            trimflag++;
        if (line.EndsWith("\t"))
            trimflag++;
    }
}

//...
void MyReadData::saveLabels(TString line)
{
    if (!line.BeginsWith("#LABEL"))
        return;

    line.ReplaceAll("#LABEL:", "");
    line.ReplaceAll("#LABEL", "");
    trimLine(line);
    line += sep;

    label.clear();
    while (line.Length() > 0 && line.Index(sep)>0)
    {
        TString head = line;
        line.Remove(0, line.Index(sep) + 1);
        head.Remove(head.Index(sep), head.Length());
        trimLine(head);
        if (head.Length() > 0)
            label.push_back(head);
        if(line.Index(sep)==-1)
            label.push_back(line);
    }
}

//...
void MyReadData::saveRecord(TString line)
{
    trimLine(line);
    line += sep;

    vector<double> _dat;
    while (line.Length() > 0)
    {
        TString head = line;
        head.Remove(head.Index(sep), head.Length());
        line.Remove(0, line.Index(sep) + 1);
        trimLine(head);
        if (head.Length() > 0)
            _dat.push_back(head.Atof());
    }

    if (_dat.size() == 0)
        return;

    data.resize(data.size() + 1);
    data[data.size() - 1].resize(_dat.size());
    for (int i = 0; i < (int)_dat.size(); i++)
        data[data.size() - 1][i] = _dat[i];
}

//______________________________________________
// 2. get info
void MyReadData::PrintInfo(int printlevel)
{
    cout << "MyReadData::PrintInfo: Data has read from " << fileName << endl;
    cout << "                       Data Size : [" << nrec << "][" << ncol << "]" << endl;
    cout << "--> Label: ";
    for (int i = 0; i < ncol; i++)
        cout << label[i] << "\t";
    cout << endl;

    if (printlevel <= 1)
        return;

    cout << "--> Data: \n";
    for (int i = 0; i < nrec; i++)
    {
        cout << " " << i << " -->\t";
        for (int j = 0; j < ncol; j++)
            cout << data[i][j] << ", ";
        cout << endl;
    }
}

//______________________________________________
// 3. sort data
void MyReadData::SetIndex(int _idx)
{
    idx = _idx;
    idx = (idx < 0) ? 0 : idx;
    idx = (idx >= ncol) ? 0 : idx;

    if (verbose)
        cout << "MyReadData::SetIndex: Index is set to : " << idx << "." << endl;
}

void MyReadData::SortData()
{
    SetIndex(idx);

    int m = data.size();
    for (int i = 1; i < (int)data.size(); i++)
    {
        m -= 1;
        for (int j = 0; j < m; j++)
        {
            if (data[j][idx] > data[j + 1][idx])
            {
                for (int k = 0; k < (int)data[j].size(); k++)
                {
                    double tmp;
                    tmp = data[j][k];
                    data[j][k] = data[j + 1][k];
                    data[j + 1][k] = tmp;
                }
            }
        }
    }
    if (verbose)
        cout << "MyReadData::SortData: Data is sorted according to column " << idx << "." << endl;
}

//______________________________________________
// 4. get data
double MyReadData::GetData(int xrec, int yrec) // data[x-record][y-record]
{
    if (0 <= xrec && xrec < (int)data.size())
        if (0 <= yrec && yrec < (int)data[xrec].size())
            return data[xrec][yrec];
    return 0;
}

double MyReadData::GetValue(double xval, int yrec) // in index-column find the x-record according to xval, and get [xrec][yrec]
{
    if (xval < data[0][idx])
        return data[0][yrec];
    if (xval > data[nrec - 1][idx])
        return data[nrec - 1][yrec];

    double val = 0;
    for (int i = 0; i < nrec - 1; i++)
        if (data[i][idx] <= xval && xval < data[i + 1][idx])
        {
            val = data[i][yrec] + (data[i + 1][yrec] - data[i][yrec]) * (xval - data[i][idx]) / (data[i + 1][idx] - data[i][idx]);
            break;
        }

    return val;
}

//______________________________________________
// 5. get histogram
TGraph *MyReadData::GetGraph(int xcol, int ycol)
{
    if (g != NULL)
        delete g;

    if (xcol < 0 || xcol >= ncol)
        xcol = 0;
    if (ycol < 0 || ycol >= ncol)
        ycol = ncol - 1;

    g = new TGraph(nrec);
    for (int i = 0; i < nrec; i++)
        g->SetPoint(i, data[i][xcol], data[i][ycol]);

    g->SetMarkerStyle(8);
    g->SetMarkerSize(1);
    g->GetXaxis()->SetTitle(label[xcol]);
    g->GetYaxis()->SetTitle(label[ycol]);
    return g;
}

TH1 *MyReadData::GetHistogram(int col)
{
    if (h != NULL)
        delete h;

    if (col < 0 || col >= ncol)
        col = 0;

    double xmin = data[0][col];
    double xmax = data[0][col];
    for (int i = 1; i < nrec; i++)
    {
        xmin = (xmin > data[i][col]) ? data[i][col] : xmin;
        xmax = (xmax < data[i][col]) ? data[i][col] : xmax;
    }

    h = new TH1F("h", label[col], 100, xmin - 1, xmax + 1);
    for (int i = 0; i < (int)data.size(); i++)
        h->Fill(data[i][col]);

    return h;
}

//______________________________________________
// 6. get histogram
void MyReadData::SetBranchName(int id, TString bname)
{
    if (id < 0 || id >= ncol)
        return;

    label[id] = bname;
}

void MyReadData::SaveToRootFile(TString rname)
{
    TFile f(rname, "recreate");
    TTree t1("t1", "tree");

    double *val = new double[ncol];
    for (int i = 0; i < ncol; i++)
        t1.Branch(label[i], &val[i], label[i] + "/D");

    for (int n = 0; n < nrec; n++)
    {
        for (int i = 0; i < ncol; i++)
            val[i] = data[n][i];
        t1.Fill();
    }

    t1.Write();
    delete []val;
}
