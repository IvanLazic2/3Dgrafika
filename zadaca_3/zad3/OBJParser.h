#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>
#include <sstream>

using namespace std;

class OBJParser
{
    

public:
    using vType = vector<tuple<double, double, double>>;
    using tType = vector<tuple<double, double>>;

    class Face
    {
    public:
        vType V;
        vType VN;
        tType VT;
        vector<size_t> VIndexes;
        vector<size_t> VNIndexes;
        vector<size_t> VTIndexes;

        Face(
            vType v,
            vType vn,
            tType vt,
            vector<size_t> vi,
            vector<size_t> vni,
            vector<size_t> vti)
        {
            V = v;
            VN = vn;
            VT = vt;
            VIndexes = vi;
            VNIndexes = vni;
            VTIndexes = vti;
        }
    };

private:
    ifstream objfile;
    vType _v;
    vType _vn;
    tType _vt;
    vector<OBJParser::Face> _f;

    void set_values(vType &vtype)
    {
        vector<double> values(3);
        for (size_t i = 0; i < 3; i++)
        {
            double value;
            objfile >> value;

            values[i] = value;
        }
        vtype.emplace_back(values[0], values[1], values[2]);
    }

    void set_texture_values(tType &ttype)
    {
        vector<double> values(2);
        for (size_t i = 0; i < 2; i++)
        {
            double value;
            objfile >> value;

            values[i] = value;
        }
        ttype.emplace_back(values[0], values[1]);
    }

    void read()
    {
        if (objfile.is_open())
        {
            while (objfile.good())
            {
                string type;
                objfile >> type;

                if (type == "v")
                    set_values(_v);
                else if (type == "vn")
                    set_values(_vn);
                else if (type == "vt")
                    set_texture_values(_vt);
                else if (type == "f")
                {
                    vType fV;
                    vType fVN;
                    tType fVT;
                    vector<size_t> fVIndexes;
                    vector<size_t> fVNIndexes;
                    vector<size_t> fVTIndexes;

                    for (size_t i = 0; i < 3; i++)
                    {

                        string value;
                        objfile >> value;
                        stringstream ss(value);
                        string segment;
                        vector<int> seglist;

                        while (getline(ss, segment, '/'))
                        {
                            if (segment != "")
                            {
                                seglist.push_back(stoi(segment) - 1);
                            }
                        }

                        fV.push_back(_v[seglist[0]]);
                        fVN.push_back(_vn[seglist[seglist.size() - 1]]);
                        fVIndexes.push_back(seglist[0]);
                        fVNIndexes.push_back(seglist[seglist.size() - 1]);
                        if (seglist.size() == 3)
                        {
                            fVT.push_back(_vt[seglist[1]]);
                            fVTIndexes.push_back(seglist[1]);
                        }
                    }

                    _f.emplace_back(fV, fVN, fVT, fVIndexes, fVNIndexes, fVTIndexes);
                }
            }
        }
    }

public:
    OBJParser(const string filename)
    {
        objfile.open(filename);

        this->read();
    }

    vType get_v() const { return _v; }
    vType get_vn() const { return _vn; }
    tType get_vt() const { return _vt; }
    vector<OBJParser::Face> get_f() const { return _f; }
};