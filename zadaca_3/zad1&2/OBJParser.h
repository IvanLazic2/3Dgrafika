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

    class Face
    {
    public:
        vType V;
        vType VN;
        vector<size_t> VIndexes;
        vector<size_t> VNIndexes;

        Face(
            vType v,
            vType vn,
            vector<size_t> vi,
            vector<size_t> vni)
        {
            V = v;
            VN = vn;
            VIndexes = vi;
            VNIndexes = vni;
        }
    };

private:
    ifstream objfile;
    vType _v;
    vType _vn;
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
        vtype.push_back({values[0], values[1], values[2]});
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
                else if (type == "f")
                {
                    vType fV;
                    vType fVN;
                    vector<size_t> fVIndexes;
                    vector<size_t> fVNIndexes;

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
                    }

                    _f.push_back({fV, fVN, fVIndexes, fVNIndexes});
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
    vector<OBJParser::Face> get_f() const { return _f; }
};