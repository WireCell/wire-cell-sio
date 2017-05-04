#include "WireCellSio/JsonDepoSource.h"
#include "WireCellIface/SimpleDepo.h"
#include "WireCellUtil/NamedFactory.h"

#include "WireCellUtil/Point.h"
#include "WireCellUtil/Persist.h"

WIRECELL_FACTORY(JsonDepoSource, WireCell::Sio::JsonDepoSource, WireCell::IDepoSource, WireCell::IConfigurable);

#include <iostream>
#include <string>
#include <locale>               // for std::tolower

using namespace std;
using namespace WireCell;

Sio::JsonDepoSource::JsonDepoSource()
    : m_eos(false)
{
}

Sio::JsonDepoSource::~JsonDepoSource()
{
}


bool Sio::JsonDepoSource::operator()(IDepo::pointer& out)
{
    if (m_eos) {
	return false;
    }

    if (m_depos.empty()) {
	m_eos = true;
	out = nullptr;
	return true;
    }

    out = m_depos.back();
    m_depos.pop_back();
    return true;
}

WireCell::Configuration Sio::JsonDepoSource::default_configuration() const
{
    Configuration cfg;
    cfg["filename"] = "";       // json file name
    cfg["jsonpath"] = "depos";  // path to depo list in json data
    cfg["model"] = "";    // model to produce electrons from "q" value
                                // Each model may have its own parameters.
                                // See the structs below.
    return cfg;
}

class IRecombinationModel { // Fixme: should promote this to first class model!
public:
    virtual ~IRecombinationModel() {}
    virtual double operator()(Json::Value depo) = 0;
};
class ElectronsModel : public IRecombinationModel {
    double scale;
public:
    ElectronsModel(Json::Value cfg) {
        scale = get(cfg,"scale",-1.0*units::eplus);
        cerr << "ElectronsModel: scale by " << scale << endl;
        cerr << cfg << endl;
    }
    virtual ~ElectronsModel() {}
    virtual double operator()(Json::Value depo) {
        return scale*depo["n"].asDouble();
    }
};
class ScaledModel : public IRecombinationModel {
    double scale;
public:
    ScaledModel(Json::Value cfg) {
        scale = get(cfg,"scale",1.0);
    }
    virtual ~ScaledModel() {}
    virtual double operator()(Json::Value depo) {
        return scale*depo["q"].asDouble();
    }
};

const double k3t_default = 0.0486*units::gram/(units::MeV*units::cm2)*(units::kilovolt/units::cm);
const double B_default = 0.212*units::gram/(units::MeV*units::cm2)*(units::kilovolt/units::cm);
const double Efield_default = 273*units::volt/units::cm;
const double density_default = 1.396*units::gram/units::mm3;

struct BirksModel : public IRecombinationModel {
    double scale;
    double A3t, k3tepsilon;
public:
    // see http://lar.bnl.gov/properties/pass.html#recombination
    // A3t is unitless, k3tepsilon is in units of [distance]/[energy].
    BirksModel(Json::Value cfg) {
        scale = get(cfg,"scale",1.0);
        A3t = get(cfg,"A3t",0.8);
        const double k3t = get(cfg,"k3t",k3t_default);
        const double Efield = get(cfg,"Efield", Efield_default);
        const double density = get(cfg,"density", density_default);
        k3tepsilon = k3t*Efield*density;
    }
    virtual ~BirksModel() {}
    virtual double operator()(Json::Value depo) {
        const double dE = scale * depo["q"].asDouble();
        const double dX = depo["s"].asDouble();
        return dE * A3t / (1 + k3tepsilon * dE/dX);
    }
};    
struct BoxModel : public IRecombinationModel {
    double scale;
    double A, Boverepsilon;
public:
    // see http://lar.bnl.gov/properties/pass.html#recombination
    // A is unitless and Boverepsilon is in units of [distance]/[energy]
    BoxModel(Json::Value cfg) {
        scale = get(cfg,"scale",1.0);
        A = cfg["A"].asDouble();
        const double B = get(cfg,"A",B_default);
        const double Efield = get(cfg,"Efield", Efield_default);
        const double density = get(cfg,"density", density_default);
        Boverepsilon = B/(Efield*density);
    }
    virtual ~BoxModel() {}
    virtual double operator()(Json::Value depo) {
        const double dE = scale * depo["q"].asDouble();
        const double dX = depo["s"].asDouble();
        const double val = Boverepsilon * dE/dX;
        return dE * log(A + val) / val;
    }
};        


void Sio::JsonDepoSource::configure(const WireCell::Configuration& cfg)
{
    m_depos.clear();

    // Figure out how to turn "q" and maybe also "s" into number of electrons
    IRecombinationModel* qtoelectrons = nullptr;
    std::string model = get<string>(cfg, "model");
    {
        std::locale loc;
        for (std::string::size_type ind = 0; ind<model.size(); ++ind) {
            model[ind] = std::tolower(model[ind], loc);
        }
    }
    if (model.empty() || model == "electrons") {
        qtoelectrons = new ElectronsModel(cfg);
    }
    if (model.empty() || model == "scaled") {
        qtoelectrons = new ScaledModel(cfg);
    }
    if (model == "birks") {
        qtoelectrons = new BirksModel(cfg);
    }
    if (model == "box" || model == "modifiedbox") {
        qtoelectrons = new BoxModel(cfg);
    }
    if (!qtoelectrons) {
        cerr << "Sio::JsonDepoSource::configure: unknown model for producing deposited electrons: \"" << model << "\"\n";
        return;
    }
    cerr << "Sio::JsonDepoSource::configure: model: \"" << model << "\"\n";

    // get and load JSON file.
    string filename = get<string>(cfg,"filename");
    string dotpath = get<string>(cfg,"jsonpath","depos");
    if (filename.empty()) {
        cerr << "JsonDepoSource::configure: no JSON filename given" << endl;
        return;                 // fixme: uh, error handle much?
    }
    Json::Value top = WireCell::Persist::load(filename.c_str());

    double qtot = 0;
    auto jdepos = branch(top, dotpath);
    for (auto jdepo : jdepos) {
        const double q = (*qtoelectrons)(jdepo);
        qtot += q;
        auto depo = std::make_shared<SimpleDepo>(
            get(jdepo,"t",0.0*units::ns),
            Point(get(jdepo, "x", 0.0),
                  get(jdepo, "y", 0.0),
                  get(jdepo, "z", 0.0)),
            q);
        m_depos.push_back(depo);
    }
    std::sort(m_depos.begin(), m_depos.end(), descending_time);
    cerr << "Sio::JsonDepoSource::configure: "
         << "slurped in " << m_depos.size() << " depositions, "
         << " = " << -1*qtot/units::eplus << " electrons\n";
    delete qtoelectrons;
    qtoelectrons = nullptr;
}


