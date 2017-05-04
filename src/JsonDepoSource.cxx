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
    cfg["qunit"] = 1.0;         // multiply this to each "q" value
    cfg["recombination"] = "";  // model to produce electrons from "q" value
                                // Each model may have its own parameters.
                                // See the structs below.
    return cfg;
}

class IRecombinationModel { // Fixme: should promote this to first class model!
public:
    virtual ~IRecombinationModel() {}
    virtual double operator()(Json::Value depo) = 0;
};
class IdentityModel : public IRecombinationModel {
    double qunit;
public:
    IdentityModel(Json::Value cfg) {
        qunit = cfg["qunit"].asDouble();
    }
    virtual ~IdentityModel() {}
    virtual double operator()(Json::Value depo) {
        return depo["q"].asDouble();
    }
};
struct BirksModel : public IRecombinationModel {
    double qunit;
    double A3t, k3tepsilon;
public:
    // see http://lar.bnl.gov/properties/pass.html#recombination
    // A3t is unitless, k3tepsilon is in units of [distance]/[energy].
    BirksModel(Json::Value cfg) {
        qunit = cfg["qunit"].asDouble();
        A3t = cfg["A3t"].asDouble();
        k3tepsilon = cfg["k3tepsilon"].asDouble();
    }
    virtual ~BirksModel() {}
    virtual double operator()(Json::Value depo) {
        const double dE = qunit * depo["q"].asDouble();
        const double dX = depo["s"].asDouble();
        return dE * A3t / (1 + k3tepsilon * dE/dX);
    }
};    
struct BoxModel : public IRecombinationModel {
    double qunit;
    double A, Boverepsilon;
public:
    // see http://lar.bnl.gov/properties/pass.html#recombination
    // A is unitless and Boverepsilon is in units of [distance]/[energy]
    BoxModel(Json::Value cfg) {
        qunit = cfg["qunit"].asDouble();
        A = cfg["A"].asDouble();
        Boverepsilon = cfg["Boverepsilon"].asDouble();
    }
    virtual ~BoxModel() {}
    virtual double operator()(Json::Value depo) {
        const double dE = qunit * depo["q"].asDouble();
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
    std::string model = get<string>(cfg, "recombination");
    {
        std::locale loc;
        for (std::string::size_type ind = 0; ind<model.size(); ++ind) {
            model[ind] = std::tolower(model[ind], loc);
        }
    }
    if (model.empty() || model == "identity") {
        qtoelectrons = new IdentityModel(cfg);
    }
    if (model == "birks") {
        qtoelectrons = new BirksModel(cfg);
    }
    if (model == "box" || model == "modifiedbox") {
        qtoelectrons = new BoxModel(cfg);
    }
    if (!qtoelectrons) {
        cerr << "JsonDepoSource::configure: unknown model for producing deposited electrons: \"" << model << "\"\n";
        return;
    }


    // get and load JSON file.
    string filename = get<string>(cfg,"filename");
    string dotpath = get<string>(cfg,"jsonpath","depos");
    if (filename.empty()) {
        cerr << "JsonDepoSource::configure: no JSON filename given" << endl;
        return;                 // fixme: uh, error handle much?
    }
    Json::Value top = WireCell::Persist::load(filename.c_str());


    auto jdepos = branch(top, dotpath);
    for (auto jdepo : jdepos) {
        auto depo = std::make_shared<SimpleDepo>(
            get(jdepo,"t",0.0*units::ns),
            Point(get(jdepo, "x", 0.0),
                  get(jdepo, "y", 0.0),
                  get(jdepo, "z", 0.0)),
            (*qtoelectrons)(jdepo));
        m_depos.push_back(depo);
    }
    std::sort(m_depos.begin(), m_depos.end(), descending_time);
    delete qtoelectrons;
    qtoelectrons = nullptr;
}


