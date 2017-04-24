#include "WireCellSio/JsonDepoSource.h"
#include "WireCellIface/SimpleDepo.h"
#include "WireCellUtil/NamedFactory.h"

#include "WireCellUtil/Point.h"
#include "WireCellUtil/Persist.h"

WIRECELL_FACTORY(JsonDepoSource, WireCell::Sio::JsonDepoSource, WireCell::IDepoSource, WireCell::IConfigurable);

#include <iostream>
#include <string>

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
    return cfg;
}

void Sio::JsonDepoSource::configure(const WireCell::Configuration& cfg)
{
    m_depos.clear();

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
            get(jdepo,"q",1000.0));
        m_depos.push_back(depo);
    }
    std::sort(m_depos.begin(), m_depos.end(), descending_time);
}


