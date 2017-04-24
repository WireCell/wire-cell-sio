#ifndef WIRECELLSIO_JSONDEPOSOURCE
#define WIRECELLSIO_JSONDEPOSOURCE

#include "WireCellIface/IDepoSource.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell {
    namespace Sio {
        class JsonDepoSource : public IDepoSource, public IConfigurable {
        public:
        
            JsonDepoSource();
            virtual ~JsonDepoSource();

            /// IDepoSource
            virtual bool operator()(IDepo::pointer& out);

            /// IConfigurable
            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);

        private:
            WireCell::IDepo::vector m_depos;
            bool m_eos;


        };
    }
}
#endif

