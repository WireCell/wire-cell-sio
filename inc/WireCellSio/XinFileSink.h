/** A "Xin" file is one used in the WC prototype and some large tests in the toolkit. */

#ifndef WIRECELLSIO_XINFILESINK
#define WIRECELLSIO_XINFILESINK

#include "WireCellIface/IFrameSink.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell {
    namespace Sio {

        class XinFileSink : public IFrameSink, public IConfigurable {
        public:

            XinFileSink();
            virtual ~XinFileSink();

            /// IFrameSink
            virtual bool operator()(const IFrame::pointer& out);

            /// IConfigurable
            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);
        private:
            Configuration m_cfg;

        };
    }
}

#endif
