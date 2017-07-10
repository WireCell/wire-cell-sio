/** Sink data to a file format used by the "Magnify" GUI . 
 *
 * FIXME: currently this class TOTALLY violates the encapsulation of
 * DFP by requiring the input file in order to transfer input data out
 * of band of the flow.
 */

#ifndef WIRECELLSIO_MAGNIFYFILESINK
#define WIRECELLSIO_MAGNIFYFILESINK

#include "WireCellIface/IFrameSink.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell {
    namespace Sio {

        class MagnifySink : public IFrameSink, public IConfigurable {
        public:

            MagnifySink();
            virtual ~MagnifySink();

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
