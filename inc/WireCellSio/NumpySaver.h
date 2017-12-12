/** Save frames to Numpy array files. */

#ifndef WIRECELLSIO_NUMPYSAVER
#define WIRECELLSIO_NUMPYSAVER

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IDepoFilter.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell {
    namespace Sio {
        class NumpySaver : public WireCell::IFrameFilter,
                           public WireCell::IDepoFilter,
                           public WireCell::IConfigurable {
        public:
            NumpySaver();
            virtual ~NumpySaver();

            /// IFrameFilter
            virtual bool operator()(const WireCell::IFrame::pointer& inframe,
                                    WireCell::IFrame::pointer& outframe);

            /// IDepoFilter.  This works by buffering depos and saving
            /// them at the same time a frame is saved.
            virtual bool operator()(const WireCell::IDepo::pointer& indepo,
                                    WireCell::IDepo::pointer& outdepo);

            /// IConfigurable
            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);
        private:

            Configuration m_cfg;
            int m_save_count;   // count frames saved

            std::vector<WireCell::IDepo::pointer> m_depos;
        };
    }
}
#endif
