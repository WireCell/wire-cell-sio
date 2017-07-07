/** This only works for MB. It's used by tests and XinFileSource. */

#ifndef WIRECELLSIO_XINFILEITERATOR
#define WIRECELLSIO_XINFILEITERATOR

#include "WireCellUtil/Waveform.h"
#include "WireCellIface/IFrame.h"

#include <vector>

class TH2;
class TFile;

namespace WireCell {
    namespace Sio {


        class XinFileIterator {
            TH2* hist[3];		// per plane
            WireCell::Waveform::ChannelMaskMap ret;
            TFile *file;
        public:
            XinFileIterator(const char* filename, const char* histtype="raw");

            int plane(int ch);

            int index(int ch);

            std::vector<float> at(int ch);

            
            void clear();
  
            IFrame::pointer frame();
    
        };
    }
}

#endif



