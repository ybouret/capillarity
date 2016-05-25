#include "yocto/program.hpp"
#include "yocto/graphics/rawpix.hpp"
#include "yocto/graphics/ops/gradient.hpp"
#include "yocto/graphics/ops/histogram.hpp"
#include "yocto/graphics/image/png.hpp"
#include "yocto/graphics/image/jpeg.hpp"
#include "yocto/graphics/image/tiff.hpp"
#include "yocto/ios/ocstream.hpp"


using namespace yocto;
using namespace graphics;

YOCTO_PROGRAM_START()
{
    image &IMG = image::instance();

    IMG.declare( new png_format()  );
    IMG.declare( new jpeg_format() );
    IMG.declare( new tiff_format() );

    const image::format &PNG = IMG["PNG"];

    if(argc>1)
    {
        const string filename = argv[1];
        pixmapf      pxm( IMG.loadf(filename,NULL) );
        const unit_t w = pxm.w;
        const unit_t h = pxm.h;

        xpatches     xps;
        xpatch::create(xps, pxm, NULL);

        pixmapf grd(w,h);
        {
            pixmapf tmp(w,h);
            gradient G;
            G.compute(grd, tmp, pxm, xps, NULL);
        }
        PNG.save("lens.png", pxm, NULL);
        PNG.save("grad.png", grd, NULL);

        pixmapf edge(w,h);
        {
            histogram H;
            H.update(grd, xps, NULL);
            const size_t t = H.threshold();
            std::cerr << "threshold=" << t << std::endl;
            threshold::apply(edge, t, grd, threshold::keep_foreground);
        }
        PNG.save("edge.png",edge,NULL);

        string output = filename;
        vfs::change_extension(output, "coords");
        std::cerr << "saving to " << output << std::endl;
        ios::wcstream fp(output);
        vector<float>  I(w,as_capacity);
        vector<unit_t> J(w,as_capacity);
        for(unit_t i=0;i<w;++i)
        {
            float   Imax = edge[0][i];
            unit_t  jmax = 0;
            for(unit_t j=1;j<h;++j)
            {
                const float Itmp = edge[j][i];
                if(Itmp>Imax)
                {
                    jmax = j;
                    Imax = Itmp;
                }
            }
            std::cerr << "Imax" << i << "=" << Imax << "@" << jmax << std::endl;
            I.push_back(Imax);
            J.push_back(jmax);
        }

    }

}
YOCTO_PROGRAM_END()
