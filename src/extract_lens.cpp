#include "yocto/program.hpp"
#include "yocto/graphics/rawpix.hpp"
#include "yocto/graphics/ops/gradient.hpp"
#include "yocto/graphics/ops/histogram.hpp"
#include "yocto/graphics/ops/blobs.hpp"
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
        PNG.save("lens.png", pxm, NULL);

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
        PNG.save("grad.png", grd, NULL);

        pixmapf edges(w,h);
        {
            histogram H;
            H.update(grd, xps, NULL);
            const size_t t = H.threshold();
            std::cerr << "threshold=" << t << std::endl;
            threshold::apply(edges, t, grd, threshold::keep_foreground);
        }
        PNG.save("edges.png",edges,NULL);

        blobs B(w,h);
        B.build(edges,8);
        std::cerr << "#blobs=" << B.content.size() << std::endl;
        get_named_color<size_t> blobColors;
        PNG.save("blobs.png", B, blobColors, NULL);
        if(B.content.size()<=0)
        {
            throw exception("No main edge detected!");
        }

        pixmapf edge(w,h);
        B.content[1]->transfer(edge, pxm);
        PNG.save("edge.png",edge,NULL);


        string output = filename;
        vfs::change_extension(output, "coords");
        std::cerr << "saving to " << output << std::endl;
        ios::wcstream fp(output);
        for(unit_t i=0;i<w;++i)
        {

            for(unit_t j=0;j<h;++j)
            {
                if(edge[j][i]>0)
                {
                    fp("%g %g\n", double(i), double(j));
                    edge[j][i] = 1.0;
                    break;
                }
            }

        }

        PNG.save("edge1.png",edge,NULL);



    }

}
YOCTO_PROGRAM_END()
