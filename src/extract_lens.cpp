#include "yocto/program.hpp"
#include "yocto/gfx/rawpix.hpp"
#include "yocto/gfx/ops/edges.hpp"
#include "yocto/gfx/ops/histogram.hpp"
#include "yocto/gfx/ops/particles.hpp"
#include "yocto/gfx/ops/filter.hpp"
#include "yocto/gfx/image/png.hpp"
#include "yocto/gfx/image/jpeg.hpp"
#include "yocto/gfx/image/tiff.hpp"
#include "yocto/ios/ocstream.hpp"

#include "yocto/math/alg/shapes2d.hpp"


using namespace yocto;
using namespace gfx;
using namespace math;

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

        const unit_t  w = pxm.w;
        const unit_t  h = pxm.h;
        Filter<float> F(w,h);

        xpatches     xps(pxm,true);

        F.PseudoMedian(pxm,xps);
        PNG.save("lens2.png",pxm,NULL);


        Edges edges(w,h);
        edges.build_from(pxm,xps);

        PNG.save("edges.png",edges,0);

        pixmapf fg(w,h);
        separate(threshold::keep_foreground,fg,edges,xps);
        PNG.save("fg.png",fg,0);

        tagmap tags(w,h);
        tags.build(fg,8);
        tags.colors.shift = YGFX_WHITE;
        PNG.save("tags.png",tags,tags.colors,0);



        particles pa;
        pa.load(tags);
        if(pa.size()<=0)
            throw exception("No main edge detected!");

        pixmapf edge(w,h);
        pa[1]->transfer(edge,pxm);
        PNG.save("edge.png",edge,0);

        F.Close(edge,xps);
        PNG.save("edge2.png",edge,0);

        string output = filename;
        vfs::change_extension(output, "coords");
        std::cerr << "saving to " << output << std::endl;
        ios::wcstream fp(output);
        FitConic<double>   fit_conic;
        FitCircle<double>  fit_circle;

        for(unit_t i=0;i<w;++i)
        {

            for(unit_t j=0;j<h;++j)
            {
                if(edge[j][i]>0)
                {
                    fp("%g %g\n", double(i), double(j));
                    edge[j][i] = 1.0;
                    fit_conic.append(double(i), double(j));
                    fit_circle.append(double(i), double(j));
                    break;
                }
            }

        }

        PNG.save("edge_fin.png",edge,NULL);

        double          circle_radius = 0;
        point2d<double> circle_center;
        fit_circle.compute(circle_center,circle_radius);
        std::cerr << "circle_radius=" << circle_radius << std::endl;
        std::cerr << "circle_center=" << circle_center << std::endl;

#if 0
        pixmapf grd(w,h);
        {
            pixmapf tmp(w,h);
            gradient G;
            G.compute(grd, tmp, pxm, xps, NULL);
        }
        PNG.save("grad.png", grd, NULL);


        pixmapf edges(w,h);
        {
            Histogram H;
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

        Filter<float> F;



        string output = filename;
        vfs::change_extension(output, "coords");
        std::cerr << "saving to " << output << std::endl;
        ios::wcstream fp(output);
        fit_conic<double>  FitConic;
        fit_circle<double> FitCircle;

        for(unit_t i=0;i<w;++i)
        {

            for(unit_t j=0;j<h;++j)
            {
                if(edge[j][i]>0)
                {
                    fp("%g %g\n", double(i), double(j));
                    edge[j][i] = 1.0;
                    FitConic.append(double(i), double(j));
                    FitCircle.append(double(i), double(j));
                    break;
                }
            }

        }

        PNG.save("edge1.png",edge,NULL);

        double          circle_radius = 0;
        point2d<double> circle_center;
        FitCircle.solve(circle_radius, circle_center);
        std::cerr << "circle_radius=" << circle_radius << std::endl;
        std::cerr << "circle_center=" << circle_center << std::endl;

#endif


    }

}
YOCTO_PROGRAM_END()
