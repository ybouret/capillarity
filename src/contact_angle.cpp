#include "yocto/program.hpp"
#include "yocto/graphics/image/png.hpp"
#include "yocto/graphics/image/jpeg.hpp"
#include "yocto/graphics/image/tiff.hpp"
#include "yocto/graphics/ops/gradient.hpp"
#include "yocto/graphics/ops/blobs.hpp"
#include "yocto/graphics/ops/histogram.hpp"
#include "yocto/graphics/ops/blur.hpp"
#include "yocto/graphics/ops/mix.hpp"

#include "yocto/ios/ocstream.hpp"
#include "yocto/graphics/rawpix.hpp"

#include "yocto/math/alg/shapes.hpp"
#include "yocto/string/conv.hpp"

using namespace yocto;
using namespace graphics;
using namespace math;

YOCTO_PROGRAM_START()
{
    const RGB C1 = named_color::get("indian_red");
    const RGB C2 = named_color::get("indigo");

    if(argc<3)
    {
        throw exception("usage: %s file xmin xmax", program);
    }

    const string filename = argv[1];
    const float  xmin     = strconv::to<int>(argv[2],"xmin");
    const unit_t xmax     = strconv::to<int>(argv[3],"xmax");
    std::cerr << "xmin=" << xmin << ", xmax=" << xmax << std::endl;

    image &IMG = image::instance();

    IMG.declare( new png_format()  );
    IMG.declare( new jpeg_format() );
    IMG.declare( new tiff_format() );

    const image::format &PNG = IMG["PNG"];

    threading::engine server(true);


    pixmapf      pxm( IMG.loadf(filename,NULL) );
    PNG.save("lens.png", pxm, NULL);


    const unit_t w = pxm.w;
    const unit_t h = pxm.h;

    xpatches     xps;
    xpatch::create(xps, pxm, &server);

    pixmapf src(w,h);
    src.copy(pxm);
    if(true)
    {
        std::cerr << "-- blurring" << std::endl;
        blur blr(0.5f);
        blr.apply(src, pxm, xps, &server);
    }

    pixmapf grd(w,h);
    {
        std::cerr << "-- gradient" << std::endl;
        pixmapf tmp(w,h);
        gradient G;
        G.compute(grd, tmp, src, xps, &server);
    }
    PNG.save("grad.png", grd, NULL);

    pixmapf grd2(grd);
    {
        for(unit_t j=0;j<h;++j)
        {
            for(unit_t i=xmin;i<=xmax;++i)
            {
                grd2[j][i] = 0;
            }
        }
    }
    PNG.save("grad2.png", grd2, NULL);


    pixmapf edges(w,h);
    {
        std::cerr << "-- isolating edges..." << std::endl;
        histogram H;
        H.update(grd,xps,&server);
        const size_t t = H.threshold();
        std::cerr << "edges_threshold=" << t << std::endl;
        threshold::apply(edges, t, grd2, threshold::keep_foreground);
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
    B.content[1]->transfer(edge,edges);
    if(B.content.size()>=2)
    {
        B.content[2]->transfer(edge,edges);
    }
    PNG.save("edge.png",edge,NULL);




    fit_circle<double> FitCircle;
    pixmap3 wksp(w,h);
    for(unit_t x=0;x<w;++x)
    {
        for(unit_t y=0;y<h;++y)
        {
            const float   g = pxm[y][x];
            const uint8_t u = gist::float2byte(g);
            wksp[y][x]      = RGB(u,u,u);
            if(edge[y][x]>0)
            {
                wksp[y][x] = C1;
                FitCircle.append(x,y);
            }
        }
    }

    double          circle_radius = 0;
    point2d<double> circle_center;
    FitCircle.solve(circle_radius, circle_center);
    std::cerr << "circle_radius=" << circle_radius << " @" << circle_center << std::endl;
    const double r2 = circle_radius*circle_radius;

    pixmapf drop(w,h);
    for(unit_t x=0;x<w;++x)
    {
        const double dx = x-circle_center.x;
        const double dx2 = dx*dx;
        for(unit_t y=0;y<h;++y)
        {
            const double dy  = y-circle_center.y;
            const double dy2 = dy*dy;
            if(dx2+dy2<=r2)
            {
                wksp[y][x] = mix::blend(wksp[y][x],C2,90);
            }
            else
            {
                if(x>=xmin&&x<=xmax)
                {
                    drop[y][x] = grd[y][x];
                }
            }
        }

    }

    PNG.save("wksp.png",wksp,NULL);
    PNG.save("drop.png",drop,NULL);

    // extracting drop
    pixmapf drop_edges(w,h);
    {
        std::cerr << "-- isolating edges..." << std::endl;
        histogram H;
        H.update(drop,xps,&server);
        const size_t t = H.threshold();
        std::cerr << "drop_threshold=" << t << std::endl;
        threshold::apply(drop_edges, t, drop, threshold::keep_foreground);
    }
    PNG.save("drop_edges.png",drop_edges,NULL);

    B.build(drop_edges,8);
    PNG.save("drop_blobs.png", B, blobColors, NULL);



#if 0


    std::cerr << "-- keeping background..." << std::endl;
    pixmapf bg(w,h);
    {
        histogram H;
        H.update(pxm,xps,&server);
        const size_t t = H.threshold();
        std::cerr << "background_threshold=" << t << std::endl;
        threshold::apply(bg, t, pxm, threshold::keep_background);
        PNG.save("bg.png",bg,NULL);
    }

    std::cerr << "-- detecting gradient..." << std::endl;
    pixmapf grd(w,h);
    {
        pixmapf tmp(w,h);
        gradient G;
        G.compute(grd, tmp, bg, xps, &server);
    }
    PNG.save("grad.png", grd, NULL);

    std::cerr << "-- isolating edges..." << std::endl;
    pixmapf edges(w,h);
    {
        histogram H;
        H.update(grd,xps,&server);
        const size_t t = H.threshold();
        std::cerr << "edges_threshold=" << t << std::endl;
        threshold::apply(edges, t, grd, threshold::keep_foreground);
    }
    PNG.save("edges.png",edges,NULL);

    std::cerr << "-- isolating lens..." << std::endl;
    pixmap3 wksp(w,h);
    for(unit_t x=0;x<w;++x)
    {
        for(unit_t y=0;y<h;++y)
        {
            const float   g = pxm[y][x];
            const uint8_t u = gist::float2byte(g);
            wksp[y][x]      = RGB(u,u,u);
            if(bg[y][x]>0)
            {
                wksp[y][x] = edge_color;
            }
        }
    }
    PNG.save("wksp.png",wksp,NULL);
#endif

}
YOCTO_PROGRAM_END()
