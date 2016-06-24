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
#include "yocto/math/stat/descr.hpp"

using namespace yocto;
using namespace graphics;

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
        Histogram H;
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




    math::fit_circle<double> FitCircle;
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

    // approximate center
    double xd = (xmin+xmax)/2;
    double yd = (h+(circle_radius+circle_center.y))/2;
    double rd = ((xmax-xmin)/2)*0.6;
    double rd2 = rd*rd;


    for(unit_t x=0;x<w;++x)
    {
        const double dx  = x-circle_center.x;
        const double dx2 = dx*dx;

        const double DX  = x - xd;
        const double DX2 = DX*DX;

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
                const double DY  = y-yd;
                const double DY2 = DY*DY;
                if(x>=xmin&&x<=xmax)
                {

                    drop[y][x] = grd[y][x];
                    //std::cerr << "(" << x << "," << y << ")" << std::endl;
                    if(DX2+DY2<=rd2)
                    {
                        //std::cerr << "*" << std::endl;
                        drop[y][x] = 0.0f;
                    }
                    //drop[0][x] = 1.0f;
                }
            }
        }

    }

    std::cerr << "xd=" << xd << std::endl;
    std::cerr << "yd=" << yd << std::endl;
    std::cerr << "rd=" << rd << std::endl;


    PNG.save("wksp.png",wksp,NULL);
    PNG.save("drop.png",drop,NULL);

    // extracting drop

    pixmapf drop_edges(w,h);
    {
        std::cerr << "-- isolating edges..." << std::endl;
        Histogram H;
        H.update(drop,xps,&server);
        const size_t t = H.threshold();
        std::cerr << "drop_threshold=" << t << std::endl;
        threshold::apply(drop_edges, t, drop, threshold::keep_foreground);
    }
    PNG.save("drop_edges.png",drop_edges,NULL);

    B.build(drop_edges,8);
    PNG.save("drop_blobs.png", B, blobColors, NULL);

    pixmapf drop_edge(w,h);
    B.content[1]->transfer(drop_edge,drop_edges);
    if(B.content.size()>=2)
    {
        B.content[2]->transfer(drop_edge,drop_edges);
    }
    PNG.save("drop_edge.png",drop_edge,NULL);


    pixmapf shape(w,h);
    unit_t sy_min = h-1;
    for(unit_t j=h-1;j>=0;--j)
    {
        for(unit_t i=xd;i>=0;--i)
        {
            if(drop_edge[j][i]>0)
            {
                shape[j][i] = 1.0;
                sy_min=j;
                break;
            }
        }
        for(unit_t i=xd;i<w;++i)
        {
            if(drop_edge[j][i]>0)
            {
                shape[j][i] = 1.0;
                sy_min=j;
                break;
            }
        }
    }
    // remove sy_min
    for(unit_t i=0;i<w;++i)
    {
        shape[sy_min][i] = 0;
    }
    ++sy_min;
    PNG.save("shape.png",shape,NULL);
    std::cerr << "sy_min=" << sy_min << std::endl;

    vector<double> SX;
    vector<double> SY;
    const unit_t sy_max = (h+sy_min)/2;
    for(unit_t j=sy_min;j<sy_max;++j)
    {
        for(unit_t i=0;i<w;++i)
        {
            if(shape[j][i]>0)
            {
                SX.push_back(i);
                SY.push_back(j);
            }
        }
    }

    {
        pixmap3 shape3(wksp);
        
        for(size_t i=1;i<=SX.size();++i)
        {
            shape3[unit_t(SY[i])][unit_t(SX[i])] = RGB(0,255,0);
        }
        PNG.save("shape3.png",shape3,NULL);

    }

    {
        ios::wcstream fp("shape.dat");
        for(size_t i=1;i<=SX.size();++i)
        {
            fp("%g %g\n", SX[i], SY[i]);
        }
    }
    
    math::compute_average(xd,SX);
    std::cerr << "xd=" << xd << std::endl;
}
YOCTO_PROGRAM_END()
