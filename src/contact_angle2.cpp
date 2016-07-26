#include "yocto/program.hpp"
#include "yocto/gfx/image/png.hpp"
#include "yocto/gfx/image/jpeg.hpp"
#include "yocto/gfx/image/tiff.hpp"
#include "yocto/gfx/ops/edges.hpp"
#include "yocto/gfx/ops/particles.hpp"
#include "yocto/gfx/ops/histogram.hpp"
#include "yocto/gfx/ops/filter.hpp"

#include "yocto/ios/ocstream.hpp"
#include "yocto/gfx/rawpix.hpp"

#include "yocto/math/alg/shapes2d.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/math/stat/descr.hpp"
#include "yocto/string/conv.hpp"

#include "yocto/math/alg/shapes2d.hpp"

using namespace yocto;
using namespace gfx;
using namespace math;

static unit_t y_limit = 0;
static inline bool is_bad_vertex(const vertex &v) throw()
{
    return v.y>=y_limit;
}

static inline void load_fit( FitShape2D<double> &fit_shape, const vlist &L, const pixmap<float> &wksp)
{
    for(const vnode *node=L.head;node;node=node->next)
    {
        const double x = node->vtx.x;
        const double y = node->vtx.y;
        const double w = wksp[node->vtx];
        fit_shape.append(x,y,w);
    }
}

static inline RGB Float2RGB(const float f)
{
    const uint8_t u = gist::float2byte(f);
    return RGB(u,u,u);
}

YOCTO_PROGRAM_START()
{
    YOCTO_GFX_DECL_FORMAT(jpeg);
    YOCTO_GFX_DECL_FORMAT(png);
    YOCTO_GFX_DECL_FORMAT(tiff);

    imageIO          &IMG = image::instance();

    if(argv[1]>0)
    {
        const string  filename = argv[1];
        const pixmapf original( IMG.loadf(filename,0) );
        const unit_t  w = original.w;
        const unit_t  h = original.h;
        y_limit         = h/2;
        xpatches      xps(original, new threading::engine(4,0,true) );
        Filter<float> F(w,h);
        size_t        nf = 0;
        if(argc>2)
        {
            nf = strconv::to<size_t>(argv[2],"#filters");
        }

        //______________________________________________________________________
        //
        //
        // extracting foreground...
        //
        //______________________________________________________________________
        pixmapf  img(original);
        IMG.save("img.png",img,NULL);

        std::cerr << "-- Filtering #" << nf << std::endl;
        for(size_t i=1;i<=nf;++i)
        {
            (std::cerr << ".").flush();
            F.Median(img,xps);
        }
        if(nf>0) std::cerr << std::endl;
        IMG.save("img_f.png",img,NULL);

        std::cerr << "-- Building Edges" << std::endl;
        Edges edges(w,h);
        edges.build_from(img,xps);
        IMG.save("edges.png",edges,NULL);


        std::cerr << "-- Extracting Foreground" << std::endl;
        pixmap1 edgesU(edges,gist::float2byte,edges);
        pixmap1 fg(w,h);
        separate(threshold::keep_foreground,fg,edgesU,xps);
        IMG.save("fg.png",fg,0);


        std::cerr << "-- Extracting Tags" << std::endl;
        tagmap tags(w,h);
        tags.build(fg,8);
        IMG.save("tags.png", tags, tags.colors, NULL);

        std::cerr << "-- Extract Particles and remove Shallows..." << std::endl;
        particles pa;
        pa.load(tags);
        std::cerr << "#particles=" << pa.size() << std::endl;
        pa.remove_shallow_with(tags);
        std::cerr << "#particles=" << pa.size() << std::endl;
        IMG.save("tags_lw.png", tags, tags.colors, NULL);

        std::cerr << "-- Removing Vertices" << std::endl;
        pa.reject_all_vertices_from(tags,is_bad_vertex);
        IMG.save("tags_ok.png",tags,tags.colors,NULL);
        std::cerr << "#particles=" << pa.size() << std::endl;

        if(pa.size()<2)
        {
            throw exception("Must have at least 2 particles...");
        }
        particle &p1 = *pa[1];
        particle &p2 = *pa[2];

        pixmapf wksp(w,h);
        p1.transfer(wksp,pixel<float>::invert,original);
        p2.transfer(wksp,pixel<float>::invert,original);
        IMG.save("wksp.png",wksp,NULL);

        std::cerr << "p1.size=" << p1.size << std::endl;
        std::cerr << "p2.size=" << p2.size << std::endl;

        FitConic<double> fit_conic;
        fit_conic.data.free();
        fit_conic.data.ensure(p1.size+p2.size);

        load_fit(fit_conic,p1,wksp);
        load_fit(fit_conic,p2,wksp);

        vector<double> A(6);
        fit_conic.compute(FitConicEllipse,A);
        std::cerr << "A=" << A << std::endl;

        point2d<double> center, radius;
        matrix<double>  rotate(2);
        fit_conic.Reduce(center, radius, rotate, A);
        std::cerr << "center=" << center << std::endl;
        std::cerr << "radius=" << radius << std::endl;
        std::cerr << "rotate=" << rotate << std::endl;

        pixmap3 wfit(w,h);
        p1.transfer(wfit,Float2RGB,original);
        p2.transfer(wfit,Float2RGB,original);

        const double vc  = FitConic<double>::Eval(A,center.x,center.y);
        RGB          col = named_color::fetch(YGFX_RED);
        for(unit_t j=0;j<h;++j)
        {
            for(unit_t i=0;i<w;++i)
            {
                const double val = FitConic<double>::Eval(A,i,j);
                if(val*vc>=0)
                {
                    wfit[j][i] = pixel<RGB>::blend(wfit[j][i], col, 80);
                }
            }
        }

        IMG.save("wfit.png",wfit,NULL);




#if 0
        p1.mask(wksp, 1.0f, 255);
        p2.mask(wksp, 1.0f, 255);
        IMG.save("wsym.png",wksp,NULL);

        p1.split_using(tags);
        p2.split_using(tags);
        pixmapf outl(w,h);
        p1.mask_border(outl, 1.0f, 255);
        p2.mask_border(outl, 1.0f, 255);

        IMG.save("outline.png",outl,NULL);
#endif
        

    }


}
YOCTO_PROGRAM_END()
