#include "yocto/program.hpp"
#include "yocto/gfx/image/png.hpp"
#include "yocto/gfx/image/jpeg.hpp"
#include "yocto/gfx/image/tiff.hpp"
#include "yocto/gfx/ops/edges.hpp"
#include "yocto/gfx/ops/particles.hpp"
#include "yocto/gfx/ops/histogram.hpp"
#include "yocto/gfx/ops/filter.hpp"
#include "yocto/gfx/ops/transform.hpp"
#include "yocto/gfx/ops/blur.hpp"
#include "yocto/gfx/draw/line.hpp"

#include "yocto/ios/ocstream.hpp"
#include "yocto/gfx/rawpix.hpp"

#include "yocto/math/alg/shapes2d.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/math/stat/descr.hpp"
#include "yocto/string/conv.hpp"

#include "yocto/math/alg/shapes2d.hpp"
#include "yocto/math/opt/cgrad.hpp"
#include "yocto/sort/remove-if.hpp"

using namespace yocto;
using namespace gfx;
using namespace math;

static unit_t y_limit = 0;
static inline bool is_bad_vertex(const vertex &v) throw()
{
    return v.y>=y_limit;
}

static inline
void remove_particle(particle::ptr &p, tagmap &tags)
{
    for(const vnode *node = p->head;node;node=node->next)
    {
        tags[node->vtx] = 0;
    }
    p->clear();
}

static inline float invert(const float x) throw() { return 1.0f-x; }

static inline void zero_above_y_limit(pixmapf &px)
{
    for(unit_t j=y_limit;j<px.h;++j)
    {
        for(unit_t i=0;i<px.w;++i)
        {
            px[j][i] = 0;
        }
    }
}

static inline void rescale_all( float *arr, const size_t n)
{
    assert(arr);
    assert(n>0);
    float amax = arr[0];
    float amin = arr[0];
    for(size_t i=1;i<n;++i)
    {
        const float atmp = arr[i];
        amax = max_of(amax,atmp);
        amin = min_of(amin,atmp);
    }
    const float adel = amax-amin;
    if(adel>0)
    {
        for(size_t i=0;i<n;++i)
        {
            arr[i] = (arr[i]-amin)/adel;
        }
    }
}

static inline RGB Float2RGB(const float f)
{
    const uint8_t u = gist::float2byte(f);
    return RGB(u,u,u);
}

class Shaper
{
public:
    point2d<double> center;
    point2d<double> radius;
    matrix<double>  rotate;

    Shaper(const size_t n) : center(), radius(), rotate(2) {}

    ~Shaper() throw() {}



private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Shaper);
};

static inline void load_fit( FitShape2D<double> &f, const particle &p , const pixmapf &w)
{
    for(const vnode *node = p.head;node;node=node->next)
    {
        const double x = node->vtx.x;
        const double y = node->vtx.y;
        f.append(x,y,w[node->vtx]);
    }
}


YOCTO_PROGRAM_START()
{
    YOCTO_GFX_DECL_FORMAT(jpeg);
    YOCTO_GFX_DECL_FORMAT(png);
    YOCTO_GFX_DECL_FORMAT(tiff);

    imageIO          &IMG = image::instance();
    thread_server     srv( new threading::engine(4,0,true) );
    if(argv[1]>0)
    {
        // Load original image
        std::cerr << "-- Loading" << std::endl;
        const string  filename = argv[1];
        const pixmapf original( IMG.loadf(filename,0) );
        const unit_t  w = original.w;
        const unit_t  h = original.h;
        y_limit = h/2;
        xpatches      xps(original,srv);
        transform     TR;
        Filter<float> F(w,h);
        Edges         edges(w,h);
        Histogram     H;

        // copy
        pixmapf img(original);
        IMG.save("img.png",img,NULL);

        //----------------------------------------------------------------------
        std::cerr << "-- Edges" << std::endl;
        //----------------------------------------------------------------------
        edges.build_from(img,xps);
        IMG.save("img_edges0.png",edges,NULL);

        //----------------------------------------------------------------------
        std::cerr << "-- Histogram..." << std::endl;
        //----------------------------------------------------------------------
        H.reset();
        H.update(edges,xps);
        const size_t ht = H.threshold();
        std::cerr << "ht=" << ht << std::endl;

        //----------------------------------------------------------------------
        std::cerr << "-- Prepare" << std::endl;
        //----------------------------------------------------------------------
        threshold::apply(img,ht,edges,threshold::keep_foreground);
        IMG.save("img_edges1.png",img,NULL);

        //----------------------------------------------------------------------
        std::cerr << "-- Build Tags" << std::endl;
        //----------------------------------------------------------------------
        tagmap tags(w,h);
        tags.build(img,8);

        //----------------------------------------------------------------------
        std::cerr << "-- Load Particles v0" << std::endl;
        //----------------------------------------------------------------------
        particles pa;
        pa.load(tags);
        pa.remove_shallow_with(tags);
        IMG.save("tags0.png",tags,tags.colors,NULL);
        if(pa.size()<2)
            throw exception("No enough particles");

        //----------------------------------------------------------------------
        std::cerr << "-- Load Particles v1 (remove too small)" << std::endl;
        //----------------------------------------------------------------------
        const size_t min_count = pa[1]->size/2;
        while(pa.size()>2)
        {
            const size_t n = pa.size();
            if(pa[n]->size>min_count)
            {
                break;
            }
            remove_particle(pa[n],tags);
            pa.pop_back();
        }
        IMG.save("tags1.png",tags,tags.colors,NULL);

        //----------------------------------------------------------------------
        std::cerr << "-- Load Particles v2 (cut upper half)" << std::endl;
        //----------------------------------------------------------------------
        pa.reject_all_vertices_from(tags,is_bad_vertex);
        IMG.save("tags2.png",tags,tags.colors,NULL);

        //----------------------------------------------------------------------
        std::cerr << "-- Load Particles v3 (keep 2 biggests)" << std::endl;
        //----------------------------------------------------------------------
        while(pa.size()>2)
        {
            remove_particle( pa.back(), tags);
            pa.pop_back();
        }
        IMG.save("tags3.png",tags,tags.colors,NULL);

        particle &p1 = *pa[1];
        particle &p2 = *pa[2];

        //----------------------------------------------------------------------
        std::cerr << "-- Building Workspace from original data" << std::endl;
        //----------------------------------------------------------------------
        pixmapf wksp(w,h);
        p1.transfer(wksp,invert,original);
        p2.transfer(wksp,invert,original);
        IMG.save("wksp.png",wksp,NULL);

        const patch a1 = p1.aabb();
        const patch a2 = p2.aabb();
        std::cerr << "a1=" << a1 << std::endl;
        std::cerr << "a2=" << a2 << std::endl;

        pixmap3 wimg(wksp,Float2RGB,wksp);

        FitConic<double> fit_shape;
        fit_shape.data.ensure(p2.size+p1.size);
        load_fit(fit_shape, p1, wksp);
        load_fit(fit_shape, p2, wksp);
        vector<double> A(6);
        fit_shape.compute(FitConicEllipse,A);
        std::cerr << "A=" << A << std::endl;
        Shaper S(fit_shape.data.size());
        FitConic<double>::Reduce(S.center, S.radius, S.rotate, A);

        std::cerr << "S.center=" << S.center << std::endl;
        std::cerr << "S.radius=" << S.radius << std::endl;
        std::cerr << "S.rotate=" << S.rotate << std::endl;

        point2d<double> u0,u1,v0,v1;
        u0.x = -S.radius.x;
        u0.y = 0;
        tao::mul(v0,S.rotate,u0);
        u1.x = S.radius.x;
        u1.y = 0;
        tao::mul(v1,S.rotate,u1);

        v0 += S.center;
        v1 += S.center;

        draw_line(wimg, v0.x, v0.y, v1.x, v1.y, named_color::fetch(YGFX_RED));

        u0.x=0;u1.x=0;
        u0.y=-S.radius.y;
        u1.y= S.radius.y;
        tao::mul(v0,S.rotate,u0);
        tao::mul(v1,S.rotate,u1);
        v0 += S.center;
        v1 += S.center;

        draw_line(wimg, v0.x, v0.y, v1.x, v1.y, named_color::fetch(YGFX_BLUE));

        IMG.save("wimg.png",wimg,NULL);

#if 0
        y_limit = min_of(a1.lower.y,a2.lower.y);
        const unit_t x_lower = min_of(a1.lower.x,a2.lower.x);
        const unit_t x_upper = max_of(a1.upper.x,a2.upper.x);
        const patch  cutp( vertex(x_lower,0), vertex(x_upper,y_limit-1) );

        //----------------------------------------------------------------------
        std::cerr << "-- Processing with y limit=" << y_limit << std::endl;
        //----------------------------------------------------------------------
        const unit_t  cw = cutp.width.x;
        const unit_t  ch = cutp.width.y;
        pixmapf       cimg(cw,ch);
        xpatches      cxps(cimg,srv);
        Filter<float> cF(cw,ch);

        for(unit_t j=0;j<ch;++j)
        {
            for(unit_t i=0;i<cw;++i)
            {
                cimg[j][i] = original[j][x_lower+i];
            }
        }
        IMG.save("img_cut0.png",cimg,0);
        rescale_all((float*)cimg.entry, cimg.items);
        IMG.save("img_cut1.png",cimg,0);

        cF.Close(cimg,cxps);
        IMG.save("img_cut2.png",cimg,0);
        cF.Close(cimg,cxps);
        IMG.save("img_cut3.png",cimg,0);


        Edges cedges(cw,ch);
        cedges.build_from(cimg,cxps);
        IMG.save("edges_sub0.png",cedges,NULL);
        cF.Open(cedges,cxps);
        IMG.save("edges_sub1.png",cedges,NULL);
        cF.PseudoMedian(cedges,cxps);
        IMG.save("edges_sub2.png",cedges,NULL);

        //separate(threshold::keep_foreground,cimg,cedges,cxps);
        rescale_all((float*)cedges.entry, cedges.items);
        IMG.save("edges_sub3.png",cedges,NULL);
        cF.Close(cedges,cxps);
        IMG.save("edges_sub4.png",cedges,NULL);
#endif


    }

}
YOCTO_PROGRAM_END()