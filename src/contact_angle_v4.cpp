#include "yocto/program.hpp"

#include "yocto/gfx/image.hpp"
#include "yocto/gfx/image/jpeg.hpp"
#include "yocto/gfx/image/png.hpp"
#include "yocto/gfx/image/tiff.hpp"

#include "yocto/gfx/ops/edges.hpp"
#include "yocto/gfx/ops/filter.hpp"
#include "yocto/string/conv.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/gfx/draw/line.hpp"
#include "yocto/gfx/draw/circle.hpp"
#include "yocto/math/alg/shapes2d.hpp"
#include "yocto/sort/quick.hpp"

using namespace yocto;
using namespace gfx;
using namespace math;


static unit_t y_limit = 0;
static inline bool is_bad_vertex(const vertex &v) throw()
{
    return v.y>=y_limit;
}

static inline
int compare_by_decreasing_particle_height( const particle::ptr &lhs, const particle::ptr &rhs )
{
    const vertex L = lhs->width();
    const vertex R = rhs->width();
    return __compare(R.y,L.y);
}

#include "yocto/associative/map.hpp"
static inline
int compare_vtx_by_x(const vertex &lhs, const vertex &rhs) throw()
{
    return lhs.x-rhs.x;
}

static inline
int compare_vtx_by_y(const vertex &lhs, const vertex &rhs) throw()
{
    return lhs.y-rhs.y;
}

typedef vector<vertex>       Vertices;


static inline
void analyze_particle(Vertices       &shape,
                      const particle &P,
                      const int       side)
{
    typedef map<unit_t,Vertices> VtxMap;
    VtxMap vmap(P.size,as_capacity);
    for(const vnode *node = P.head;node;node=node->next)
    {
        const vertex v = node->vtx;
        const unit_t y = v.y;
        Vertices *V = vmap.search(y);
        if(V)
        {
            V->push_back(v);
        }
        else
        {
            Vertices vtmp(1,as_capacity);
            vtmp.push_back(v);
            if(!vmap.insert(y,vtmp))
            {
                throw exception("unexpected vmap failure");
            }
        }
    }
    std::cerr << "vmap.size=" << vmap.size() << " / " << P.size << std::endl;
    for(VtxMap::iterator i=vmap.begin();i!=vmap.end();++i)
    {
        Vertices &V = *i;
        quicksort(V,compare_vtx_by_x);
    }
    shape.free();
    shape.ensure(vmap.size());
    for(VtxMap::iterator i=vmap.begin();i!=vmap.end();++i)
    {
        const Vertices &V = *i;
        if(side<0)
        {
            shape.push_back( V.front() );
        }
        else
        {
            shape.push_back( V.back()  );
        }
    }
    quicksort(shape,compare_vtx_by_y);
}

#if 0
static inline
size_t find_min_shape_index( const Vertices &shape, const int direction )
{
    const size_t n = shape.size();
    size_t       i = max_of<size_t>(1,n/2);
    unit_t       x = shape[i].x;
    for(--i;i>0;--i)
    {
        const unit_t xtmp=shape[i].x;
        if(direction>0)
        {
            if(xtmp<x)
            {
                return i;
            }
            goto END;
        }

        if(direction<0)
        {
            if(xtmp>x)
            {
                return i;
            }
            goto END;
        }

    END:
        x=xtmp;
    }
    return 1;
}
#endif

static inline unit_t distance(const vertex &a, const vertex &b )
{
    return square_of(b.x-a.x) + square_of(a.y-b.y);
}

static inline
void find_min_lower_distance(const Vertices &A, size_t &ia,
                             const Vertices &B, size_t &ib)
{
    const size_t na = max_of<size_t>(1,A.size()/2);
    const size_t nb = max_of<size_t>(2,B.size()/2);
    ia=1;
    ib=1;
    unit_t dmin = distance(A[ia], B[ia]);
    for(size_t a=1;a<=na;++a)
    {
        for(size_t b=1;b<=nb;++b)
        {
            const unit_t dtmp = distance(A[a],B[b]);
            if(dtmp<=dmin)
            {
                ia=a;
                ib=b;
                dmin = dtmp;
            }
        }
    }
}

typedef  point2d<double> Point;
#include "yocto/sort/index.hpp"
class Shaper
{
public:
    Point           center;
    Point           radius;
    matrix<double>  rotate;

    vector<double>   x;
    vector<double>   y;
    vector<double>   Theta;
    vector<double>   radii;
    FitConic<double> fc;
    vector<double>   param;
    vector<double>   ratio;

    explicit Shaper() : center(), radius(), rotate(2), x(), y(), Theta(), radii(), fc(), param(6), ratio()
    {

    }

    virtual ~Shaper() throw()
    {

    }

    inline double ThetaFor(const double xx, const double yy) const
    {
        const Point v(xx,yy);
        const Point dv=v-center;
        Point       V;
        tao::mul_trn(V,rotate,dv);
        return Atan2(V.x,V.y);
    }

    void computeTheta()
    {
        const size_t n = x.size();
        for(size_t i=n;i>0;--i)
        {
            const Point v(x[i],y[i]);
            const Point dv=v-center;
            Point       V;
            tao::mul_trn(V,rotate,dv);
            Theta[i] = Atan2(V.x,V.y);
            radii[i] = dv.norm();
        }

        vector<size_t> idx(n);
        vector<double> tmp(n);
        make_index(Theta, idx, __compare<double> );
        make_rank(x, idx, tmp);
        make_rank(y, idx, tmp);
        make_rank(radii, idx, tmp );
        make_rank(Theta, idx, tmp );
        ratio.make(n);
        for(size_t i=n;i>0;--i)
        {
            const double th = Theta[i];
            ratio[i] = radii[i]/computeR(th) - 1.0;
        }
    }

    inline double computeR(const double th) const throw()
    {
        return radius.x*radius.y / sqrt( square_of( radius.x * sin(th) ) + square_of( radius.y * cos(th) ) );
    }

    inline double getSin(const size_t k)
    {
        double ans = 0;
        const size_t n = x.size();
        for(size_t i=1;i<n;++i)
        {
            const double h1 = ratio[i]   * sin(k*Theta[i]);
            const double h2 = ratio[i+1] * sin(k*Theta[i+1]);
            ans += (Theta[i+1]-Theta[i]) * (h2+h1) * 0.5;
        }
        return ans/numeric<double>::pi;
    }

    inline double getCos(const size_t k)
    {
        double ans = 0;
        const size_t n = x.size();
        for(size_t i=1;i<n;++i)
        {
            const double h1 = ratio[i]   * cos(k*Theta[i]);
            const double h2 = ratio[i+1] * cos(k*Theta[i+1]);
            ans += (Theta[i+1]-Theta[i]) * (h2+h1) * 0.5;
        }
        return ans/numeric<double>::pi;
    }




private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Shaper);
};

YOCTO_PROGRAM_START()
{
    YOCTO_GFX_DECL_FORMAT(jpeg);
    YOCTO_GFX_DECL_FORMAT(png);
    YOCTO_GFX_DECL_FORMAT(tiff);

    imageIO          &IMG = image::instance();
    float             sig = 1.4f;
    if(argc>1)
    {
        const string  filename = argv[1];
        std::cerr << "-- Loading " << filename << std::endl;
        const pixmap3 source( IMG.load3(filename,0) );
        const unit_t w = source.w;
        const unit_t h = source.h;
        std::cerr << "-- Prepare Engine" << std::endl;
        xpatches xps(source,true);
        y_limit = h/2;

        //______________________________________________________________________
        //
        // hard copy
        //______________________________________________________________________
        const pixmapf img0(source,RGB::to_float,source);
        IMG.save("img0.png",img0,0);

        pixmapf       img(img0);

        //______________________________________________________________________
        //
        // Filtering if wanted
        //______________________________________________________________________
        if(argc>2)
        {
            sig = strconv::to_float(argv[2],"sig");
        }

        if(sig>0)
        {
            std::cerr << "-- Gaussian Filtering, sigma=" << sig << std::endl;
            const stencil_gauss   g5(2,sig);
            stencil::dispatcher   dsp(w,h);
            dsp(g5,img,img0,xps);
        }
        else
        {
            std::cerr << "-- Using Raw Image..." << std::endl;
        }
        IMG.save("img.png",img,0);


        //______________________________________________________________________
        //
        // detect edges
        //______________________________________________________________________
        std::cerr << "-- Detecting Edges..." << std::endl;
        EdgeDetector      ED(w,h);
        tagmap           &tags  = ED.tags;
        particles        &edges = ED.edges;
        stencil_scharr_x5 Gx;
        stencil_scharr_y5 Gy;

        ED.build_from(img,Gx,Gy,xps,1.0f);
        ED.tags.colors.shift = YGFX_RED;
        IMG.save("img-grad.png", ED, 0 );
        IMG.save("img-tags.png", tags, tags.colors, 0);
        std::cerr << "#edges=" << edges.size() << std::endl;

        //______________________________________________________________________
        //
        // copy edges over original picture
        //______________________________________________________________________
        pixmap3 tgt(source);
        for(size_t i=edges.size();i>0;--i)
        {
            const particle &pa = *edges[i];
            pa.mask(tgt, named_color::fetch(pa.tag+tags.colors.shift), 255);
        }
        IMG.save("img-edges.png", tgt, 0);

        //______________________________________________________________________
        //
        // close shapes
        //______________________________________________________________________
        std::cerr << "-- Closing..." << std::endl;
        Filter<float> F(w,h);
        F.Close(ED,xps);
        IMG.save("img-close.png", ED, 0);

        //______________________________________________________________________
        //
        // detect final blobs
        //______________________________________________________________________
        std::cerr << "-- Final Blobs..." << std::endl;
        tags.build(ED,8);
        edges.load(tags);
        IMG.save("img-blobs.png", tags, tags.colors, 0);


        //______________________________________________________________________
        //
        // Remove upper part
        //______________________________________________________________________
        std::cerr << "-- Removing Vertices" << std::endl;
        edges.reject_all_vertices_from(tags,is_bad_vertex);
        std::cerr << "#final_edges=" << edges.size() << std::endl;
        IMG.save("img-final.png", tags, tags.colors, 0);
        if( edges.size() < 2 )
        {
            throw exception("cannot find 2 particles!");
        }

        //______________________________________________________________________
        //
        // sort by Y extension
        //______________________________________________________________________
        quicksort(edges, compare_by_decreasing_particle_height);
        particle::ptr A( edges[1] );
        particle::ptr B( edges[2] );
        if( A->center().x >= B->center().x )
        {
            bswap(A,B);
        }

        //______________________________________________________________________
        //
        // A is LEFT part, B is RIGHT part
        //______________________________________________________________________
        const RGB Acolor = named_color::fetch(YGFX_RED);
        const RGB Bcolor = named_color::fetch(YGFX_GREEN);
        tgt.copy(source);
        A->mask(tgt, Acolor, 255);
        B->mask(tgt, Bcolor, 255);
        IMG.save("img-wksp.png", tgt, 0 );

        //______________________________________________________________________
        //
        // analyze shapes from particles
        //______________________________________________________________________
        Vertices As;
        analyze_particle(As,*A,1);

        Vertices Bs;
        analyze_particle(Bs,*B,-1);

        tgt.copy(source);

        for(size_t i=1;i<=As.size();++i)
        {
            tgt[As[i]] = Acolor;
        }
        for(size_t i=1;i<=Bs.size();++i)
        {
            tgt[Bs[i]] = Bcolor;
        }

        {
            ios::wcstream fp("a_shape.dat");
            for(size_t i=1;i<=As.size();++i)
            {
                fp("%ld %ld\n", As[i].x, As[i].y);
            }
        }

        {
            ios::wcstream fp("b_shape.dat");
            for(size_t i=1;i<=Bs.size();++i)
            {
                fp("%ld %ld\n", Bs[i].x, Bs[i].y);
            }
        }

        //______________________________________________________________________
        //
        // Evaluate intersection
        //______________________________________________________________________
        size_t ia=1,ib=1;
        find_min_lower_distance(As,ia,Bs,ib);
        {
            const vertex mav = As[ia];
            const vertex mbv = Bs[ib];
            draw_line(tgt, mav.x, mav.y, mbv.x, mbv.y, named_color::fetch(YGFX_WHITE_SMOKE));
        }
        IMG.save("img-shape.png", tgt, 0 );

        //______________________________________________________________________
        //
        // bounding box
        //______________________________________________________________________
        unit_t xmin = As[ia].x;
        unit_t xmax = xmin;
        unit_t ymin = As[ia].y;
        unit_t ymax = ymin;

        for(size_t i=As.size();i>=ia;--i)
        {
            const vertex v = As[i];
            xmin = min_of(v.x,xmin);
            xmax = max_of(v.x,xmax);
            ymin = min_of(v.y,ymin);
            ymax = max_of(v.y,ymax);
        }

        for(size_t i=Bs.size();i>=ib;--i)
        {
            const vertex v = Bs[i];
            xmin = min_of(v.x,xmin);
            xmax = max_of(v.x,xmax);
            ymin = min_of(v.y,ymin);
            ymax = max_of(v.y,ymax);
        }
        //const unit_t delta_x = xmax-xmin+1;
        //const unit_t delta_y = ymax-ymin+1;


        std::cerr << "Fitting Ellipse..." << std::endl;
        Shaper S;
        FitConic<double> &fc = S.fc;

        const size_t nA   = 1+As.size()-ia;
        const size_t nB   = 1+Bs.size()-ib;
        const size_t ntot = nA+nB;
        S.x.free();
        S.y.free();
        S.x.ensure(ntot);
        S.y.ensure(ntot);
        S.Theta.make(ntot);
        S.radii.make(ntot);

        fc.data.free();
        fc.data.ensure(ntot);

        for(size_t i=As.size();i>=ia;--i)
        {
            const vertex v = As[i];
            fc.append(v.x, v.y, ED[v]);
            S.x.push_back(v.x);
            S.y.push_back(v.y);
        }
        for(size_t i=Bs.size();i>=ib;--i)
        {
            const vertex v = Bs[i];
            fc.append(v.x, v.y, ED[v]);
            S.x.push_back(v.x);
            S.y.push_back(v.y);
        }


        array<double> &param = S.param;
        fc.compute(FitConicEllipse,param);
        point2d<double>     &center = S.center;
        point2d<double>     &radius = S.radius;
        matrix<double>      &rotate = S.rotate;
        fc.Reduce(center,radius,rotate,param);
        std::cerr << "center=" << center << std::endl;
        std::cerr << "radius=" << radius << std::endl;
        std::cerr << "rotate=" << rotate << std::endl;

        S.computeTheta();
        const unit_t cx = unit_t(floor(center.x+0.5));
        const unit_t cy = unit_t(floor(center.y+0.5));
        const vertex cc(cx,cy);
        const RGB    Scolor = named_color::fetch(YGFX_CYAN);

        for(size_t i=1;i<=ntot;++i)
        {
            const double theta   = S.Theta[i];
            const double r       = S.radii[i];
            const Point    V( cos(theta) * r, sin(theta) * r );
            point2d<double> v;
            tao::mul(v, rotate, V);
            v += center;
            const vertex p(unit_t(floor(v.x+0.5)),unit_t(floor(v.y+0.5)));
            draw_line(tgt, cx, cy, p.x, p.y, named_color::fetch(YGFX_WHITE), 100);
        }

        for(double th=0;th<360;++th)
        {
            const double   theta = Deg2Rad(th);
            const double   r     = S.computeR(theta);
            const Point    V( cos(theta) * r, sin(theta) * r );
            point2d<double> v;
            tao::mul(v, rotate, V);
            v += center;
            const vertex p(unit_t(floor(v.x+0.5)),unit_t(floor(v.y+0.5)));
            if(tgt.has(p))
            {
                tgt[p] = Scolor;
            }
        }

        draw_disk(tgt, cx, cy, 2, named_color::fetch(YGFX_YELLOW_OCHRE) );
        IMG.save("img-shape.png", tgt, 0 );


        {
            ios::wcstream fp("r_theta.dat");
            for(size_t i=1;i<=ntot;++i)
            {
                const double th  = S.Theta[i];
                fp("%g %g %g %g\n", th, S.radii[i], S.computeR(th), S.ratio[i]);
            }
        }

        vector<double> tA(nA,as_capacity);
        vector<double> rA(nA,as_capacity);

        for(size_t i=ia;i<=As.size();++i)
        {
            const Point v(As[i].x,As[i].y);
            const Point dv=v-center;
            Point       V;
            tao::mul_trn(V,rotate,dv);
            tA.push_back(Atan2(V.x,V.y));
            rA.push_back(dv.norm()/S.computeR(tA.back())-1);
        }
        co_qsort(tA,rA);
        {
            ios::wcstream fp("ratio_A.dat");
            for(size_t i=1;i<=nA;++i)
            {
                fp("%g %g\n", tA[i], rA[i]);
            }
        }

        vector<double> tB(nB,as_capacity);
        vector<double> rB(nB,as_capacity);

        for(size_t i=ib;i<=Bs.size();++i)
        {
            const Point v(Bs[i].x,Bs[i].y);
            const Point dv=v-center;
            Point       V;
            tao::mul_trn(V,rotate,dv);
            tB.push_back(Atan2(V.x,V.y));
            rB.push_back(dv.norm()/S.computeR(tB.back())-1);
        }
        co_qsort(tB,rB);
        {
            ios::wcstream fp("ratio_B.dat");
            for(size_t i=1;i<=nB;++i)
            {
                fp("%g %g\n", tB[i], rB[i]);
            }
        }


    }

}
YOCTO_PROGRAM_END()

