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
#include "yocto/sort/index.hpp"
#include "yocto/math/fit/glsf-spec.hpp"
#include "yocto/math/fcn/drvs.hpp"

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

class Shaper
{
public:
    vector<double> x;
    vector<double> y;
    vector<double> w;
    vector<double> theta;
    vector<double> rho;

    Point  center;
    double radius;

    Shaper() :
    x(),
    y(),
    w(),
    center(),
    radius(0)
    {
    }

    void ApproxCircle()
    {
        const size_t n = x.size();
        {
            FitCircle<double> fc;
            for(size_t i=x.size();i>0;--i)
            {
                fc.append(x[i], y[i], w[i]);
            }
            fc.compute(center,radius);
        }
        std::cerr << "center=" << center << std::endl;
        std::cerr << "radius=" << radius << std::endl;
        theta.make(n);
        rho.make(n);

        for(size_t i=n;i>0;--i)
        {
            const Point v(x[i],y[i]);
            const Point dv = v - center;
            theta[i] = Atan2(-dv.y,dv.x);
            rho[i]   = dv.norm();
        }

        {
            ios::wcstream fp("polar.dat");
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g\n", theta[i], rho[i]/radius-1);
            }
        }
    }

    ~Shaper() throw()
    {
    }


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Shaper);
};

#include "yocto/math/fcn/zfind.hpp"

class Fitter
{
public:
    vector<double> aorg;
    const Point    center;
    const double   radius;
    mutable double yvalue;

    explicit Fitter(const size_t nv, const Point &cc, const double rr) :
    aorg(nv),
    center(cc),
    radius(rr),
    yvalue(0)
    {
    }

    double computeR(const double angle) const
    {
        const double u  = sin(angle);
        const double rf = (1+_GLS::Polynomial<double>::Eval(u,aorg)) * radius;
        return rf;
    }

    double computeX(const double angle) const
    {
        const double rf = computeR(angle);
        return center.x + rf * sin(angle);
    }

    double computeY(const double angle) const
    {
        const double rf = computeR(angle);
        return center.y - rf * cos(angle);
    }

    Point computeP(const double angle) const
    {
        const double rf = computeR(angle);
        return Point(center.x + rf * sin(angle),center.y - rf * cos(angle));
    }

    vertex computeV(const double angle) const
    {
        const Point p = computeP(angle);
        return vertex( floor(p.x+0.5), floor(p.y+0.5) );
    }

    double compute_dY(double angle) const
    {
        return computeY(angle) - yvalue;
    }


    virtual ~Fitter() throw()
    {
    }

private:

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


        size_t degree = 6;
        if(argc>3)
        {
            degree = strconv::to_size(argv[3],"degree");
        }

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
        const vertex mav = As[ia];
        const vertex mbv = Bs[ib];
        {

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


        const size_t nA = 1+As.size()-ia;
        const size_t nB = 1+Bs.size()-ib;
        const size_t ntot = nA+nB;

        Shaper shape;
        shape.x.ensure(ntot);
        shape.y.ensure(ntot);
        shape.w.ensure(ntot);

        for(size_t i=As.size();i>=ia;--i)
        {
            const vertex v = As[i];
            xmin = min_of(v.x,xmin);
            xmax = max_of(v.x,xmax);
            ymin = min_of(v.y,ymin);
            ymax = max_of(v.y,ymax);

            shape.x.push_back(v.x);
            shape.y.push_back(v.y);
            shape.w.push_back(ED[v]);


        }

        for(size_t i=Bs.size();i>=ib;--i)
        {
            const vertex v = Bs[i];
            xmin = min_of(v.x,xmin);
            xmax = max_of(v.x,xmax);
            ymin = min_of(v.y,ymin);
            ymax = max_of(v.y,ymax);

            shape.x.push_back(v.x);
            shape.y.push_back(v.y);
            shape.w.push_back(ED[v]);
        }
        assert(ntot==shape.x.size());
        
        shape.ApproxCircle();
        const unit_t cx = floor(shape.center.x+0.5);
        const unit_t cy = floor(shape.center.y+0.5);

#if 1
        draw_circle(tgt, cx, cy, floor(shape.radius+0.5), RGB(255,255,255), 255);
        for(size_t i=shape.x.size();i>0;--i)
        {
            const double r = shape.rho[i];
            Point        v( r*sin(shape.theta[i]), -r*cos(shape.theta[i]) );
            v += shape.center;
            draw_line(tgt, cx, cy, v.x, v.y,  RGB(255,255,255), 127);
        }
#endif

        IMG.save("img-shape.png", tgt, 0 );

        vector<double> XX(ntot,as_capacity);
        vector<double> YY(ntot,as_capacity);
        vector<double> Theta(ntot,as_capacity);
        vector<double> Param(ntot,as_capacity);
        vector<double> Ratio(ntot,as_capacity);
        for(size_t i=shape.x.size();i>0;--i)
        {
            if( Fabs(shape.theta[i])<numeric<double>::pi/2 )
            {
                XX.push_back(shape.x[i]);
                YY.push_back(shape.y[i]);
                Theta.push_back(shape.theta[i]);
                Param.push_back(sin(shape.theta[i]));
                Ratio.push_back(shape.rho[i]/shape.radius-1);
            }
        }

        const size_t nfit = XX.size();
        {
            vector<size_t> idx(nfit);
            vector<double> tmp(nfit);
            make_index(Param,idx,__compare<double>);
            make_rank(Param, idx, tmp);
            make_rank(Theta, idx, tmp);
            make_rank(XX, idx, tmp);
            make_rank(YY, idx, tmp);
            make_rank(Ratio, idx, tmp);
        }
        assert(nfit==Param.size());
        vector<double> Rfit(nfit);

        {
            ios::wcstream fp("tofit.dat");
            for(size_t i=1;i<=nfit;++i)
            {
                fp("%g %g\n", Param[i],Ratio[i]);
            }
        }

        Fitter                ff(1+degree,shape.center,shape.radius);
        array<double>        &aorg = ff.aorg;
        GLS<double>::Samples  samples;
        GLS<double>::Sample  &sample = samples.append(Param, Ratio, Rfit);
        if( !_GLS::Polynomial<double>::Start(sample, aorg) )
        {
            throw exception("unexpected polyfit start failure!!!");
        }
        std::cerr << "aorg=" << aorg << std::endl;

        {
            const size_t np  = 1000;
            const double ulo = Param[1];
            const double uhi = Param[nfit];
            const double du  = uhi-ulo;
            ios::wcstream fp("rfit.dat");
            for(size_t i=0;i<=np;++i)
            {
                const double u  = ulo + (i*du)/np;
                const double rf = _GLS::Polynomial<double>::Eval(u,aorg);
                fp("%g %g\n", u, rf);
            }
        }

        tgt.copy(source);
        {
            draw_line(tgt, mav.x, mav.y, mbv.x, mbv.y, named_color::fetch(YGFX_MAGENTA),200);
        }
        for(size_t i=1;i<=nfit;++i)
        {
            const vertex v(XX[i],YY[i]);
            if(tgt.has(v))
            {
                tgt[v] = named_color::fetch(YGFX_RED);
            }
        }

        {
            const double dth = numeric<double>::pi/2000;
            for(double th=Theta[1];th<=Theta[nfit];th+=dth)
            {
                //const double rf    = ff.computeR(th);
                //vertex v( shape.center.x + rf*sin(th), shape.center.y-rf*cos(th) );
                const vertex v = ff.computeV(th);
                if(tgt.has(v))
                {
                    tgt[v] = pixel<RGB>::blend(tgt[v], named_color::fetch(YGFX_CYAN), 0xff);
                }

            }
        }

        IMG.save("img-fitted.png",tgt,0);

        tgt.copy(source);
        zfind<double> solve( Deg2Rad(0.01) );
        ff.yvalue = mav.y;
        numeric<double>::function zfn( &ff, & Fitter::compute_dY);
        const double angleA = solve(zfn,Theta[1],0);
        std::cerr << "angleA=" << Rad2Deg(angleA) << std::endl;

        ff.yvalue = mbv.y;
        const double angleB = solve(zfn,0,Theta[nfit]);
        std::cerr << "angleB=" << Rad2Deg(angleB) << std::endl;

        const vertex touchA = ff.computeV(angleA);
        const vertex touchB = ff.computeV(angleB);

        draw_disk(tgt,touchA.x, touchA.y, 2, named_color::fetch(YGFX_YELLOW), 100);
        draw_disk(tgt,touchB.x, touchB.y, 2, named_color::fetch(YGFX_YELLOW), 100);


        derivative<double> drvs;
        numeric<double>::function dfx( &ff, & Fitter::computeX);
        numeric<double>::function dfy( &ff, & Fitter::computeY);

        const double dfxA = drvs(dfx,angleA,1e-4);
        const double dfyA = drvs(dfy,angleA,1e-4);
        const double thA  = Atan2(dfxA,dfyA);
        std::cerr << "thA=" << Rad2Deg(thA) << std::endl;
        const Point  arrA = (ff.radius/2)*Point(cos(thA+numeric<double>::pi),sin(thA+numeric<double>::pi));
        draw_line(tgt,touchA.x, touchA.y, touchA.x + unit_t(arrA.x), touchA.y+unit_t(arrA.y), named_color::fetch(YGFX_MAGENTA), 200 );


        const double dfxB = drvs(dfx,angleB,1e-4);
        const double dfyB = drvs(dfy,angleB,1e-4);
        const double thB  = Atan2(dfxB,dfyB);
        std::cerr << "thB=" << Rad2Deg(thB) << std::endl;
        const Point  arrB = (ff.radius/2)*Point(cos(thB),sin(thB));
        draw_line(tgt,touchB.x, touchB.y, touchB.x + unit_t(arrB.x), touchB.y+unit_t(arrB.y), named_color::fetch(YGFX_MAGENTA), 200 );

        const double wetB = Rad2Deg(numeric<double>::pi-thB);
        const double wetA = Rad2Deg(numeric<double>::pi+thA);
        std::cerr << "WettingA=" << wetA << std::endl;
        std::cerr << "WettingB=" << wetB << std::endl;

        IMG.save("img-result.png",tgt,0);



    }
    
}
YOCTO_PROGRAM_END()

