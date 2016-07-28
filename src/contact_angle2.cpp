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
#include "yocto/math/opt/cgrad.hpp"

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

#include "yocto/code/ipower.hpp"
#include "yocto/sort/index.hpp"

class Shaper
{
public:
    point2d<double> center;
    point2d<double> radius;
    matrix<double>  rotate;
    vector<double>  X;
    vector<double>  Y;
    vector<double>  alpha; //!< angle/center
    vector<double>  rho;   //!< radis/center
    vector<double>  W;

    inline void load_particle(const vlist &L, const pixmap<float> &wksp)
    {
        for(const vnode *node = L.head;node;node=node->next)
        {
            const double x = double(node->vtx.x);
            const double y = double(node->vtx.y);
            const double w = wksp[node->vtx];

            const point2d<double> delta_org(x-center.x,y-center.y);
            point2d<double>       delta;
            tao::mul_trn(delta,rotate,delta_org);
            const double dx =  delta.x;
            const double dy = -delta.y;
            const double rr = Hypotenuse(dx,dy);
            const double at = atan(dx/(dy+rr));
            const double theta = at+at;
            if(Fabs(theta)>numeric<double>::pi/2)
                continue;
            
            X.push_back(x);
            Y.push_back(y);
            W.push_back(w);
            alpha.push_back( theta );
            rho.push_back(rr);
        }
    }


    void sort()
    {
        const size_t   n = X.size();
        vector<double> tmp(n);
        vector<size_t> idx(n);
        make_index(alpha,idx,__compare<float>);
        make_rank(X,idx,tmp);
        make_rank(Y,idx,tmp);
        make_rank(W,idx,tmp);
        make_rank(alpha,idx,tmp);
        make_rank(rho,idx,tmp);
    }

    inline Shaper(size_t n) :
    center(),
    radius(),
    rotate(2),
    X(n,as_capacity),
    Y(n,as_capacity),
    alpha(n,as_capacity),
    rho(n,as_capacity),
    W(n,as_capacity)
    {
    }

    inline ~Shaper() throw()
    {
    }

    double ellradius(const double theta)
    {
        const double a    = radius.x;
        const double b    = radius.y;
        const double ab   = a*b;
        return ab/Hypotenuse(a*cos(theta),b*sin(theta));
    }

    double F(const double theta,const array<double> &params)
    {
        const double arg   = theta*(2/numeric<double>::pi);
        double ans = params[1];
        for(size_t i=2;i<=params.size();++i)
        {
            ans += params[i] * ipower(arg,i-1);
        }
        return ans;
    }

    double rhoFit(const double theta, const array<double> &params)
    {
        return ellradius(theta)*F(theta,params);
    }

    double H(const array<double> &params)
    {
        double ans = 0;

        const size_t n = X.size();
        for(size_t i=n;i>0;--i)
        {
            const double theta = alpha[i];
            const double rf    = rhoFit(theta,params);
            ans += Fabs(W[i]*(rf-rho[i]));
        }

        return ans/n;
    }

    bool HCB(const array<double> &params)
    {
        std::cerr << "Try " << params << std::endl;
        return true;
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

        // workspace is inverted: black is more weight!
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

        Shaper shape(p1.size+p2.size);

        fit_conic.Reduce(shape.center, shape.radius, shape.rotate, A);
        std::cerr << "center=" << shape.center << std::endl;
        std::cerr << "radius=" << shape.radius << std::endl;
        std::cerr << "rotate=" << shape.rotate << std::endl;



        // graphical init
        pixmap3 wfit(w,h);
        p1.transfer(wfit,Float2RGB,original);
        p2.transfer(wfit,Float2RGB,original);

        const double vc  = FitConic<double>::Eval(A,shape.center.x,shape.center.y);
        RGB          col = named_color::fetch(YGFX_RED);
        for(unit_t j=0;j<h;++j)
        {
            for(unit_t i=0;i<w;++i)
            {
                const double val = FitConic<double>::Eval(A,i,j);
                if(val*vc>=0)
                {
                    wfit[j][i] = pixel<RGB>::blend(wfit[j][i], col, 50);
                }
            }
        }

        const double a = shape.radius.x;
        const double b = shape.radius.y;
        const double ab = a*b;
        col = named_color::fetch(YGFX_GREEN);
        for(double theta=-1.57;theta<=1.57;theta+=0.001)
        {
            const double rho = ab/Hypotenuse(a*cos(theta),b*sin(theta));
            const point2d<double> P(rho*sin(theta),-rho*cos(theta));
            point2d<double> p;
            tao::mul(p,shape.rotate, P);
            p += shape.center;
            vertex v( gist::float2unit(p.x), gist::float2unit(p.y) );
            if(wfit.has(v))
            {
                wfit[v] = col;
            }
        }

        IMG.save("wfit.png",wfit,NULL);
        
        
        // preparing data
        shape.load_particle(p1,wksp);
        shape.load_particle(p2,wksp);

        shape.sort();

        {
            ios::wcstream fp("polar.dat");
            for(size_t i=1;i<=shape.X.size();++i)
            {
                const double theta = shape.alpha[i];
                //const double rho   = shape.rho[i];
                //fp("%g %g %g\n", alpha, rho, shape.ellradius(shape.alpha[i]) );
                const double rho = shape.ellradius(theta);
                const point2d<double> P(rho*sin(theta),-rho*cos(theta));
                point2d<double> p;
                tao::mul(p,shape.rotate, P);
                p += shape.center;
                fp("%g %g %g %g %g %g %g\n", shape.X[i], shape.Y[i], p.x, p.y, theta, shape.rho[i], rho);
            }
        }

        const size_t   nvar = 3;
        vector<double> aorg(nvar,0);
        vector<bool>   used(nvar,true);
        vector<double> scal(nvar,1e-4);
        aorg[1] = 1;
        //aorg[2] = 1;

        cgrad<double>               CG;
        cgrad<double>::scalar_field H(   &shape, & Shaper::H   );
        cgrad<double>::callback     HCB( &shape, & Shaper::HCB );
        std::cerr << "Optimizing Shape" << std::endl;
        if(CG.run(H, aorg, used, scal, 1e-7, &HCB))
        {
            std::cerr << "SUCCESS: " << aorg << std::endl;
        }
        else
        {
            std::cerr << "FAILURE" << std::endl;
            return 1;
        }

        {
            ios::wcstream fp("polar_fit.dat");
            for(size_t i=1;i<=shape.X.size();++i)
            {
                const double theta = shape.alpha[i];
                const double rho   = shape.rhoFit(theta,aorg);
                const point2d<double> P(rho*sin(theta),-rho*cos(theta));
                point2d<double> p;
                tao::mul(p,shape.rotate, P);
                p += shape.center;
                fp("%g %g %g %g %g %g %g\n", shape.X[i], shape.Y[i], p.x, p.y, theta, shape.rho[i], rho);
            }
        }

    }
    
    
}
YOCTO_PROGRAM_END()
