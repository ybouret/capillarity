#include "yocto/program.hpp"
#include "yocto/gfx/image/png.hpp"
#include "yocto/gfx/image/jpeg.hpp"
#include "yocto/gfx/image/tiff.hpp"
#include "yocto/gfx/ops/stencil.hpp"
#include "yocto/gfx/ops/edges.hpp"
#include "yocto/gfx/ops/filter.hpp"
#include "yocto/gfx/draw/line.hpp"
#include "yocto/gfx/draw/circle.hpp"
#include "yocto/container/matrix.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/alg/shapes2d.hpp"
#include "yocto/math/fit/glsf-spec.hpp"
#include "yocto/math/fcn/zfind.hpp"

using namespace yocto;
using namespace gfx;
using namespace math;

typedef point2d<double> V2D;

static inline
void trim_particle( particle &p, const double y ) throw()
{
    vlist stk;
    while(p.size)
    {
        vnode *node = p.pop_front();
        if(node->vtx.y<=y)
        {
            stk.push_back(node);
        }
        else
        {
            delete node;
        }
    }
    p.swap_with(stk);
}


class Geometry
{
public:
    V2D            center;
    double         radius;
    vector<double> aorg;

    Geometry() : center(), radius(0), aorg(3)
    {
    }

    ~Geometry() throw()
    {
    }

    V2D Coord( const double alpha ) const
    {
        const double rr = radius+_GLS::Polynomial<double>::Eval(alpha,aorg);
        return V2D(center.x + rr * sin(alpha),
                   center.y - rr * cos(alpha) );
    }


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Geometry);
};

// little class to compute intersection
class FindInter
{
public:
    vector<double> aorg;
    double         radius;
    double         deltaY;

    FindInter() : aorg(3), radius(0), deltaY(0)
    {
    }

    ~FindInter() throw()
    {
    }

    double Compute( const double alpha )
    {
        return (radius+_GLS::Polynomial<double>::Eval(alpha,aorg)) * cos(alpha) - deltaY;
    }


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(FindInter);
};


YOCTO_PROGRAM_START()
{
    YOCTO_GFX_DECL_FORMAT(png);
    YOCTO_GFX_DECL_FORMAT(jpeg);
    YOCTO_GFX_DECL_FORMAT(tiff);

    imageIO &IMG = image::instance();
    float    sig = 1.4f;

    if(argc>1)
    {
        if(!(argc>2))
        {
            throw exception("usage: %s filename [right|left]", program);
        }
        bool            is_right = true;
        bool            is_left  = false;

        const string side = argv[2];
        if( side == "right" )
        {
            // do nothing
        }
        else
        {
            if( side == "left" )
            {
                is_left  =  true;
                is_right = false;
            }
            else
            {
                throw exception("unexpected side='%s'", side.c_str() );
            }
        }


        //______________________________________________________________________
        //
        // loading original file, in colors, and keep pipette on the left
        // to work on the right side, easier...
        //______________________________________________________________________
        const string  filename = argv[1];
        const pixmap3 origin( IMG.load3(filename,NULL) );
        if(is_left)
        {
            std::cerr << "-- Flipping Picture to get data on the right side" << std::endl;
            ((pixmap3 &)origin).flip_horizontal();
        }
        IMG.save("img-origin.png",origin,NULL);


        //______________________________________________________________________
        //
        // making a grayscale picture
        //______________________________________________________________________
        pixmapf      img0(origin,RGB::to_float,origin);

        //______________________________________________________________________
        //
        // preparing
        //______________________________________________________________________


        const unit_t h = origin.h;
        const unit_t w = origin.w;
        xpatches     xps(origin,false);
        std::cerr << "parallel=" << xps.size() << std::endl;
        pixmapf img(w,h);

        //______________________________________________________________________
        //
        // image, raw or filtered if sigma>0
        //______________________________________________________________________
        if(sig>0)
        {
            std::cerr << "-- Gauss " << sig << std::endl;
            stencil_gauss g5(2,sig);
            g5.apply(img,img0,xps);
        }
        else
        {
            img.copy(img0);
        }
        IMG.save("img-filter.png",img,NULL);

        //______________________________________________________________________
        //
        // edge detections
        //______________________________________________________________________
        std::cerr << "-- Edges Detection" << std::endl;
        EdgeDetector      ED(w,h);
        tagmap           &tags  = ED.tags;
        particles        &edges = ED.edges;
        stencil_scharr_x5 gx;
        stencil_scharr_y5 gy;
        ED.build_from(img,gx,gy,xps);
        IMG.save("img-grad.png", ED, 0 );
        IMG.save("img-tags.png",tags,tags.colors,NULL);
        std::cerr << "#edges=" << edges.size() << std::endl;

        //______________________________________________________________________
        //
        // copy edges over original picture
        //______________________________________________________________________
        pixmap3 tgt(origin);
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
        std::cerr << "#edges=" << edges.size() << std::endl;

        edges.sort_by_extension();

        while( edges.size() > 2 ) edges.pop_back();
        if(edges.size()<=1)
            throw exception("not enough edges, image is not fit!!!");

        tgt.copy(origin);
        for(size_t i=edges.size();i>0;--i)
        {
            const particle &pa = *edges[i];
            pa.mask(tgt, named_color::fetch(pa.tag+tags.colors.shift), 255);
            //std::cerr << "#points=" << pa.size << std::endl;
        }
        IMG.save("img-final.png", tgt, 0);


        //______________________________________________________________________
        //
        //
        // Ok, so where is the pipette, left or right ?
        //
        //______________________________________________________________________
        std::cerr << "-- Finding Right Box..." << std::endl;
        patch box_left  = edges[1]->aabb();
        patch box_right = edges[2]->aabb();
        if(box_left.lower.x>box_right.lower.x)
        {
            bswap(box_left,box_right);
            bswap(edges[1],edges[2]);
        }

        draw_patch(tgt,box_left,  named_color::fetch( YGFX_GREEN ), 127 );
        draw_patch(tgt,box_right, named_color::fetch( YGFX_RED   ), 127 );
        IMG.save("img-final.png", tgt, 0);


        particle    &pR = *edges[2];
        const patch &aR = box_right;

        particle    &pL = *edges[1];
        const patch &aL = box_left;


        //______________________________________________________________________
        //
        //
        // Now we are working...
        //
        //______________________________________________________________________
        std::cerr << "-- Building Workspace" << std::endl;
        // keep working space
        trim_particle(pL,(aL.lower.y+aL.upper.y)/2);
        trim_particle(pR,(aR.lower.y+aR.upper.y)/2);

        if(pL.size<=1||pR.size<=1)
        {
            throw exception("particles are too small");
        }

        pR.mask(tgt,  named_color::fetch(YGFX_YELLOW), 255);
        pL.mask(tgt, named_color::fetch(YGFX_MAGENTA), 255);

        IMG.save("img-final.png", tgt, 0);


        //______________________________________________________________________
        //
        // Guess the lens limit
        //______________________________________________________________________
        const vnode *p1    = pR.head;
        const vnode *p2    = pL.head;
        unit_t       min_d = vertex(p1->vtx,p2->vtx).norm2();

        for(const vnode *n1 = pR.head; n1; n1=n1->next)
        {
            for(const vnode *n2 = pL.head; n2; n2=n2->next)
            {
                const unit_t tmp_d = vertex(n1->vtx,n2->vtx).norm2();
                if(tmp_d<min_d)
                {
                    p1    = n1;
                    p2    = n2;
                    min_d = tmp_d;
                }
            }
        }

        draw_line(tgt,p1->vtx,p2->vtx, named_color::fetch(YGFX_CYAN), 0xff);
        IMG.save("img-final.png", tgt, 0);

        //______________________________________________________________________
        //
        // study geometry
        //______________________________________________________________________
        Geometry Geom;
        V2D           &center = Geom.center;
        double        &radius = Geom.radius;
        array<double> &aorg   = Geom.aorg;

        //______________________________________________________________________
        //
        // ok, now we need the points and do something...
        //______________________________________________________________________
        FitCircle<double> fc;

        const unit_t y_low = p1->vtx.y;
        vector<double> X(pR.size,as_capacity);
        vector<double> Y(pR.size,as_capacity);
        vector<double> A(pR.size,as_capacity);
        vector<double> R(pR.size,as_capacity);

        for(const vnode *n1 = pR.head; n1; n1=n1->next)
        {
            const vertex v = n1->vtx;
            if(v.y>=y_low)
            {
                X.push_back(v.x);
                Y.push_back(v.y);
                fc.append(v.x,v.y);
            }
        }

        // process all the points
        const size_t N = X.size();

        fc.compute(center,radius);

        std::cerr << "center=" << center << std::endl;
        std::cerr << "radius=" << radius << std::endl;

        for(size_t i=1;i<=N;++i)
        {
            const double dx = X[i] - center.x;
            const double dy = center.y - Y[i];
            const double aa = Atan2(dy,dx);
            A.push_back(aa);
            R.push_back( Hypotenuse(dx,dy) );
        }



        {
            ios::wcstream fp("shape.dat");
            for(size_t i=1;i<=N;++i)
            {
                fp("%g %g %g %g %g\n", X[i], Y[i], Rad2Deg(A[i]), R[i], R[i] - radius);
            }
        }


        co_qsort(A,R);

        // TODO: define delta in pixels and sweep angle
        const double Delta = 2;
        const double Sweep = 30;
        const double Aini  = A[1];
        const double Aend  = min_of<double>( Aini + Deg2Rad(Sweep), 90 );

        vector<double> AA(pR.size,as_capacity);
        vector<double> RR(pR.size,as_capacity);
        vector<double> XX(pR.size,as_capacity);
        vector<double> YY(pR.size,as_capacity);

        for(size_t i=1;i<=N;++i)
        {
            const double aa = A[i];
            const double rr = R[i];
            const double xx = rr*sin(aa)+center.x;
            const double yy = center.y-rr*cos(aa);

            if(yy>=y_low+Delta&&aa<=Aend)
            {
                AA.push_back(aa);
                RR.push_back(rr);
                XX.push_back(xx);
                YY.push_back(yy);
            }

        }
        const size_t NN = AA.size();
        vector<double> RF(NN);
        vector<double> R0(NN);
        for(size_t i=1;i<=NN;++i) R0[i] = RR[i] - radius;

        GLS<double>::Samples samples;
        GLS<double>::Sample &sample  = samples.append(AA,R0,RF);
        if(!_GLS::Polynomial<double>::Start(sample,aorg))
        {
            throw exception("Unable to fit polynomial");
        }

        {
            ios::wcstream fp("shapefit.dat");
            for(size_t i=1;i<=NN;++i)
            {
                fp("%g %g %g %g %g %g\n", XX[i], YY[i], Rad2Deg(AA[i]), RR[i], RR[i] - radius, RF[i]);
            }
        }


        std::cerr << "y_low=" << y_low << std::endl;


        //______________________________________________________________________
        //
        // Now we have an approximation of the shape...
        //______________________________________________________________________



#if 0
        particle       *pWork1 = 0;
        const patch    *pArea1 = 0;
        particle       *pWork2 = 0;
        const patch    *pArea2 = 0;

        if(is_left)
        {
            std::cerr << "\tuse left side, pipette is on the right!" << std::endl;
            pWork1 = & *edges[1];
            pArea1 = & box_left;

            pWork2 = & *edges[2];
            pArea2 = & box_right;
        }
        else
        {
            std::cerr << "\tuse right side, pipette is on the left!" << std::endl;
            pWork1 = & *edges[2];
            pArea1 = & box_right;

            pWork2 = & *edges[1];
            pArea2 = &  box_left;
        }

        pWork1->mask(tgt, named_color::fetch(YGFX_ORANGE), 255);
        IMG.save("img-final.png", tgt, 0);

        //______________________________________________________________________
        //
        //
        // Now we are working...
        //
        //______________________________________________________________________
        std::cerr << "-- Building Workspace" << std::endl;
        // keep working space
        trim_particle(*pWork1,(pArea1->lower.y+pArea1->upper.y)/2);
        trim_particle(*pWork2,(pArea2->lower.y+pArea2->upper.y)/2);

        if(pWork1->size<=1||pWork2->size<=1)
        {
            throw exception("particles are too small");
        }

        pWork1->mask(tgt,  named_color::fetch(YGFX_YELLOW), 255);
        pWork2->mask(tgt, named_color::fetch(YGFX_MAGENTA), 255);

        IMG.save("img-final.png", tgt, 0);

        std::cerr << "-- Guessing Intersection" << std::endl;

        //______________________________________________________________________
        //
        // Guess the lens limit
        //______________________________________________________________________
        const vnode *p1    = pWork1->head;
        const vnode *p2    = pWork2->head;
        unit_t       min_d = vertex(p1->vtx,p2->vtx).norm2();

        for(const vnode *n1 = pWork1->head; n1; n1=n1->next)
        {
            for(const vnode *n2 = pWork2->head; n2; n2=n2->next)
            {
                const unit_t tmp_d = vertex(n1->vtx,n2->vtx).norm2();
                if(tmp_d<min_d)
                {
                    p1    = n1;
                    p2    = n2;
                    min_d = tmp_d;
                }
            }
        }

        draw_line(tgt,p1->vtx,p2->vtx, named_color::fetch(YGFX_CYAN), 0xff);
        IMG.save("img-final.png", tgt, 0);

        //______________________________________________________________________
        //
        // ok, now we need the points and do something...
        //______________________________________________________________________
        FitCircle<double> fc;

        const unit_t y_low = p1->vtx.y;
        vector<double> X(pWork1->size,as_capacity);
        vector<double> Y(pWork1->size,as_capacity);
        vector<double> A(pWork1->size,as_capacity);
        vector<double> R(pWork1->size,as_capacity);

        for(const vnode *n1 = pWork1->head; n1; n1=n1->next)
        {
            const vertex v = n1->vtx;
            if(v.y>=y_low)
            {
                X.push_back(v.x);
                Y.push_back(v.y);
                fc.append(v.x,v.y);
            }
        }

        // process all the points
        const size_t N = X.size();

        V2D    center;
        double radius = 0;
        fc.compute(center,radius);

        std::cerr << "center=" << center << std::endl;
        std::cerr << "radius=" << radius << std::endl;

        for(size_t i=1;i<=N;++i)
        {
            const double dx = X[i] - center.x;
            const double dy = center.y - Y[i];
            const double aa = Atan2(dy,dx);
            A.push_back(aa);
            R.push_back( Hypotenuse(dx,dy) );
        }




        {
            ios::wcstream fp("shape.dat");
            for(size_t i=1;i<=N;++i)
            {
                fp("%g %g %g %g %g\n", X[i], Y[i], Rad2Deg(A[i]), R[i], R[i] - radius);
            }
        }

        //draw_circle(tgt,unit_t(center.x), unit_t(center.y), unit_t(radius), named_color::fetch(YGFX_IVORY), 0xff);
        //IMG.save("img-final.png", tgt, 0);


        // take the fitting section
        if(is_left)
        {
            for(size_t i=1;i<=N;++i)
            {
                A[i] = -A[i];
            }
        }

        co_qsort(A,R);

        // TODO: define delta in pixels and sweep angle
        const double Delta = 2;
        const double Sweep = 30;
        const double Aini  = A[1];
        const double Aend  = min_of<double>( Aini + Deg2Rad(Sweep), 90 );

        vector<double> AA(pWork1->size,as_capacity);
        vector<double> RR(pWork1->size,as_capacity);
        vector<double> XX(pWork1->size,as_capacity);
        vector<double> YY(pWork1->size,as_capacity);

        for(size_t i=1;i<=N;++i)
        {
            const double aa = A[i];
            const double rr = R[i];
            const double xx = rr*sin(aa)+center.x;
            const double yy = center.y-rr*cos(aa);

            if(yy>=y_low+Delta&&aa<=Aend)
            {
                AA.push_back(aa);
                RR.push_back(rr);
                XX.push_back(xx);
                YY.push_back(yy);
            }

        }
        const size_t NN = AA.size();
        vector<double> RF(NN);
        vector<double> R0(NN);
        for(size_t i=1;i<=NN;++i) R0[i] = RR[i] - radius;

        GLS<double>::Samples samples;
        GLS<double>::Sample &sample  = samples.append(AA,R0,RF);

        FindInter      inter;
        array<double> &aorg = inter.aorg;
        inter.radius        = radius;
        inter.deltaY        = center.y - y_low;

        if(!_GLS::Polynomial<double>::Start(sample,aorg))
        {
            throw exception("Unable to fit polynomial");
        }

        {
            ios::wcstream fp("shapefit.dat");
            for(size_t i=1;i<=NN;++i)
            {
                fp("%g %g %g %g %g %g\n", XX[i], YY[i], Rad2Deg(AA[i]), RR[i], RR[i] - radius, RF[i]);
            }
        }


        std::cerr << "y_low=" << y_low << std::endl;

        // find the interception angle
        zfind<double>             solver(1e-5);
        numeric<double>::function zfn( &inter, & FindInter::Compute );
        triplet<double>           zAlpha = { 0, 0, Aend };
        triplet<double>           zValue = { zfn(zAlpha.a), 0, zfn(zAlpha.c) };
        if(zValue.a*zValue.c>0)
        {
            throw exception("Cannot find interception, corrupted picture?");
        }

        std::cerr << "alpha_ini=" << Rad2Deg(Aini) << std::endl;
        std::cerr << "alpha_end=" << Rad2Deg(Aend) << std::endl;



        const double alpha0 = solver.run(zfn, zAlpha, zValue);
        std::cerr << "alpha0=" << Rad2Deg(alpha0) << std::endl;





        vector<double> drvs( aorg.size() );
        _GLS::Polynomial<double>::ComputeDrvs(drvs,aorg);
        std::cerr << "polynomial=" << aorg << std::endl;
        std::cerr << "derivative=" << drvs << std::endl;

        const double rr0   = radius+_GLS::Polynomial<double>::Eval(alpha0,aorg);
        const double rp0   = _GLS::Polynomial<double>::Eval(alpha0,drvs);
        const double ca0   = cos(alpha0);
        const double sa0   = sin(alpha0);
        const double num   = rr0 * sa0 - rp0 * ca0;
        const double den   = rr0 * ca0 + rp0 * sa0;
        const double theta = Atan(num/den);
        std::cerr << "theta=" << (180.0-Rad2Deg(theta)) << std::endl;

        const double xx0 = center.x + sa0 * rr0;
        tgt.copy(origin);
        {
            const double alpha_step = 0.001;
            for(double alpha=Aini;alpha<=Aend;alpha+=alpha_step)
            {
                const double aa = is_right ? alpha : -alpha;
                const double rr = radius+_GLS::Polynomial<double>::Eval(alpha,aorg);
                const double ca = cos(aa);
                const double sa = sin(aa);
                const double xx = center.x + rr * sa;
                const double yy = center.y - rr * ca;
                const point2d<unit_t> pp(xx,yy);
                if(tgt.has(pp))
                {
                    //std::cerr << "plot " << pp << std::endl;
                    //pixel<RGB>::blend(tgt[pp],named_color::fetch(YGFX_YELLOW),0xff);
                    tgt[pp] = named_color::fetch(YGFX_YELLOW);
                }
            }
        }

        const double beta = is_right ? theta : numeric<double>::pi-theta;
        const double xx1  = is_right ? xx0   : center.x - sa0*rr0;
        draw_disk(tgt,unit_t(xx1), unit_t(y_low), 2, named_color::fetch(YGFX_FIREBRICK), 127);
        draw_line(tgt,
                  unit_t(xx1), unit_t(y_low),
                  unit_t(xx1+radius*cos(beta)),unit_t(y_low+radius*sin(beta)),
                  named_color::fetch(YGFX_FIREBRICK), 127
                  );
        IMG.save("img-angle.png", tgt, 0);

        {
            ios::wcstream fp("resfit.dat");
            for(double alpha=alpha0;alpha<=Aend;alpha+=0.01)
            {
                const double aa = alpha;
                const double rr = _GLS::Polynomial<double>::Eval(aa,aorg)+radius;
                const double xx = rr*sin(aa)+center.x;
                const double yy = center.y-rr*cos(aa);
                fp("%g %g\n", xx,yy);
            }
        }
#endif


    }


}
YOCTO_PROGRAM_END()
