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


class FindInter
{
public:
    vector<double> aorg;
    double         radius;
    double         deltaY;

    FindInter() : aorg(6), radius(0), deltaY(0)
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
        // loading original file, in colors
        //______________________________________________________________________
        const string  filename = argv[1];
        const pixmap3 origin( IMG.load3(filename,NULL) );
        IMG.save("img-origin.png",origin,NULL);

        //______________________________________________________________________
        //
        // making a grayscale picture
        //______________________________________________________________________
        pixmapf      img0(origin,RGB::to_float,origin);
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

        // TODO: sort edges by extension !!!!
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
        std::cerr << "-- Finding Pipette" << std:: endl;
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

#if 0
        //______________________________________________________________________
        //
        // we compute the circularity index
        //______________________________________________________________________
        FitCircle<double> fcL, fcR;

        for( const vnode *node = edges[1]->head; node; node = node->next )
        {
            fcL.append(node->vtx.x,node->vtx.y);
        }

        for( const vnode *node = edges[2]->head; node; node = node->next )
        {
            fcR.append(node->vtx.x,node->vtx.y);
        }

        V2D centerL, centerR;
        double radiusL=0, radiusR=0;
        const double rmsL = fcL.compute(centerL,radiusL);
        const double rmsR = fcR.compute(centerR,radiusR);

        std::cerr << "rmsL=" << rmsL << std::endl;
        std::cerr << "rmsR=" << rmsR << std::endl;
#endif

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

        draw_circle(tgt,unit_t(center.x), unit_t(center.y), unit_t(radius), named_color::fetch(YGFX_IVORY), 0xff);
        IMG.save("img-final.png", tgt, 0);


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

        {
            ios::wcstream fp("resfit.dat");
            for(double alpha=0;alpha<=Aend;alpha+=0.01)
            {
                const double aa = alpha;
                const double rr = _GLS::Polynomial<double>::Eval(aa,aorg)+radius;
                const double xx = rr*sin(aa)+center.x;
                const double yy = center.y-rr*cos(aa);
                fp("%g %g\n", xx,yy);
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

        const double alpha0 = solver.run(zfn, zAlpha, zValue);
        std::cerr << "alpha0=" << Rad2Deg(alpha0) << std::endl;
        


#if 0
        // now use elliptic approximation
        FitConic<double> ell;
        for(size_t i=1;i<=NN;++i)
        {
            ell.append(XX[i],YY[i]);
        }
        vector<double> params(6);
        ell.compute(FitConicEllipse,params);
        std::cerr << "params=" << params << std::endl;

        V2D            Center;
        V2D            Radius;
        matrix<double> Rotate(2);

        ell.Reduce(Center,Radius,Rotate,params);
        std::cerr << "Center=" << Center << std::endl;
        std::cerr << "Radius=" << Radius << std::endl;
        std::cerr << "Rotate=" << Rotate << std::endl;

        {
            ios::wcstream fp("shape_ell.dat");
            for(double theta=0;theta<=6.3;theta+=0.01)
            {
                const double  c  = cos(theta);
                const double  s  = sin(theta);
                const V2D     rr(Radius.x*c,Radius.y*s);
                V2D           v;
                tao::mul(v, Rotate, rr);
                v += Center;
                fp("%g %g\n", v.x, v.y);
            }
        }
        const double y1 = p1->vtx.y;
        std::cerr << "y1=" << y1 << std::endl;
        find_x_for(y1,params);
#endif


    }


}
YOCTO_PROGRAM_END()
