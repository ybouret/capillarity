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

YOCTO_PROGRAM_START()
{
    YOCTO_GFX_DECL_FORMAT(png);
    YOCTO_GFX_DECL_FORMAT(jpeg);
    YOCTO_GFX_DECL_FORMAT(tiff);

    imageIO &IMG = image::instance();
    float    sig = 1.4f;

    if(argc>1)
    {
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
        // Ok, so where is the pipette, left or right ?
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

        particle pp(0);
        for( const vnode *node = edges[1]->head; node; node = node->next )
        {
            pp.push_back( new vnode(*node) );
        }
        for( const vnode *node = edges[2]->head; node; node = node->next )
        {
            pp.push_back( new vnode(*node) );
        }

        const size_t np = pp.size;
        V2D G;
        for(const vnode *node = pp.head;node;node=node->next)
        {
            G.x += double(node->vtx.x);
            G.y += double(node->vtx.y);
        }
        G.x /= np;
        G.y /= np;
        std::cerr << "G=" << G << std::endl;
        double torque = 0;
        for(const vnode *node = pp.head;node;node=node->next)
        {
            V2D Q(node->vtx.x,node->vtx.y);
            V2D GQ(G,Q);
            const double mass = Q.y; // virtual mass
            torque -= mass*GQ.x;
        }
        torque /= np;
        std::cerr << "torque=" << torque << std::endl;

        particle       *pWork1 = 0;
        const patch    *pArea1 = 0;
        particle       *pWork2 = 0;
        const patch    *pArea2 = 0;
        bool            is_right = true;
        bool            is_left  = false;
        if(torque<0)
        {
            std::cerr << "\tpipette is on the right!" << std::endl;
            pWork1 = & *edges[1];
            pArea1 = & box_left;

            pWork2 = & *edges[2];
            pArea2 = & box_right;
        }
        else
        {
            std::cerr << "\tpipette is on the left!" << std::endl;
            pWork1 = & *edges[2];
            pArea1 = & box_right;

            pWork2 = & *edges[1];
            pArea2 = &  box_left;
            is_right = false;
            is_left  = true;

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


    }


}
YOCTO_PROGRAM_END()
