#include "yocto/program.hpp"
#include "yocto/gfx/image/png.hpp"
#include "yocto/gfx/image/jpeg.hpp"
#include "yocto/gfx/image/tiff.hpp"
#include "yocto/gfx/ops/stencil.hpp"
#include "yocto/gfx/ops/edges.hpp"
#include "yocto/gfx/ops/filter.hpp"
#include "yocto/gfx/draw/line.hpp"
#include "yocto/container/matrix.hpp"

using namespace yocto;
using namespace gfx;

typedef point2d<double> V2D;

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

        const particle *pWork = 0;
        const patch    *pArea = 0;
        if(torque<0)
        {
            std::cerr << "pipette is on the right!" << std::endl;
            pWork = & *edges[1];
            pArea = & box_left;
        }
        else
        {
            std::cerr << "pipette is on the left!" << std::endl;
            pWork = & *edges[2];
            pArea = & box_right;
        }

        pWork->mask(tgt, named_color::fetch(YGFX_ORANGE), 255);
        IMG.save("img-final.png", tgt, 0);

        //______________________________________________________________________
        //
        // Now we are working...
        //______________________________________________________________________
        


    }


}
YOCTO_PROGRAM_END()
