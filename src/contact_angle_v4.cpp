#include "yocto/program.hpp"

#include "yocto/gfx/image.hpp"
#include "yocto/gfx/image/jpeg.hpp"
#include "yocto/gfx/image/png.hpp"
#include "yocto/gfx/image/tiff.hpp"

#include "yocto/gfx/ops/edges.hpp"
#include "yocto/gfx/ops/filter.hpp"
#include "yocto/string/conv.hpp"

using namespace yocto;
using namespace gfx;
using namespace math;


static unit_t y_limit = 0;
static inline bool is_bad_vertex(const vertex &v) throw()
{
    return v.y>=y_limit;
}

static inline
int compare_by_decreasing_height( const particle::ptr &lhs, const particle::ptr &rhs )
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
}


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

        // sort by Y extension
        quicksort(edges, compare_by_decreasing_height);
        particle::ptr A( edges[1] );
        particle::ptr B( edges[2] );
        if( A->center().x >= B->center().x )
        {
            bswap(A,B);
        }
        // A is LEFT part, B is RIGHT part
        tgt.copy(source);
        A->mask(tgt, named_color::fetch(YGFX_RED),   255);
        B->mask(tgt, named_color::fetch(YGFX_GREEN), 255);
        IMG.save("img-wksp.png", tgt, 0 );

        // analyze shapes
        Vertices A_shape;
        analyze_particle(A_shape,*A,1);

        Vertices B_shape;
        analyze_particle(B_shape,*B,-1);

        tgt.copy(source);
        for(size_t i=1;i<=A_shape.size();++i)
        {
            tgt[A_shape[i]] = named_color::fetch(YGFX_RED);
        }
        for(size_t i=1;i<=B_shape.size();++i)
        {
            tgt[B_shape[i]] = named_color::fetch(YGFX_GREEN);
        }
        IMG.save("img-shape.png", tgt, 0 );

    }

}
YOCTO_PROGRAM_END()

