#include "yocto/program.hpp"
#include "yocto/gfx/image/png.hpp"
#include "yocto/gfx/image/jpeg.hpp"
#include "yocto/gfx/image/tiff.hpp"
#include "yocto/gfx/ops/stencil.hpp"
#include "yocto/gfx/ops/edges.hpp"
#include "yocto/gfx/ops/filter.hpp"
#include "yocto/gfx/draw/stencil.hpp"
#include "yocto/container/matrix.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/alg/shapes2d.hpp"
#include "yocto/math/fit/glsf-spec.hpp"
#include "yocto/math/fcn/zfind.hpp"
#include "yocto/fs/local-fs.hpp"
#include "yocto/lingua/pattern/matcher.hpp"
#include "yocto/lingua/pattern/regexp.hpp"
#include "yocto/lua/lua-state.hpp"

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
    double         deltaY;

    Geometry() : center(), radius(0), aorg(3), deltaY(0)
    {
    }

    ~Geometry() throw()
    {
    }

    vertex Coord( const double alpha ) const throw()
    {
        const double rr = radius+_GLS::Polynomial<double>::Eval(alpha,aorg);
        return vertex( unit_t(center.x + rr * sin(alpha)),
                      unit_t(center.y - rr * cos(alpha))
                      );
    }


    double ComputeZero( const double alpha ) throw()
    {
        return (radius+_GLS::Polynomial<double>::Eval(alpha,aorg)) * cos(alpha) - deltaY;
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Geometry);
};

// To fit the ridge intensity profile
class Ridge
{
public:
    static const size_t NVAR = 4;

    Ridge()
    {
    }

    ~Ridge()
    {
    }

    double Compute( const double alpha, const array<double> &param )
    {
        assert(NVAR==param.size());
        const double level = param[1];
        const double split = param[2];
        const double slope = param[3];
        const double curv  = param[4];
        if(alpha<=split)
        {
            return level;
        }
        else
        {
            const double da = alpha-split;
            return level + da * slope + da*da * curv;
        }
    }


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Ridge);
};


static double DeltaPixel = 2;
static double SweepAngle = 30;
static float  GaussSigma = 1.4f;

static inline
double process_file( const string &filename, const string &side )
{

    imageIO &IMG = image::instance();

    // preparing data base infos
    const string base_name = vfs::get_base_name(filename);
    string       root_name = base_name;
    vfs::remove_extension(root_name);
    std::cerr << "base_name=" << base_name << std::endl;
    std::cerr << "root_name=" << root_name << std::endl;
    string       dir_name  = vfs::get_file_dir(filename);
    std::cerr << "dir_name =" << dir_name << std::endl;
    for(size_t i=0;i<dir_name.size()-1;++i)
    {
        if(dir_name[i]=='/')
        {
            dir_name[i] = '_';
        }
    }
    const string db_name = "db/" + dir_name;
    std::cerr << "db_name  =" << db_name << std::endl;

    vfs & fs = local_fs::instance();
    fs.create_sub_dir(db_name);


    bool            is_right = true;
    bool            is_left  = false;

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
    const pixmap3 origin( IMG.load3(filename,NULL) );
    if(is_left)
    {
        std::cerr << "-- Flipping Picture to get data on the right side" << std::endl;
        ((pixmap3 &)origin).flip_horizontal();
    }
    //IMG.save("img-origin.png",origin,NULL);


    //______________________________________________________________________
    //
    // making a grayscale picture
    //______________________________________________________________________
    const pixmapf img0(origin,RGB::to_float,origin);

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
    if(GaussSigma>0)
    {
        std::cerr << "-- Gauss " << GaussSigma << std::endl;
        stencil_gauss g5(2,GaussSigma);
        g5.apply(img,img0,xps);
    }
    else
    {
        img.copy(img0);
    }
    //IMG.save("img-filter.png",img,NULL);

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
    //IMG.save("img-grad.png", ED, 0 );
    //IMG.save("img-tags.png",tags,tags.colors,NULL);
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
    //IMG.save("img-edges.png", tgt, 0);

    //______________________________________________________________________
    //
    // close shapes
    //______________________________________________________________________
    std::cerr << "-- Closing..." << std::endl;
    Filter<float> F(w,h);
    F.Close(ED,xps);
    //IMG.save("img-close.png", ED, 0);

    //______________________________________________________________________
    //
    // detect final blobs
    //______________________________________________________________________
    std::cerr << "-- Final Blobs..." << std::endl;
    tags.build(ED,8);
    edges.load(tags);
    //IMG.save("img-blobs.png", tags, tags.colors, 0);
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
    const string img_final = db_name + root_name + "-final.png";
    //IMG.save(img_final, tgt, 0);


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
    //IMG.save(img_final, tgt, 0);


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

    //IMG.save(img_final, tgt, 0);


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
    //IMG.save(img_final, tgt, 0);

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

    if(false)
    {
        const string shape_name = db_name + root_name + "-shape.dat";
        ios::wcstream fp(shape_name);
        for(size_t i=1;i<=N;++i)
        {
            fp("%g %g %g %g %g\n", X[i], Y[i], Rad2Deg(A[i]), R[i], R[i] - radius);
        }
    }


    co_qsort(A,R);

    // TODO: define delta in pixels and sweep angle ??
    const double Aini  = A[1];
    const double Aend  = min_of<double>( Aini + Deg2Rad(SweepAngle), Deg2Rad(90.0) );

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

        if(yy>=y_low+DeltaPixel&&aa<=Aend)
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

    if(false)
    {
        const string shape_fit_name = db_name + root_name + "-shape.dat";
        ios::wcstream fp(shape_fit_name);
        for(size_t i=1;i<=NN;++i)
        {
            fp("%g %g %g %g %g %g\n", XX[i], YY[i], Rad2Deg(AA[i]), RR[i], RR[i] - radius, RF[i]);
        }
    }


    std::cerr << "y_low=" << y_low << std::endl;


    //______________________________________________________________________
    //
    // Now we have an approximation of the shape...let's scan the line!
    // Let us find the starting angle...then scan down...
    //______________________________________________________________________
    Geom.deltaY = center.y - y_low;
    zfind<double>             solver(1e-5);
    numeric<double>::function zfn( &Geom, & Geometry::ComputeZero );
    triplet<double>           zAlpha = { 0, 0, Aend };
    triplet<double>           zValue = { zfn(zAlpha.a), 0, zfn(zAlpha.c) };
    if(zValue.a*zValue.c>0)
    {
        throw exception("Cannot find interception, corrupted picture?");
    }


    const double alpha0 = solver.run(zfn, zAlpha, zValue);
    std::cerr << "alpha0=" << Rad2Deg(alpha0) << std::endl;

    //alpha0 += Deg2Rad(2.0);

    const double alphaScan = Deg2Rad(5.0);
    const double alpha_min = max_of<double>(0,alpha0-alphaScan);
    const double alpha_max = min_of<double>(Aend,alpha0+alphaScan);
    const double alpha_del = alpha_max - alpha_min;


    //______________________________________________________________________
    //
    // Let's find all the coordinates around the intersection
    //______________________________________________________________________
    const double da      = 1.0/(2.0*radius); //should have half a pixel of difference...
    const size_t na     = size_t(ceil(alpha_del/da)+1);
    std::cerr << "na=" << na << std::endl;
    vector<vertex>  va(na,as_capacity);
    vector<double>  ak(na,as_capacity);

    for(size_t i=0;i<=na;++i)
    {
        const double alpha = alpha_min + (i*alpha_del)/double(na);
        const vertex pos   = Geom.Coord(alpha);
        if(tgt.has(pos))
        {
            tgt[pos] = named_color::fetch(YGFX_DEEP_PINK);
            bool append = true;
            for(size_t j=va.size();j>0;--j)
            {
                if( pos == va[j] )
                {
                    append = false;
                    break;
                }
            }
            if(append)
            {
                va.push_back(pos);
                ak.push_back(alpha);
            }
        }
    }
    //IMG.save(img_final, tgt, 0);
    std::cerr << "#vtx=" << va.size() << "/" << na << std::endl;

    const size_t ns = va.size();



    //______________________________________________________________________
    //
    // Let's find the change in intensity
    //______________________________________________________________________
    vector<double>       ridgeI(ns); // ridge intensity
    vector<double>       ridgeF(ns); // ridge fit
    pixmapf             &levelI = img; // filtered

    for(size_t i=1;i<=ns;++i)
    {
        ridgeI[i] = levelI[ va[i] ];
    }

    GLS<double>::Samples ridgeSamples;
    ridgeSamples.append(ak,ridgeI,ridgeF);


    vector<double> ridgeA(Ridge::NVAR);
    vector<bool>   ridgeU(Ridge::NVAR,true);
    vector<double> ridgeE(Ridge::NVAR);

    ridgeSamples.prepare(Ridge::NVAR);

    Ridge                 ridge;
    GLS<double>::Function ridgeFit( &ridge, & Ridge::Compute );


    // initialize variables
    ridgeA[1] = ridgeI[1]; // level
    ridgeA[2] = alpha0;    // split
    ridgeA[3] = (ridgeI[ns] - ridgeA[1])/alphaScan; // slope
    ridgeA[4] = 0; // curv

    ridgeU[4] = false;

    if(false)
    {
        ios::wcstream fp("scan.dat");
        for(size_t i=1;i<=ns;++i)
        {
            fp("%g %g %g\n", ak[i], ridgeI[i], ridgeFit(ak[i],ridgeA));
        }
    }


    //______________________________________________________________________
    //
    // fit level 1, linear part only
    //______________________________________________________________________

    if( ! ridgeSamples.fit_with(ridgeFit, ridgeA, ridgeU, ridgeE) )
    {
        throw exception("couldn't fit level-1");
    }

    std::cerr << "RIDGE 1/2: " << std::endl;
    GLS<double>::display(std::cerr, ridgeA, ridgeE);
    if(false)
    {
        const string scanfit_name = db_name + root_name + "-scanfit1.dat";
        ios::wcstream fp(scanfit_name);
        for(size_t i=1;i<=ns;++i)
        {
            fp("%g %g %g\n", ak[i], ridgeI[i], ridgeF[i]);
        }
    }

    //______________________________________________________________________
    //
    // fit level 2, with quadratic part
    //______________________________________________________________________
    ridgeU[4] = true;

    if( ! ridgeSamples.fit_with(ridgeFit, ridgeA, ridgeU, ridgeE) )
    {
        throw exception("couldn't fit level-2");
    }

    std::cerr << "RIDGE 2/2: " << std::endl;
    GLS<double>::display(std::cerr, ridgeA, ridgeE);

    if(false)
    {
        const string scanfit_name = db_name + root_name + "-scanfit2.dat";
        ios::wcstream fp(scanfit_name);
        for(size_t i=1;i<=ns;++i)
        {
            fp("%g %g %g\n", ak[i], ridgeI[i], ridgeF[i]);
        }
    }

    //______________________________________________________________________
    //
    // ok, we found the intersection angle
    //______________________________________________________________________
    const double alphaI = ridgeA[2];

    //______________________________________________________________________
    //
    // so we compute the intersection coordinate
    //______________________________________________________________________
    const vertex Q      = Geom.Coord(alphaI);

    //______________________________________________________________________
    //
    // then we compute the derivative of the extrapolation
    //______________________________________________________________________
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
    const double beta  = Atan(num/den);
    const double theta = numeric<double>::pi - beta;
    std::cerr << "theta=" << Rad2Deg(theta) << std::endl;


    // final drawing
    tgt.copy(origin);
    vector<vertex> sten;


    const unit_t thick = 2;

    // horiontal
    sten.free();
    collect_disk(sten, Q.x, Q.y, thick);
    extrude_stencil(sten, vertex(p2->vtx.x-Q.x,0) );

    draw_stencil(sten, tgt, vertex(0,0), named_color::fetch(YGFX_BLUE_LIGHT), 0xff);

    // profile
    sten.free();
    {
        double aa = alpha_min;
        while(aa<Aend)
        {
            const vertex pos = Geom.Coord(aa);
            sten.push_back(pos);
            aa += da;
        }
        clean_stencil(sten);
        extrude_stencil(sten, vertex(-2,0) );
        draw_stencil(sten, tgt, vertex(0,0), named_color::fetch(YGFX_YELLOW), 0xff);
    }



    //tangent
    sten.free();
    collect_disk(sten, Q.x, Q.y, thick);
    clean_stencil(sten);
    const double length = radius/2.0;
    extrude_stencil(sten, vertex(length*cos(beta),length*sin(beta) ) );
    draw_stencil(sten, tgt, vertex(0,0), named_color::fetch(YGFX_RED), 200);

    //intersection
    //draw_disk(tgt, Q.x, Q.y, 6, named_color::fetch(YGFX_MAGENTA), 127 );

#if 0
    if(false)
    {
        const double length = radius/5.0;
        draw_line(tgt,
                  Q.x,Q.y,
                  unit_t(Q.x+length*cos(beta)),unit_t(Q.y+length*sin(beta)),
                  named_color::fetch(YGFX_FIREBRICK), 0xff
                  );
    }
#endif

    const string img_angle = db_name+root_name+"-angle.png";
    IMG.save(img_angle, tgt, 0);


    const string info_name = db_name+root_name+"-theta.txt";
    {
        ios::wcstream fp(info_name);
        fp("%.3f\n",Rad2Deg(theta));
    }

    return Rad2Deg(theta);
}


#define __SHOW(VAR) std::cerr << #VAR " = " << VAR << std::endl

YOCTO_PROGRAM_START()
{
    YOCTO_GFX_DECL_FORMAT(png);
    YOCTO_GFX_DECL_FORMAT(jpeg);
    YOCTO_GFX_DECL_FORMAT(tiff);

    std::cerr << "-- Reading Parameters" << std::endl;
    Lua::State L;
    L.DoFile("angle.lua");

    DeltaPixel = max_of<double>(0,L.Get<double>("DeltaPixel" ));
    SweepAngle = clamp<double>(5,L.Get<double>("SweepAngle" ),90);
    GaussSigma = max_of<float>(0,L.Get<double>("GaussSigma") );

    __SHOW(DeltaPixel);
    __SHOW(SweepAngle);
    __SHOW(GaussSigma);


    if(argc<=2)
    {
        throw exception("usage: %s folder [left|right]", program);
    }

    const string folder = argv[1];
    const string side   = argv[2];

    vfs &fs = local_fs::instance();

    //-- extract files to process
    lingua::matcher        pm("theta(A|R)");
    auto_ptr<vfs::scanner> scan( fs.new_scanner( folder )     );
    for( const vfs::entry *ep = scan->next(); ep ; ep = scan->next() )
    {
        //std::cerr << '<' << ep->attr << '>' << ep->base_name << " [" << ep->path << "]" << std::endl;
        const string bname = ep->base_name;
        if(pm.partial_match(bname))
        {
            std::cerr << "will process " << ep->path << std::endl;
            const string filename = ep->path;
            const double theta    = process_file(filename,side);
            std::cerr << "theta=" << theta << std::endl;
        }
    }
    
}
YOCTO_PROGRAM_END()
