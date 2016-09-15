#include "yocto/program.hpp"
#include "yocto/gfx/image/png.hpp"
#include "yocto/gfx/image/jpeg.hpp"
#include "yocto/gfx/image/tiff.hpp"

using namespace yocto;
using namespace gfx;

YOCTO_PROGRAM_START()
{
    YOCTO_GFX_DECL_FORMAT(png);
    YOCTO_GFX_DECL_FORMAT(jpeg);
    YOCTO_GFX_DECL_FORMAT(tiff);

    imageIO &IMG = image::instance();

    if(argc>1)
    {
        const string  filename = argv[1];
        const pixmap3 origin( IMG.load3(filename,NULL) );
        IMG.save("img.png",origin,NULL);
    }

}
YOCTO_PROGRAM_END()
