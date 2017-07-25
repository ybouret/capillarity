#include "yocto/program.hpp"
#include "yocto/fs/local-fs.hpp"
#include "yocto/ptr/auto.hpp"

using namespace yocto;

YOCTO_PROGRAM_START()
{

    vfs &fs = local_fs::instance();

    for(int iarg = 1; iarg < argc; ++iarg )
    {
        string dirname = argv[iarg];
        (void)vfs::as_directory(dirname);
        std::cerr << "processing " << dirname << std::endl;
        auto_ptr<vfs::scanner> scan( fs.new_scanner( dirname ) );
        for( const vfs::entry *ep = scan->next(); ep ; ep = scan->next() )
        {
            //std::cerr << '\t' << '<' << ep->attr << '>' << ep->base_name << " [" << ep->path << "]" << std::endl;
            if( ep->is_regular() )
            {
                const string ext = ep->extension;
                if( "txt" == ext )
                {
                    std::cerr << "\t[" << ep->path << "]" << std::endl;
                }
            }
        }

    }

}
YOCTO_PROGRAM_END()
