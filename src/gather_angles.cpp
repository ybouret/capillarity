#include "yocto/program.hpp"
#include "yocto/fs/local-fs.hpp"

using namespace yocto;

YOCTO_PROGRAM_START()
{

    for(int iarg = 1; iarg < argc; ++iarg )
    {
        string dirname = argv[iarg];
        (void)vfs::as_directory(dirname);
        std::cerr << "processing " << dirname << std::endl;
    }

}
YOCTO_PROGRAM_END()
