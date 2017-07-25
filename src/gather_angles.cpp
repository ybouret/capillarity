#include "yocto/program.hpp"
#include "yocto/fs/local-fs.hpp"
#include "yocto/ptr/auto.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/ios/icstream.hpp"
#include "yocto/ios/ocstream.hpp"

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
        vector<string> lines;

        // load content and format lines
        for( const vfs::entry *ep = scan->next(); ep ; ep = scan->next() )
        {
            //std::cerr << '\t' << '<' << ep->attr << '>' << ep->base_name << " [" << ep->path << "]" << std::endl;
            if( ep->is_regular() )
            {
                const string fn  = ep->base_name;
                if(fn=="theta.txt")
                {
                    continue;
                }
                const string ext = ep->extension;
                if( "txt" == ext )
                {
                    string id = ep->base_name;
                    id.trim(10);
                    std::cerr << "\t[" << ep->path << "] -> " << id << std::endl;
                    string line;
                    ios::icstream fp(ep->path);
                    if( !fp.gets(line) )
                    {
                        throw exception("couldn't get line from '%s'", id.c_str() );
                    }
                    std::cerr << "\t\t" << line << std::endl;
                    line += ' ';
                    line += id;

                    lines.push_back(line);
                }
            }
        }


        // prepare output
        const string output = dirname + "theta.txt";
        std::cerr << "saving into " << output << std::endl;
        {
            ios::wcstream fp(output);
            for(size_t i=1;i<=lines.size();++i)
            {
                fp << lines[i] << "\r\n";
            }
        }

    }

}
YOCTO_PROGRAM_END()
