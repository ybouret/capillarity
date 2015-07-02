#include "window.hpp"

Window:: ~Window() throw() {}

Window:: Window(const threading::context &ctx,  DataFile &df ) :
threading::window(ctx, df.N, 1 ),
data(df)
{
}
