Rvalue=60
PROG=../bin/cbridge4
SWAP=true  $PROG vierge_tirage100mu_R81mm.txt*   R0=$Rvalue
SWAP=true  $PROG vierge_tirage20mu_R81mm.txt*    R0=$Rvalue
SWAP=true  $PROG vierge_tirage50mu_R81mm.txt*    R0=$Rvalue
SWAP=true  $PROG vierge_tirage5mu_R81mm.txt*     R0=$Rvalue
SWAP=false $PROG vierge_tirage5mu_R81mm_bis.txt* R0=$Rvalue

