#include "diaplexis/spacetime.h"

int main(int argc,char** argv)
{
  DIAPLEXIS::Spacetime cosmos(argv[1]);
  cosmos.evolve();

  return 0;
}
