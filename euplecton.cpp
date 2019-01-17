#include <boost/timer/timer.hpp>
#include "diaplexis/spacetime.h"

int main(int argc,char** argv)
{
  if (argc != 2) {
    std::cerr << "Usage: ./euplecton <parameter file>" << std::endl;
    return 1;
  }
  std::string filename(argv[1]);
  boost::timer::auto_cpu_timer timer(3);

  DIAPLEXIS::Spacetime cosmos(filename);
  cosmos.evolve();

  return 0;
}
