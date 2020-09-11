#ifdef POLYPHYLLON
#include <diaplexis/polyphyllon/spacetime.h>
#else
#include <diaplexis/monophyllon/spacetime.h>
#endif

int main(int argc,char** argv)
{
  if (argc != 2) {
    std::cerr << "Usage: ./euplecton <parameter file>" << std::endl;
    return 1;
  }
  std::string filename(argv[1]);

  DIAPLEXIS::Spacetime cosmos(filename);
  cosmos.evolve();

  return 0;
}
