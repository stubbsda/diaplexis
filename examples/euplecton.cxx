#ifdef POLYPHYLLON
#include <diaplexis/polyphyllon/spacetime.h>
#else
#include <diaplexis/monophyllon/spacetime.h>
#endif

int main(int argc,char** argv)
{
  if (argc != 3) {
    std::cerr << "Usage: ./euplecton --type=[discrete,continuous] <parameter file>" << std::endl;
    return 1;
  }

  // Parse the value of argv[1], which should have the form '--type=discrete' or '--type=continuous'
  bool discrete = true;
  std::string pname,pvalue,ftype(argv[1]);
  unsigned int eq = ftype.find('=');
  if (eq == std::string::npos) {
    std::cerr << "Usage: ./euplecton --type=[discrete,continuous] <parameter file>" << std::endl;
    return 1;
  }
  pname = ftype.substr(0,eq);
  pvalue = ftype.substr(eq+1,ftype.length());
  if (pname == "--type") {
    if (pvalue == "discrete") {
      discrete = true;
    }
    else if (pvalue == "continuous") {
      discrete = false;
    }
    else {
      std::cerr << "Usage: ./euplecton --type=[discrete,continuous] <parameter file>" << std::endl;
      return 1;
    }
  }
  else {
    std::cerr << "Usage: ./euplecton --type=[discrete,continuous] <parameter file>" << std::endl;
    return 1;
  }

  std::string filename(argv[2]);
  if (discrete) {
    DIAPLEXIS::Spacetime<SYNARMOSMA::UINT64,SYNARMOSMA::INT64> cosmos(filename);
    cosmos.evolve();
  }
  else {
    DIAPLEXIS::Spacetime<double,double> cosmos(filename);
    cosmos.evolve();
  }

  return 0;
}
