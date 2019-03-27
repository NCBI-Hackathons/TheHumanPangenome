#include <cstdlib>
#include <iostream>

#include "ggff.h"

int main(int argc, char **argv) noexcept {
  if (argc != 3) {
    std::cerr << "ggfftools ggff1 ggff2\n";
    return EXIT_FAILURE;
  }

  ggff::GGFFIndex idx("test.g2f2i");
  idx.add_ggff(argv[0]);
  idx.add_ggff(argv[1]);
  return EXIT_SUCCESS;
}
