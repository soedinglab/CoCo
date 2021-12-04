# CoCo

CoCo is a software suite for different **Co**nsensus **Co**rrection applications using spaced k-mer count profiles of short reads or contigs. CoCo is open source GPL-licensed software implemented in C++.

### Compile from source
To compile CoCo `git`, a `C++/11` capable compiler (e.g. `gcc 4.7+`, `clang 3.5+`) and `cmake` (3.10 or higher) are required. Afterwards, the CoCo binary will be located in the `build/` directory.

      git clone https://github.com/soedinglab/CoCo.git
      cd CoCo
      git submodule update --init
      mkdir build && cd build
      cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
      make -j 4
