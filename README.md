# digitizationpp

** work in progress **

Digitization code in C++

## Dependencies

* ROOT [compiled with the C++17 standard]


## Installation

`git clone https://github.com/CYGNUS-RD/digitizationpp.git`

`cd digitizationpp`

`export CXX="/path/to/your/c++17/compiler"`

`mkdir build-dir && cd build-dir`

`cmake ..`

`cmake --build .`

## Suggested usage

Put all the MC `.root` files you want to digitize in the `input/` folder, then:

`cd build-dir`

`./digitizationpp ../config/ConfigFile_new.txt`

The output file will be saved in the `OutDir/` folder.


