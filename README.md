# digitizationpp

** work in progress **

Digitization code in C++

## Dependencies

* ROOT [compiled with the C++17 standard]
* ROOTANA [https://bitbucket.org/tmidas/rootana/src/master/]
* OPENCV [https://docs.opencv.org/4.x/d7/d9f/tutorial_linux_install.html]

Before compiling, set the variables ROOTANASYS and OPENCVSYS in your environment:

`export ROOTANASYS="/path/to/rootana/installation/"`
`export OPENCVSYS="/path/to/opencv/installation/"`


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


