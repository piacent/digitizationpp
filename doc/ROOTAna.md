# ROOTAna installation on MAC M1
Last update: 17 / 06 / 2024

Rootana version: 17 / 06 / 2024 (master branch)

Link to the repository: https://bitbucket.org/tmidas/rootana/src/master/

Kerenel version:

`Darwin Kernel Version 22.6.0: Wed Jul  5 22:22:52 PDT 2023; root:xnu-8796.141.3~6/RELEASE_ARM64_T8103 arm64`

## ROOTAna built from source with the C++20 standard


Create source and install directories:

`mkdir rootana`

Clone the code:

`git clone https://bitbucket.org/tmidas/rootana.git rootana`

Edit the Makefile and change the following line:

`CXXFLAGS = -std=c++11 -O2 -g -Wall -Wuninitialized -I./include`

into 

`CXXFLAGS = -std=c++20 -O2 -g -Wall -Wuninitialized -I./include`

Build rootana:

`cd rootana`

`make`

Create `ROOTANASYS` variable, adding this line to your `.profile`:

`export ROOTANASYS=/path/to/rootana/`
