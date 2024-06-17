# OpenCV installation on MAC M1
Last update: 17 / 06 / 2024

OpenCV version: 4.x

Link to the installation instructions: https://docs.opencv.org/4.x/d7/d9f/tutorial_linux_install.html 

Kerenel version:

`Darwin Kernel Version 22.6.0: Wed Jul  5 22:22:52 PDT 2023; root:xnu-8796.141.3~6/RELEASE_ARM64_T8103 arm64`

## OpenCV built from source with the C++20 standard

Create install directory:

`mkdir opencv-install`

Clone the code:

`git clone --branch 4.x https://github.com/opencv/opencv.git`

Create build directory:

`cd opencv && mkdir build && cd build`

`cmake -DCMAKE_CXX_STANDARD=20 -DCMAKE_INSTALL_PREFIX=$(cd ../../opencv-install && pwd) ../`

Build the code:

`cmake --build . --target install`

Create `OPENCVSYS` variable, adding this line to your `.profile`:

`export OPENCVSYS="/path/to/opencv-install/"`
