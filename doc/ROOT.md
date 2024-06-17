# ROOT installation on MAC M1
Last update: 03 / 06 / 2024

Root version: v6-32-00

Kernel version:

`Darwin Kernel Version 22.6.0: Wed Jul  5 22:22:52 PDT 2023; root:xnu-8796.141.3~6/RELEASE_ARM64_T8103 arm64`

## ROOT built from source with C++20 standard

```
brew uninstall llvm && brew install llvm@16
brew install zstd
brew install libedit
```

Create source, build, and install directories:

`mkdir root_src && mkdir root-built && mkdir root-install`

Clone the code:

`git clone --branch v6-32-00 --depth 1 https://github.com/root-project/root.git root_src`

Added following lines to `core/CMakeFileLists.txt` inside the `if (runtime_cxxmodules)` loop:

```
  list(APPEND core_implicit_modules "-mByproduct" "std_vector")
  list(APPEND core_implicit_modules "-mByproduct" "std_functional")
  list(APPEND core_implicit_modules "-mByproduct" "std_sstream")
  list(APPEND core_implicit_modules "-mByproduct" "std_numeric")
  list(APPEND core_implicit_modules "-mByproduct" "std_fstream")
  list(APPEND core_implicit_modules "-mByproduct" "std_iostream")
  list(APPEND core_implicit_modules "-mByproduct" "std_map")
  list(APPEND core_implicit_modules "-mByproduct" "std_unordered_set")
  list(APPEND core_implicit_modules "-mByproduct" "std_chrono")
  list(APPEND core_implicit_modules "-mByproduct" "std_complex")
  list(APPEND core_implicit_modules "-mByproduct" "std_forward_list")
  list(APPEND core_implicit_modules "-mByproduct" "std_list")
  list(APPEND core_implicit_modules "-mByproduct" "std_condition_variable")
  list(APPEND core_implicit_modules "-mByproduct" "std_thread")
  list(APPEND core_implicit_modules "-mByproduct" "std_queue")
  list(APPEND core_implicit_modules "-mByproduct" "std_span")
```


Setup cmake:

`cd root-build` 

`cmake -DCMAKE_CXX_STANDARD=20 -DCMAKE_INSTALL_PREFIX=$(cd ../root-install && pwd) $(cd ../root_src && pwd) -DCMAKE_CXX_COMPILER=$(which clang++) -DCMAKE_C_COMPILER=$(which clang) -DPython3_ROOT_DIR=$(which python) -Dbuiltin_glew=ON -Dcocoa=OFF -Dbuiltin_tbb=OFF -Dx11=ON -Dpyroot=OFF`

Build parallely on `N` cores:

`cmake --build . --target install - jN`

Add the following line to your `~/.profile`:

`source /path/to/your/root-install/bin/thisroot.sh`
