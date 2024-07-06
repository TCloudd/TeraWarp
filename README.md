# TeraWarp: Fast Terabyte-Scale Volumetric Image Warping Tool
![Image](figure1.png)

# Contents
- [System requirements](#System-requirments)
  - [Required device](#--Required-device)
  - [Optional device](#--Optional-device)
- [Installation](#Installation)
  - [Windows](#--Windows)
  - [Linux](#--Linux)
- [Run TeraWarp](#Run-TeraWarp)
  - [Interface](#--Interface)
  - [Commandline](#--Commandline)

# System requirments
TeraWarp works on desktop computers, for terabyte-scale data, our recommended computer configurations for efficiency gains are as follows:

## - Required device
RAM: 64GB or larger

CPU: 2.3GHz Intel Xeon Gold or better

Hard disk: Three times the size of the original image or larger

## - Optional device
GPU: NVIDIA GeForce RTX 3090 GPU or better

(Note: The time required for image warping is closely related to GPU, and we recommend using GPU to process the image.)

# Installation
TeraWarp is accessible as an open-source plugin for the Vaa3D platform. The tool supports Windows 64-bit systems and popular Linux distributions. The basis for successfully compiling TeraWarp is that you have installed Vaa3D.


## - Windows

1. Download Qt 4.8.6 (qt-opensource-windows-x86-vs2010-4.8.6.exe) from [QT official website](https://download.qt.io/archive/qt/4.8/4.8.6/). We recommend using the Visual Studio method of installation, here we are using Visual Studio 2013.

2. Download CUDA Toolkit 11.1 from [NVIDIA CUDA Toolkit download page](https://developer.nvidia.com/cuda-downloads)

3. Download OpenCV 3.1.0.

4. Download this project to (vaa3d project path)/vaa3d_tools/hackathon and edit `“plugin_Terawarp.pro”`, replace `CUDA_DIR` and `OpenCV_DIR` with the paths where CUDA and OpenCV are located on the current system.

5. Open command line Terminal in the Visual Studio 2013 folder. This will give you a command line window. Then run the following command:
```bash
qmake
nmake -f Makefile.Release
```


## - Linux

1. Download Qt 4.8.6. Compile the “Release” version and install it by run the following commands in a Terminal window.
```bash
cd <unzipped_qt_source_code_folder>
./configure
make
sudo make install
```

2. Download CUDA Toolkit 11.1 from [NVIDIA CUDA Toolkit download page](https://developer.nvidia.com/cuda-downloads) and run the installer:
```bash
sudo sh cuda_version_linux.run
```

3. Download OpenCV 3.1.0. 

4. Download this project to (vaa3d project path)/vaa3d_tools/hackathon and edit `“plugin_Terawarp.pro”`, replace `CUDA_DIR` and `OpenCV_DIR` with the paths where CUDA and OpenCV are located on the current system.

5. Open a terminal window in the current project path and run the following command:
```bash
qmake
make
```


# Run TeraWarp
We provide two ways to run TeraWarp, one through the interface and one through the command line.

## - Interface
Open Vaa3D, select `Plug-In -> Tc -> TeraWarp` in the upper left corner of the interface, parameter selection via interface controls.

## - Commandline









