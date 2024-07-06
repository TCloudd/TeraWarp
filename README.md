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
- [License](#license)

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

3. Install OpenCV 3.1.0. 

4. Download this project to (vaa3d project path)/vaa3d_tools/hackathon and edit `“plugin_Terawarp.pro”`, replace `CUDA_DIR` and `OpenCV_DIR` with the paths where CUDA and OpenCV are located on the current system.

5. Open a terminal window in the current project path and run the following command:
```bash
qmake
make
```



# Run TeraWarp





# License




