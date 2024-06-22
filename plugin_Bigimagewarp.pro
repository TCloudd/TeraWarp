
TEMPLATE	= lib
CONFIG	+= qt plugin warn_off
#CONFIG	+= x86_64
VAA3DPATH = ../../../../v3d_external
INCLUDEPATH	+= $$VAA3DPATH/v3d_main/basic_c_fun
INCLUDEPATH  += $$VAA3DPATH/v3d_main/common_lib/include
INCLUDEPATH  += $$VAA3DPATH/v3d_main/jba/newmat11
INCLUDEPATH  += $$VAA3DPATH/v3d_main/jba/c++


win32{
QMAKE_CXXFLAGS+=/openmp

LIBS += -L$$VAA3DPATH/v3d_main/common_lib/winlib64 -llibtiff
LIBS += -L$$VAA3DPATH/v3d_main/common_lib/winlib64 -llibnewmat

# cuda
LIBS += -laff
LIBS += -lchazhi
LIBS += -lextendornormal
LIBS += -lGET_A_i
LIBS += -lGETA
LIBS += -lqr
LIBS += -lweiyijisuan
LIBS += -lXNX4C
}



unix{
#CONFIG += console c++11
QMAKE_CXXFLAGS += -fopenmp

LIBS += -L$$VAA3DPATH/v3d_main/common_lib/lib -lv3dtiff
LIBS += -L$$VAA3DPATH/v3d_main/jba/c++ -lv3dnewmat

OTHER_FILES += GET_A_i.cu
OTHER_FILES += aff.cu
OTHER_FILES += aff.cu
OTHER_FILES += chazhi.cu
OTHER_FILES += extendornormal.cu
OTHER_FILES += GETA.cu
OTHER_FILES += weiyijisuan.cu
OTHER_FILES += XNX4C.cu


CUDA_SOURCES += GET_A_i.cu
CUDA_SOURCES += aff.cu
CUDA_SOURCES += aff.cu
CUDA_SOURCES += chazhi.cu
CUDA_SOURCES += extendornormal.cu
CUDA_SOURCES += GETA.cu
CUDA_SOURCES += weiyijisuan.cu
CUDA_SOURCES += XNX4C.cu

CUDA_SDK = "/usr/local/cuda-11.1"
CUDA_DIR = "/usr/local/cuda-11.1"
SYSTEM_NAME = ubuntu
SYSTEM_TYPE = 64
CUDA_ARCH = sm_50
NVCC_OPTIONS = --use_fast_math

INCLUDEPATH += $$CUDA_DIR/include
QMAKE_LIBDIR += $$CUDA_DIR/lib64/

CUDA_OBJECTS_DIR = ./
CUDA_LIBS = -lcudart_static -lcusolver_static -lcublas_static -lcublasLt_static -lcudadevrt -lculibos

CUDA_INC = $$join(INCLUDEPATH,'" -I"','-I"','"')
LIBS += $$CUDA_LIBS


# Configuration of the Cuda compiler
    cuda.input = CUDA_SOURCES
    cuda.output = $$CUDA_OBJECTS_DIR/${QMAKE_FILE_BASE}_cuda.o
    cuda.commands = $$CUDA_DIR/bin/nvcc $$NVCC_OPTIONS $$CUDA_INC $$LIBS -Xcompiler -fPIC --machine $$SYSTEM_TYPE -arch=$$CUDA_ARCH -c -o ${QMAKE_FILE_OUT} ${QMAKE_FILE_NAME}
    cuda.dependency_type = TYPE_C
    QMAKE_EXTRA_COMPILERS += cuda
#Set the additional inclusion directory required for the program to run

LIBS += -fopenmp

}




FORMS         = q_paradialog_bigimagewarp.ui

HEADERS      += $$VAA3DPATH/v3d_main/basic_c_fun/v3d_message.h
HEADERS      += $$VAA3DPATH/v3d_main/basic_c_fun/basic_surf_objs.h
#HEADERS      += $$VAA3DPATH/v3d_main/basic_c_fun/basic_memory.h
HEADERS	     += $$VAA3DPATH/v3d_main/basic_c_fun/stackutil.h
HEADERS      += $$VAA3DPATH/v3d_main/basic_c_fun/mg_image_lib.h
HEADERS      += $$VAA3DPATH/v3d_main/basic_c_fun/mg_utilities.h
HEADERS      += $$VAA3DPATH/v3d_main/jba/c++/jba_mainfunc.h
HEADERS      += $$VAA3DPATH/v3d_main/jba/c++/jba_match_landmarks.h
HEADERS      += q_paradialog_bigimagewarp.h
HEADERS      += plugin_Bigimagewarp.h
HEADERS      += dir_make.h
HEADERS      += dirent_win.h
HEADERS      += RawFmtMngr.h
HEADERS      += IM_config.h
HEADERS      += q_warp_affine_tps.h
HEADERS      += Bigwarp.h
HEADERS      += q_imgwarp_tps_quicksmallmemory.h
HEADERS      += q_littleQuickWarp_common.h


SOURCES      += $$VAA3DPATH/v3d_main/basic_c_fun/v3d_message.cpp
SOURCES      += $$VAA3DPATH/v3d_main/basic_c_fun/stackutil.cpp
SOURCES      += $$VAA3DPATH/v3d_main/basic_c_fun/basic_surf_objs.cpp
#SOURCES      += $$VAA3DPATH/v3d_main//basic_c_fun/basic_memory.cpp
SOURCES      += $$VAA3DPATH/v3d_main/basic_c_fun/mg_image_lib.cpp
SOURCES      += $$VAA3DPATH/v3d_main/basic_c_fun/mg_utilities.cpp
SOURCES      += q_paradialog_bigimagewarp.cpp
SOURCES      += plugin_Bigimagewarp.cpp
SOURCES      += dir_make.cpp
SOURCES      += RawFmtMngr.cpp
SOURCES      += IM_config.cpp
SOURCES      += q_warp_affine_tps.cpp
SOURCES      += Bigwarp.cpp
SOURCES      += q_imgwarp_tps_quicksmallmemory.cpp

TARGET	= $$qtLibraryTarget(Bigimagewarp)
DESTDIR	= $$VAA3DPATH/bin/plugins/BigImageWarp/
