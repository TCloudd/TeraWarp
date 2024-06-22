#ifdef _WIN32
#include "dirent_win.h"
#else
#include <dirent.h>
#endif
//#include "q_warp_affine_tps.h"
//#include "Bigwarp.h"
//void read_lryic(char imgSrcFile[], stringstream &img_path);

bool read_RES_dir(char* &imgSrcFile, string *** &path_Z);

void eachnumber(long long *sz);

bool read_2Dtif_dir(const char* imgSrcFile, string * &path_tif, long long &DEPTH_TIF);//2DTIF_dirent 2022.6.14