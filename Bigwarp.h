#include <QtGui>

//#include "q_warp_affine_tps.h"
#include "basic_memory.cpp"//note: should not include .h file, since they are template functions
#include <string>
#include <sstream>
#include "stackutil.h"
#include "dir_make.h"
#include "RawFmtMngr.h"
#include <vector>
#include <algorithm>

void MkdirWithPath(const string& FolderPath);

bool read_2Dtif_BLOCK(string * &path_tif, unsigned char *&p_img_read_block, unsigned char ****&p_img_read_4d_block, long long * &sz_img_read_block, const long long *sz_img_ori,
	long long * &sz_img_read_size, long long per_sz_block, long long * &sz_img_read_end_block);



bool conversion_point(const long long *warp_pos, long long * &read_pos, const Matrix &x4x4_affinematrix, const int type_mode);

/*
 *ֱ�Ӷ�ȡһ���ȫ������,ȱ���ǵ���һ��XY�ϴ�ʱ��ռ�úܴ��ڴ��޷�����
*/

bool q_imagewarp_stitch_2DRAW_new(const char* imgSrcFile, string imgSrcFile_save, const long long *sz_img_ori, long long per_block_size, const long long *sz_offset_x, const long long *sz_offset_y,
	const long long Z_Slice_control);
/*
2022.10.18
ʹ��Z_Slice_controlȥ����ÿһ���ȡ�Ĳ���
*/

void bubbleSort(vector<long long> &q);


bool loadImage_cyf(char imgSrcFile[], unsigned char *& data1d, V3DLONG * &sz, int & datatype, int start_x, int start_y, int end_x, int end_y);//220612_CYF

//bool loadImage_PJ2DTIF(char imgSrcFile[], unsigned char *& data1d, V3DLONG * &sz, int & datatype, long long Z_slice);//20220719_CYF

bool saveImage_new(const char filename[], const unsigned char * data1d, const V3DLONG * sz, const int datatype);//20220912_CYF(������ע��)

int saveStack2Raw_cyf(const char * filename, const unsigned char * img, const V3DLONG * sz, int datatype);

bool loadImage_PJ(string *&path_TIF, int block_index, int Y_SIZE, int X_SIZE, unsigned char ****p_img_sub_all_4d, const long long *sz_offset_x, const long long *sz_offset_y);

bool saveImage_PJ(long long Z_Slice_max, int block_index, string imgSrcFile_save, long long per_block_size, unsigned char ****p_img_sub_all_4d, V3DLONG *sz_img_sub_save_2d);

bool loadImage_PJ_control(long long Z_Slice_cur, int Z_s, long long Z_Slice_control, string *&path_TIF, int block_index, int Y_SIZE, int X_SIZE, unsigned char ****p_img_sub_all_4d, const long long *sz_offset_x, const long long *sz_offset_y);

bool saveImage_PJ_control(long long Z_Slice_cur, int Z_s, long long Z_Slice_control, int block_index, string imgSrcFile_save, long long per_block_size, unsigned char ****p_img_sub_all_4d, V3DLONG *sz_img_sub_save_2d);

