// q_warp_affine_tps.h
// warp pointset and image based on given matched pairs
// by Lei Qu
// 2010-03-22

#ifndef __Q_WARP_AFFINE_TPS_H__
#define __Q_WARP_AFFINE_TPS_H__


#include <vector>
using namespace std;
#define WANT_STREAM
#include "newmatap.h"
#include "newmatio.h"
#include "basic_surf_objs.h"

class Coord3D_PCM
{
public:
	double x,y,z;
	Coord3D_PCM(double x0,double y0,double z0) {x=x0;y=y0;z=z0;}
	Coord3D_PCM() {x=y=z=0;}
};



//affine image warp
/*bool q_imagewarp_affine(const vector<Coord3D_PCM> &vec_ctlpt_tar,const vector<Coord3D_PCM>  &vec_ctlpt_sub,
		const unsigned char *p_img_sub, const long long *sz_img_sub, const long long *sz_img_affine,
		unsigned char *&p_img_affine);*/
bool q_imagewarp_affine(const vector<Coord3D_PCM> &vec_ctlpt_tar, const vector<Coord3D_PCM>  &vec_ctlpt_sub,
	const unsigned char *p_img_sub, const long long *sz_img_sub, const long long *sz_img_affine,
	unsigned char *&p_img_affine, int type_gpu, int start_block_x, int start_block_y, int start_block_z);

//affine points warp
bool q_ptswarp_affine(const vector<Coord3D_PCM> &vec_ctlpt_tar, const vector<Coord3D_PCM>  &vec_ctlpt_sub,
	vector<Coord3D_PCM> &vec_ctlpt_subtar_affine);


//Read matched-pair index file
//	output vec2D_sub2tar_matchind is a 2D (n*2) vector
//		vec2D_sub2tar_matchind[i][0]: sub index of i-th matched pair
//		vec2D_sub2tar_matchind[i][1]: tar index of i-th matched pair
//bool q_readMatchInd_file(const QString qs_filename,vector< vector<long> > &vec2D_sub2tar_matchind);


//centrilize and scale the point set ���Ļ������ŵ㼯
bool q_normalize_points_2D(const vector<Coord3D_PCM> vec_input,vector<Coord3D_PCM> &vec_output,Matrix &x3x3_normalize);
bool q_normalize_points_3D(const vector<Coord3D_PCM> vec_input,vector<Coord3D_PCM> &vec_output,Matrix &x4x4_normalize);
//compute the affine matraix
//B=T*A
bool q_affine_compute_affinmatrix_2D(const vector<Coord3D_PCM> &arr_A,const vector<Coord3D_PCM> &arr_B,Matrix &x3x3_affinematrix);
bool q_affine_compute_affinmatrix_3D(const vector<Coord3D_PCM> &arr_A,const vector<Coord3D_PCM> &arr_B,Matrix &x4x4_affinematrix);

//CYF2022.6.14
//changelog:2022.6.16
//changelog:2022.6.21
//bool read_2Dtif_BLOCK(string * &path_tif, unsigned char *&p_img_read_block, unsigned char ****&p_img_read_4d_block, long long * &sz_img_read_block, const long long *sz_img_ori,
//	long long * &sz_img_read_size, long long per_sz_block, long long * &sz_img_read_end_block);
//CYF2022.6.15
//�Ѳ�ֵ���ַŵ�һ��������
//bool q_imagewarp_affine_chazhi_CYF(int type_gpu, const long long *sz_img_output, const Matrix &x4x4_affinematrix, long long x_offset, long long y_offset,
//	long long z_offset, const long long *sz_img_sub, unsigned char ****&p_img_sub_4d, unsigned char ****&p_img_sub2tar_4d);
//change2022.6.16
//change2022.6.21
bool q_imagewarp_affine_chazhi_CYF(int type_gpu, const long long *sz_img_output, const Matrix &x4x4_affinematrix, long long x_offset, long long y_offset, long long z_offset, const long long *sz_img_sub,
	unsigned char ****&p_img_sub_4d, unsigned char ****&p_img_sub2tar_4d, const long long *sz_img_read_block, const long long *sz_img_read_size, const long long *sz_img_sort, const int cal_mode,
	const unsigned char *p_img_sub, unsigned char *p_img_affine);
//���ڷֿ�warp�ĺ���2022.6.16
bool q_imagewarp_affine_cyf(const QList<ImageMarker> &ql_marker_tar, const QList<ImageMarker> &ql_marker_sub, string *path_TIF, long long sz_img_ori[4], string TeraWarp_save_base, int type_gpu);
//���ڽ������ķֿ�warp����
bool q_imagewarp_affine_cyf_downsample(const QList<ImageMarker> &ql_marker_tar, const QList<ImageMarker> &ql_marker_sub, string *path_TIF, long long sz_img_ori[4], string TeraWarp_save_base, int type_gpu);
// all image warp
bool q_imagewarp_affine_image(const QList<ImageMarker> &ql_marker_tar, const QList<ImageMarker> &ql_marker_sub, unsigned char *p_img_sub, long long *sz_img_sub, string TeraWarp_save_base, int type_gpu);

#endif

