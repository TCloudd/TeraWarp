
#ifndef __Q_IMGWARP_TPS_QUICKSMALLMEMORY_H__
#define __Q_IMGWARP_TPS_QUICKSMALLMEMORY_H__


#include "basic_surf_objs.h"
#include "basic_memory.cpp"
#include "newmatap.h"
#include "newmatio.h"
#include "q_littleQuickWarp_common.h"



//linear interpolate the SubDFBlock to DFBlock
//use 3d or 4d pointer instead of 1d, since generating 3d or 4d pointer from 1d is time consuming
bool q_dfblcokinterp_linear(DisplaceFieldF3D ***&pppSubDF,
	const V3DLONG szBlock_x, const V3DLONG szBlock_y, const V3DLONG szBlock_z,
	const V3DLONG substart_x, const V3DLONG substart_y, const V3DLONG substart_z,
	DisplaceFieldF3D ***&pppDFBlock);

//bspline interpolate the DF block
//use 3d or 4d pointer instead of 1d, since generating 3d or 4d pointer from 1d is time consuming
bool q_dfblcokinterp_bspline(DisplaceFieldF3D ***&pppSubDF, const Matrix &x_bsplinebasis,
	const V3DLONG sz_gridwnd, const V3DLONG substart_x, const V3DLONG substart_y, const V3DLONG substart_z,
	DisplaceFieldF3D ***&pppDFBlock);

//warp image block based on given DF
//use 3d or 4d pointer instead of 1d, since generating 3d or 4d pointer from 1d is time consuming
//template <class T>
//bool q_imgblockwarp(T ****&p_img_sub_4d, const V3DLONG *sz_img_sub, DisplaceFieldF3D ***&pppDFBlock,
//	const V3DLONG szBlock_x, const V3DLONG szBlock_y, const V3DLONG szBlock_z, const int i_interpmethod_img,
//	const V3DLONG substart_x, const V3DLONG substart_y, const V3DLONG substart_z,
//	T ****&p_img_warp_4d const V3DLONG *sz_img_sub_ori);


//TPS_linear_blockbyblock image warping
//	i_interp_method_df:  0-trilinear, 1-bspline
//	i_interp_method_img: 0-trilinear, 1-nearest neighbor
template <class T>
bool imgwarp_smallmemory(T *p_img_sub, const V3DLONG *sz_img_sub,
	const QList<ImageMarker> &ql_marker_tar, const QList<ImageMarker> &ql_marker_sub,
	V3DLONG szBlock_x, V3DLONG szBlock_y, V3DLONG szBlock_z, int i_interpmethod_df, int i_interpmethod_img,
	T *&p_img_warp, const V3DLONG *sz_img_sub_ori);


bool interpolate_coord_linear(MYFLOAT_JBA * interpolatedVal, Coord3D_JBA *c, V3DLONG numCoord,
	MYFLOAT_JBA *** templateVol3d, V3DLONG tsz0, V3DLONG tsz1, V3DLONG tsz2,
	V3DLONG tlow0, V3DLONG tup0, V3DLONG tlow1, V3DLONG tup1, V3DLONG tlow2, V3DLONG tup2);

Vol3DSimple<DisplaceFieldF3D> * compute_df_tps_subsampled_volume(const vector <Coord3D_JBA> & matchTargetPos, const vector <Coord3D_JBA> & matchSubjectPos, V3DLONG sz0, V3DLONG sz1, V3DLONG sz2,
	V3DLONG gfactor_x, V3DLONG gfactor_y, V3DLONG gfactor_z);


Vol3DSimple<DisplaceFieldF3D> * compute_df_stps_subsampled_volume_4bspline(const vector <Coord3D_JBA> & matchTargetPos, const vector <Coord3D_JBA> & matchSubjectPos, V3DLONG sz0, V3DLONG sz1, V3DLONG sz2,
	V3DLONG gfactor_x, V3DLONG gfactor_y, V3DLONG gfactor_z, int gpu_type);

bool q_nonrigid_ini_bsplinebasis_3D(const long n, Matrix &BxBxB);

Vol3DSimple <MYFLOAT_JBA> * linearinterp_regularmesh_3d(V3DLONG sz0, V3DLONG sz1, V3DLONG sz2, Vol3DSimple <MYFLOAT_JBA> * df_regular_grid);

//Created by CYF
//2022.9.27
//
bool compute_df_stps_subsampled_volume_4bspline_per(const vector <Coord3D_JBA> & matchTargetPos, const vector <Coord3D_JBA> & matchSubjectPos, Matrix &x4x4_d, Matrix &xnx4_c, float * &H_X, float * &H_Y,
	float * &H_Z, int nCpt, Image2DSimple<MYFLOAT_JBA> * &cpt_subject, int gpu_type);

Vol3DSimple<DisplaceFieldF3D> * compute_df_stps_subsampled_volume_4bspline_block(int nCpt, Matrix x4x4_d, Matrix xnx4_c, V3DLONG sz0, V3DLONG sz1, V3DLONG sz2, V3DLONG gfactor_x,
	V3DLONG gfactor_y, V3DLONG gfactor_z, float * H_X, float * H_Y, float * H_Z, long long x_offset, long long y_offset, long long z_offset, int gpu_mode);

bool imgwarp_smallmemory_CYF(const QList<ImageMarker> &ql_marker_tar, const QList<ImageMarker> &ql_marker_sub,
	V3DLONG szBlock_x, V3DLONG szBlock_y, V3DLONG szBlock_z, int i_interpmethod_df, int i_interpmethod_img, string *path_TIF, long long sz_img_ori[4], string TeraWarp_save_base, int GPU_type);

bool imgwarp_smallmemory_image(const QList<ImageMarker> &ql_marker_tar, const QList<ImageMarker> &ql_marker_sub,
                             V3DLONG szBlock_x, V3DLONG szBlock_y, V3DLONG szBlock_z, int i_interpmethod_df, int i_interpmethod_img, unsigned char *p_img_sub, V3DLONG *sz_img_sub, string TeraWarp_save_base, int GPU_type);

//bool STPS_interpolate_CYF(const int GPU_type, Vol3DSimple<DisplaceFieldF3D> *pSubDF, DisplaceFieldF3D ***pppSubDF, Matrix x_bsplinebasis,
//	const V3DLONG sz_gridwnd, DisplaceFieldF3D ***&pppDFBlock, unsigned char ****&p_img_sub_4d, const V3DLONG *sz_img_sub, const V3DLONG szBlock_x,
//	const V3DLONG szBlock_y, const V3DLONG szBlock_z, const int i_interpmethod_img, unsigned char ****&p_img_warp_4d, const long long *sz_img_sub_ori,
//	const unsigned char *p_img_sub, unsigned char *p_img_warp, long long x_offset, long long y_offset, long long z_offset, const long long *sz_img_ori,
//	const long long *sz_img_read_block);

bool q_imgblockwarp_CYF(unsigned char ****&p_img_sub_4d, const V3DLONG *sz_img_sub, DisplaceFieldF3D ***&pppDFBlock,
	const V3DLONG szBlock_x, const V3DLONG szBlock_y, const V3DLONG szBlock_z, const int i_interpmethod_img,
	const V3DLONG substart_x, const V3DLONG substart_y, const V3DLONG substart_z, unsigned char ****&p_img_warp_4d,
	const V3DLONG *sz_img_sub_read, long long x_offset, long long y_offset, long long z_offset, const long long *sz_img_ori,
	const long long *sz_img_read_block);


bool STPS_interpolate_sort_CYF(const int GPU_type, Vol3DSimple<DisplaceFieldF3D> *pSubDF, DisplaceFieldF3D ***pppSubDF, Matrix x_bsplinebasis,
	const V3DLONG sz_gridwnd, DisplaceFieldF3D ***&pppDFBlock, const V3DLONG *sz_img_sub, const V3DLONG szBlock_x, const V3DLONG szBlock_y,
	const V3DLONG szBlock_z, const int i_interpmethod_img, long long x_offset, long long y_offset, long long z_offset, const long long *sz_img_ori,
	float * &sort_x, float * &sort_y, float * &sort_z);

#endif