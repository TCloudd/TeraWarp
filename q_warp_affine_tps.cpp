// q_warp_affine_tps.cpp
// warp pointset and image based on given matched pairs
// by Lei Qu
// 2010-03-22

#include <QtGui>
#include "basic_memory.cpp"
#include "q_warp_affine_tps.h"
#include "stackutil.h"
#include "Bigwarp.h"
#include "basic_memory.cpp"//note: should not include .h file, since they are template functions
#include "basic_surf_objs.h"
#include <string>
#include <sstream>
#include "dir_make.h"
#include "RawFmtMngr.h"
#include <vector>
#include <algorithm>
#define EPS 0.0001
#define GPU_LIMIT 300
extern "C" bool gpu_interpolation_new(const int mode, const long long gsz2, const long long gsz1, const long long gsz0, const Matrix &x4x4_affinematrix, unsigned char ****&p_img_sub_4d,
	unsigned char ****&p_img_sub2tar_4d, const long long *sz_img_sub, const long long gsA2, const long long gsA1, const long long gsA0, const long long start_block_x, const long long start_block_y,
	const long long start_block_z, const long long x_read_offset, const long long y_read_offset, const long long z_read_offset, const long long gs_ori2, const long long gs_ori1, const long long gs_ori0,
	const unsigned char *p_img_sub, unsigned char *p_img_affine);
extern "C" bool gpu_interpolation_affine(const long long gsz2, const long long gsz1, const long long gsz0, const Matrix &x4x4_affinematrix, unsigned char ****&p_img_sub_4d, unsigned char ****&p_img_sub2tar_4d, const long long *sz_img_sub,
                                         const long long gsA2, const long long gsA1, const long long gsA0);

//extern "C" int gpu_A_i_new(int ncpt, const Matrix &A, Matrix &A_i);
//extern "C" Matrix matrixMultiply(const int m, const int n, const int k, Matrix &A, Matrix &B);

//affine points warp


bool q_ptswarp_affine(const vector<Coord3D_PCM> &vec_ctlpt_tar, const vector<Coord3D_PCM>  &vec_ctlpt_sub,
	vector<Coord3D_PCM> &vec_ctlpt_subtar_affine)
{
	//check parameters
	if (vec_ctlpt_tar.size() == 0 || vec_ctlpt_sub.size() == 0 || vec_ctlpt_tar.size() != vec_ctlpt_sub.size())
	{
		printf("ERROR: target or subject control points is invalid!\n");
		return false;
	}
	if (vec_ctlpt_subtar_affine.size() != 0) vec_ctlpt_subtar_affine.clear();

	//------------------------------------------------------------------------------------------------------------------------------------
	//estimate the affine matrix
	Matrix x4x4_affinematrix;
	if (!q_affine_compute_affinmatrix_3D(vec_ctlpt_sub, vec_ctlpt_tar, x4x4_affinematrix))	//B=T*A
	{
		printf("ERROR: q_affine_compute_affinmatrix_2D() return false.\n");
		return false;
	}
	//print affine matrix
	for (long long row = 1; row <= x4x4_affinematrix.nrows(); row++)
	{
		printf("\t");
		for (V3DLONG col = 1; col <= x4x4_affinematrix.ncols(); col++)
			printf("%.3f\t", x4x4_affinematrix(row, col));
		printf("\n");
	}

	//------------------------------------------------------------------------------------------------------------------------------------
	//affine points warping
	Matrix x_pt_sub2tar_homo(4, 1), x_pt_sub_homo(4, 1);
	for (unsigned int i = 0; i < vec_ctlpt_sub.size(); i++)
	{
		//compute the inverse affine projected coordinate in subject image
		x_pt_sub2tar_homo(1, 1) = vec_ctlpt_sub[i].x;
		x_pt_sub2tar_homo(2, 1) = vec_ctlpt_sub[i].y;
		x_pt_sub2tar_homo(3, 1) = vec_ctlpt_sub[i].z;
		x_pt_sub2tar_homo(4, 1) = 1.0;
		x_pt_sub_homo = x4x4_affinematrix*x_pt_sub2tar_homo;

		Coord3D_PCM tmp;
		tmp.x = x_pt_sub_homo(1, 1);
		tmp.y = x_pt_sub_homo(2, 1);
		tmp.z = x_pt_sub_homo(3, 1);
		vec_ctlpt_subtar_affine.push_back(tmp);
	}

	return true;
}

//affine image warp
bool q_imagewarp_affine(const vector<Coord3D_PCM> &vec_ctlpt_tar,const vector<Coord3D_PCM>  &vec_ctlpt_sub,
	const unsigned char *p_img_sub, const long long *sz_img_sub, const long long *sz_img_affine,
	unsigned char *&p_img_affine, int type_gpu, int start_block_x, int start_block_y, int start_block_z)
{
	//check parameters
	if (vec_ctlpt_tar.size() == 0 || vec_ctlpt_sub.size() == 0 || vec_ctlpt_tar.size() != vec_ctlpt_sub.size())
	{
		printf("ERROR: target or subject control points is invalid!\n");
		return false;
	}
	if (p_img_sub == 0 || sz_img_sub == 0)
	{
		printf("ERROR: p_img_sub or sz_img_sub is invalid.\n");
		return false;
	}
	if (p_img_affine)
	{
		printf("WARNNING: output image pointer is not null, original memory it point to will lost!\n");
		p_img_affine = 0;
	}

	//------------------------------------------------------------------------------------------------------------------------------------
	//assign output/warp image size
	long long sz_img_output[4] = { 0 };
	if (sz_img_output[0] == 0)		sz_img_output[0] = sz_img_affine[0];
	if (sz_img_output[1] == 0)		sz_img_output[1] = sz_img_affine[1];
	if (sz_img_output[2] == 0)		sz_img_output[2] = sz_img_affine[2];
	sz_img_output[3] = sz_img_sub[3];

	//allocate memory
	p_img_affine = new unsigned char[sz_img_output[0] * sz_img_output[1] * sz_img_output[2] * sz_img_output[3]]();
	if (!p_img_affine)
	{
		printf("ERROR: Fail to allocate memory for p_img_sub2tar.\n");
		return false;
	}
	unsigned char ****p_img_sub_4d = 0, ****p_img_sub2tar_4d = 0;
	if (!new4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3], p_img_sub) ||
		!new4dpointer(p_img_sub2tar_4d, sz_img_output[0], sz_img_output[1], sz_img_output[2], sz_img_output[3], p_img_affine))
	{
		printf("ERROR: Fail to allocate memory for the 4d pointer of image.\n");
		if (p_img_affine) 		{ delete[]p_img_affine;		p_img_affine = 0; }
		if (p_img_sub_4d) 		{ delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
		if (p_img_sub2tar_4d) 	{ delete4dpointer(p_img_sub2tar_4d, sz_img_output[0], sz_img_output[1], sz_img_output[2], sz_img_output[3]); }
		return false;
	}

	//------------------------------------------------------------------------------------------------------------------------------------
	//estimate the affine matrix
	Matrix x4x4_affinematrix;
	//clock_t affine_mat;
	//affine_mat = clock();
	if (!q_affine_compute_affinmatrix_3D(vec_ctlpt_tar, vec_ctlpt_sub, x4x4_affinematrix))	//B=T*A
	{
		printf("ERROR: q_affine_compute_affinmatrix_2D() return false.\n");
		if (p_img_sub_4d) 		{ delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
		if (p_img_sub2tar_4d) 	{ delete4dpointer(p_img_sub2tar_4d, sz_img_output[0], sz_img_output[1], sz_img_output[2], sz_img_output[3]); }
		if (p_img_affine) 		{ delete[]p_img_affine;		p_img_affine = 0; }
		return false;
	}
	//------------------------------------------------------------------------------------------------------------------------------------
	
	//print affine matrix
/*	for (long long row = 1; row <= x4x4_affinematrix.nrows(); row++)
	{
		printf("\t");
		for (long long col = 1; col <= x4x4_affinematrix.ncols(); col++)
			printf("%.3f\t", x4x4_affinematrix(row, col));
		printf("\n");
	}*/
	/*	SVD���Դ���
	Matrix A(3,3), U(3,3), V(3,3), A_I(3,3);
	A(1, 1) = 1.0;
	A(1, 2) = 3.0;
	A(1, 3) = 4.0;
	A(2, 1) = 2.0;
	A(2, 2) = 5.0;
	A(2, 3) = 2.0;
	A(3, 1) = 3.0;
	A(3, 2) = 4.0;
	A(3, 3) = 5.0;
	//printf("%.3f\t", A(2, 1));
	for (long long row = 1; row <= A.nrows(); row++)
	{
		printf("\t");
		for (long long col = 1; col <= A.ncols(); col++)
			printf("%.3f\t", A(row, col));
		printf("\n");
	}

	DiagonalMatrix D;
	SVD(A, D, U, V);
	for (long long row = 1; row <= V.nrows(); row++)
	{
		printf("\t");
		for (long long col = 1; col <= V.ncols(); col++)
			printf("%.3f\t", V(row, col));
		printf("\n");
	}
	Matrix Dinv = D.i();
	A_I = V*Dinv*U.t();
	for (long long row = 1; row <= A_I.nrows(); row++)
	{
		printf("\t");
		for (long long col = 1; col <= A_I.ncols(); col++)
			printf("%.3f\t", A_I(row, col));
		printf("\n");
	}
*/

	//------------------------------------------------------------------------------------------------------------------------------------
	//affine image warping

	/*for (long long a = 0; a < sz_img_sub[2]; a++)
	{
		for (long long b = 0; b < sz_img_sub[1]; b++)
		{
			for (long long c = 0; c < sz_img_sub[0]; c++)
			{

				if (p_img_sub_4d[0][a][b][c] != 0)printf("H: %d\n", p_img_sub_4d[0][a][b][c]);
			}
		}
	}*/
	if (type_gpu == 0)
	{
#pragma omp parallel for
		for (long long x = 0; x < sz_img_output[0]; x++)
		{
			printf("affine: [%d/%d]\n", sz_img_output[0], x);
			for (long long y = 0; y < sz_img_output[1]; y++)
				for (long long z = 0; z < sz_img_output[2]; z++)
				{
					Matrix x_pt_sub2tar_homo(4, 1), x_pt_sub_homo(4, 1);
					//compute the inverse affine projected coordinate in subject image���������ͶӰ����������ͼ��
					x_pt_sub2tar_homo(1, 1) = x;
					x_pt_sub2tar_homo(2, 1) = y;
					x_pt_sub2tar_homo(3, 1) = z;
					x_pt_sub2tar_homo(4, 1) = 1.0;
					x_pt_sub_homo = x4x4_affinematrix*x_pt_sub2tar_homo;

					//------------------------------------------------------------------
					//linear interpolate���Բ���
					//coordinate in subject image����ͼ���е�����
					double cur_pos[3];//x,y,z
					cur_pos[0] = x_pt_sub_homo(1, 1) - start_block_x;
					cur_pos[1] = x_pt_sub_homo(2, 1) - start_block_y;
					cur_pos[2] = x_pt_sub_homo(3, 1) - start_block_z;

					//if interpolate pixel is out of subject image region, set to -inf�����ֵ���س�������ͼ����������Ϊ-inf
					if (cur_pos[0]<0 || cur_pos[0]>sz_img_sub[0] - 1 ||
						cur_pos[1]<0 || cur_pos[1]>sz_img_sub[1] - 1 ||
						cur_pos[2]<0 || cur_pos[2]>sz_img_sub[2] - 1)
					{
						p_img_sub2tar_4d[0][z][y][x] = 0.0;
						continue;
					}

					//find 8 neighbor pixels boundary�ҳ�8���������صı߽�
					long long x_s, x_b, y_s, y_b, z_s, z_b;
					x_s = floor(cur_pos[0]);		x_b = ceil(cur_pos[0]);
					y_s = floor(cur_pos[1]);		y_b = ceil(cur_pos[1]);
					z_s = floor(cur_pos[2]);		z_b = ceil(cur_pos[2]);

					//compute weight for left and right, top and bottom -- 4 neighbor pixel's weight in a slice�������ҡ��ϡ��µ�Ȩ�ء���һ����Ƭ��4���ھ����ص�Ȩ��
					double l_w, r_w, t_w, b_w;
					l_w = 1.0 - (cur_pos[0] - x_s);	r_w = 1.0 - l_w;//x_d
					t_w = 1.0 - (cur_pos[1] - y_s);	b_w = 1.0 - t_w;//y_d
					//compute weight for higher slice and lower slice����ϸ���Ƭ�ͽϵ���Ƭ��Ȩ��
					double u_w, d_w;
					u_w = 1.0 - (cur_pos[2] - z_s);	d_w = 1.0 - u_w;//z_d

					//linear interpolate each channel
					for (long long c = 0; c < sz_img_output[3]; c++)
					{
						//linear interpolate in higher slice [t_w*(l_w*lt+r_w*rt)+b_w*(l_w*lb+r_w*rb)]
						double higher_slice;
						higher_slice = t_w*(l_w*p_img_sub_4d[c][z_s][y_s][x_s] + r_w*p_img_sub_4d[c][z_s][y_s][x_b]) +
							b_w*(l_w*p_img_sub_4d[c][z_s][y_b][x_s] + r_w*p_img_sub_4d[c][z_s][y_b][x_b]);
						//linear interpolate in lower slice [t_w*(l_w*lt+r_w*rt)+b_w*(l_w*lb+r_w*rb)]
						double lower_slice;
						lower_slice = t_w*(l_w*p_img_sub_4d[c][z_b][y_s][x_s] + r_w*p_img_sub_4d[c][z_b][y_s][x_b]) +
							b_w*(l_w*p_img_sub_4d[c][z_b][y_b][x_s] + r_w*p_img_sub_4d[c][z_b][y_b][x_b]);
						//linear interpolate the current position [u_w*higher_slice+d_w*lower_slice]
						p_img_sub2tar_4d[c][z][y][x] = u_w*higher_slice + d_w*lower_slice;
					}
				}
		}
	}
	else if (type_gpu==1)
	{
		//GPU����
		long long gsA2 = sz_img_sub[2];
		long long gsA1 = sz_img_sub[1];
		long long gsA0 = sz_img_sub[0];

		long long gsz2 = sz_img_output[2];
		long long gsz1 = sz_img_output[1];
		long long gsz0 = sz_img_output[0];
		printf("gsz2:%d,gsz1:%d,gsz0:%d\n", gsz2, gsz1, gsz0);
		clock_t stps_interpolation;
		stps_interpolation = clock();

		//GPU���ٺ���
		//gpu_interpolation_noblock(gsz2, gsz1, gsz0, x4x4_affinematrix, p_img_sub_4d, p_img_sub2tar_4d, sz_img_output, gsA2, gsA1, gsA0);
		printf("\t>>xnxn_K time consume %.2f s\n", (float)(clock() - stps_interpolation) / CLOCKS_PER_SEC);

		
	}
/*#pragma omp parallel for
	for (long long x = 0; x < sz_img_output[0]; x++)
	{
		printf("affine: [%d/%d]\n", sz_img_output[0], x);
		for (long long y = 0; y < sz_img_output[1]; y++)
			for (long long z = 0; z < sz_img_output[2]; z++)
			{
				Matrix x_pt_sub2tar_homo(4, 1), x_pt_sub_homo(4, 1);
				//compute the inverse affine projected coordinate in subject image���������ͶӰ����������ͼ��
				x_pt_sub2tar_homo(1, 1) = x;
				x_pt_sub2tar_homo(2, 1) = y;
				x_pt_sub2tar_homo(3, 1) = z;
				x_pt_sub2tar_homo(4, 1) = 1.0;
				x_pt_sub_homo = x4x4_affinematrix*x_pt_sub2tar_homo;

				//------------------------------------------------------------------
				//linear interpolate���Բ���
				//coordinate in subject image����ͼ���е�����
				double cur_pos[3];//x,y,z
				cur_pos[0] = x_pt_sub_homo(1, 1);
				cur_pos[1] = x_pt_sub_homo(2, 1);
				cur_pos[2] = x_pt_sub_homo(3, 1);

				//if interpolate pixel is out of subject image region, set to -inf�����ֵ���س�������ͼ����������Ϊ-inf
				if (cur_pos[0]<0 || cur_pos[0]>sz_img_sub[0] - 1 ||
					cur_pos[1]<0 || cur_pos[1]>sz_img_sub[1] - 1 ||
					cur_pos[2]<0 || cur_pos[2]>sz_img_sub[2] - 1)
				{
					p_img_sub2tar_4d[0][z][y][x] = 0.0;
					continue;
				}

				//find 8 neighbor pixels boundary�ҳ�8���������صı߽�
				long long x_s, x_b, y_s, y_b, z_s, z_b;
				x_s = floor(cur_pos[0]);		x_b = ceil(cur_pos[0]);
				y_s = floor(cur_pos[1]);		y_b = ceil(cur_pos[1]);
				z_s = floor(cur_pos[2]);		z_b = ceil(cur_pos[2]);

				//compute weight for left and right, top and bottom -- 4 neighbor pixel's weight in a slice�������ҡ��ϡ��µ�Ȩ�ء���һ����Ƭ��4���ھ����ص�Ȩ��
				double l_w, r_w, t_w, b_w;
				l_w = 1.0 - (cur_pos[0] - x_s);	r_w = 1.0 - l_w;
				t_w = 1.0 - (cur_pos[1] - y_s);	b_w = 1.0 - t_w;
				//compute weight for higher slice and lower slice����ϸ���Ƭ�ͽϵ���Ƭ��Ȩ��
				double u_w, d_w;
				u_w = 1.0 - (cur_pos[2] - z_s);	d_w = 1.0 - u_w;

				//linear interpolate each channel
				for (long long c = 0; c < sz_img_output[3]; c++)
				{
					//linear interpolate in higher slice [t_w*(l_w*lt+r_w*rt)+b_w*(l_w*lb+r_w*rb)]
					double higher_slice;
					higher_slice = t_w*(l_w*p_img_sub_4d[c][z_s][y_s][x_s] + r_w*p_img_sub_4d[c][z_s][y_s][x_b]) +
						b_w*(l_w*p_img_sub_4d[c][z_s][y_b][x_s] + r_w*p_img_sub_4d[c][z_s][y_b][x_b]);
					//linear interpolate in lower slice [t_w*(l_w*lt+r_w*rt)+b_w*(l_w*lb+r_w*rb)]
					double lower_slice;
					lower_slice = t_w*(l_w*p_img_sub_4d[c][z_b][y_s][x_s] + r_w*p_img_sub_4d[c][z_b][y_s][x_b]) +
						b_w*(l_w*p_img_sub_4d[c][z_b][y_b][x_s] + r_w*p_img_sub_4d[c][z_b][y_b][x_b]);
					//linear interpolate the current position [u_w*higher_slice+d_w*lower_slice]
					p_img_sub2tar_4d[c][z][y][x] = u_w*higher_slice + d_w*lower_slice;
				}
			}
	}*/
	
	//------------------------------------------------------------------------------------------------------------------------------------
	printf("6. free memory. \n");
	if (p_img_sub_4d) 		{ delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
	if (p_img_sub2tar_4d) 	{ delete4dpointer(p_img_sub2tar_4d, sz_img_output[0], sz_img_output[1], sz_img_output[2], sz_img_output[3]); }

	return true;
}

//Read matched-pair index file
//	output vec2D_sub2tar_matchind is a 2D (n*2) vector
//		vec2D_sub2tar_matchind[i][0]: sub index of i-th matched pair
//		vec2D_sub2tar_matchind[i][1]: tar index of i-th matched pair
bool q_readMatchInd_file(const QString qs_filename,vector< vector<long> > &vec2D_sub2tar_matchind)
{
	//check parameters
	if(qs_filename.isEmpty())
	{
		fprintf(stderr,"ERROR: Invalid input file name! \n");
		return false;
	}
	if(vec2D_sub2tar_matchind.size()!=0)
	{
		vec2D_sub2tar_matchind.clear();
	}

	vector<long> vec1D_sub2tar_matchind(2,-1);		//0: sub_index; 1: tar_index

	QFile qf(qs_filename);
	if(!qf.open(QIODevice::ReadOnly | QIODevice::Text))
	{
		printf("ERROR: open file [%s] fail.\n",qPrintable(qs_filename));
		return false;
	}

	long k=0;
    while(!qf.atEnd())
    {
		char curline[2000];
        qf.readLine(curline,sizeof(curline));
		k++;
		{
			if(curline[0]=='#' || curline[0]=='\n' || curline[0]=='\0') continue;

			QStringList qsl=QString(curline).trimmed().split(":");
			int qsl_count=qsl.size();
			if(qsl_count!=2)
			{
				printf("WARNING: invalid format found in line %d.\n",k);
				continue;
			}

			if(qsl[0].toLong()>=0 && qsl[1].toLong()>=0)
			{
				vec1D_sub2tar_matchind[0]=qsl[0].toLong();	//sub
				vec1D_sub2tar_matchind[1]=qsl[1].toLong();	//tar
				vec2D_sub2tar_matchind.push_back(vec1D_sub2tar_matchind);
			}
			else
				printf("WARNING: invalid matched pair index (index<0) found in line %d.\n",k);
		}
    }

    return true;
}


//centrilize and scale the point set
//	xn = T*x;
//	x: every column represent a point [2/3*N]
bool q_normalize_points_2D(const vector<Coord3D_PCM> vec_input,vector<Coord3D_PCM> &vec_output,Matrix &x3x3_normalize)
{
	//check parameters
	if(vec_input.size()<=0)
	{
		fprintf(stderr,"ERROR: Input array is null! \n");
		return false;
	}
	if(!vec_output.empty())
		vec_output.clear();
	vec_output=vec_input;
	if(x3x3_normalize.nrows()!=3 || x3x3_normalize.ncols()!=3)
	{
		x3x3_normalize.ReSize(3,3);
	}

	//compute the centriod of input point set
	Coord3D_PCM cord_centroid;
	int n_point=vec_input.size();
	for(int i=0;i<n_point;i++)
	{
		cord_centroid.x+=vec_input[i].x;
		cord_centroid.y+=vec_input[i].y;
	}
	cord_centroid.x/=n_point;
	cord_centroid.y/=n_point;
	//center the point set
	for(int i=0;i<n_point;i++)
	{
		vec_output[i].x-=cord_centroid.x;
		vec_output[i].y-=cord_centroid.y;
	}

	//compute the average distance of every point to the origin
	double d_point2o=0,d_point2o_avg=0;
	for(int i=0;i<n_point;i++)
	{
		d_point2o=sqrt(vec_output[i].x*vec_output[i].x+vec_output[i].y*vec_output[i].y);
		d_point2o_avg+=d_point2o;
	}
	d_point2o_avg/=n_point;
	//compute the scale factor
	double d_scale_factor=1.0/d_point2o_avg;
	//scale the point set
	for(int i=0;i<n_point;i++)
	{
		vec_output[i].x*=d_scale_factor;
		vec_output[i].y*=d_scale_factor;
	}

	//compute the transformation matrix
	// 1 row
	x3x3_normalize(1,1)=d_scale_factor;
	x3x3_normalize(1,2)=0;
	x3x3_normalize(1,3)=-d_scale_factor*cord_centroid.x;
	// 2 row
	x3x3_normalize(2,1)=0;
	x3x3_normalize(2,2)=d_scale_factor;
	x3x3_normalize(2,3)=-d_scale_factor*cord_centroid.y;
	// 3 row
	x3x3_normalize(3,1)=0;
	x3x3_normalize(3,2)=0;
	x3x3_normalize(3,3)=1;

	return true;
}

bool q_normalize_points_3D(const vector<Coord3D_PCM> vec_input,vector<Coord3D_PCM> &vec_output,Matrix &x4x4_normalize)
{
	//check parameters
	if(vec_input.size()<=0)
	{
		fprintf(stderr,"ERROR: Input array is null! \n");
		return false;
	}
	if(!vec_output.empty())
		vec_output.clear();
	vec_output=vec_input;
	if(x4x4_normalize.nrows()!=4 || x4x4_normalize.ncols()!=4)
	{
		x4x4_normalize.ReSize(4,4);
	}

	//compute the centriod of input point set ��������㼯������
	Coord3D_PCM cord_centroid;
	int n_point=vec_input.size();
	for(int i=0;i<n_point;i++)
	{
		cord_centroid.x+=vec_input[i].x;
		cord_centroid.y+=vec_input[i].y;
		cord_centroid.z+=vec_input[i].z;
	}
	cord_centroid.x/=n_point;
	cord_centroid.y/=n_point;
	cord_centroid.z/=n_point;
	//center the point set ʹ�㼯����
	for(int i=0;i<n_point;i++)
	{
		vec_output[i].x-=cord_centroid.x;
		vec_output[i].y-=cord_centroid.y;
		vec_output[i].z-=cord_centroid.z;
	}

	//compute the average distance of every point to the origin ����ÿ���㵽ԭ���ƽ������  
	double d_point2o=0,d_point2o_avg=0;
	for(int i=0;i<n_point;i++)
	{
		d_point2o=sqrt(vec_output[i].x*vec_output[i].x+vec_output[i].y*vec_output[i].y+vec_output[i].z*vec_output[i].z);
		d_point2o_avg+=d_point2o;
	}
	d_point2o_avg/=n_point;
	//compute the scale factor �����������
	double d_scale_factor=1.0/d_point2o_avg;
	//scale the point set �Ե㼯��������
	for(int i=0;i<n_point;i++)
	{
		vec_output[i].x*=d_scale_factor;
		vec_output[i].y*=d_scale_factor;
		vec_output[i].z*=d_scale_factor;
	}

	//compute the transformation matrix ����任����
	// 1 row
	x4x4_normalize(1,1)=d_scale_factor;
	x4x4_normalize(1,2)=0;
	x4x4_normalize(1,3)=0;
	x4x4_normalize(1,4)=-d_scale_factor*cord_centroid.x;
	// 2 row
	x4x4_normalize(2,1)=0;
	x4x4_normalize(2,2)=d_scale_factor;
	x4x4_normalize(2,3)=0;
	x4x4_normalize(2,4)=-d_scale_factor*cord_centroid.y;
	// 3 row
	x4x4_normalize(3,1)=0;
	x4x4_normalize(3,2)=0;
	x4x4_normalize(3,3)=d_scale_factor;
	x4x4_normalize(3,4)=-d_scale_factor*cord_centroid.z;
	// 4 row
	x4x4_normalize(4,1)=0;
	x4x4_normalize(4,2)=0;
	x4x4_normalize(4,3)=0;
	x4x4_normalize(4,4)=1;

	return true;
}

//compute the affine matraix
//	B=T*A
bool q_affine_compute_affinmatrix_2D(const vector<Coord3D_PCM> &vec_A,const vector<Coord3D_PCM> &vec_B,Matrix &x3x3_affinematrix)
{
	//check parameters
	if(vec_A.size()<3 || vec_A.size()!=vec_B.size())
	{
		fprintf(stderr,"ERROR: Invalid input parameters! \n");
		return false;
	}
	if(x3x3_affinematrix.nrows()!=3 || x3x3_affinematrix.ncols()!=3)
	{
		x3x3_affinematrix.ReSize(3,3);
	}

	//normalize point set
	vector<Coord3D_PCM> vec_A_norm,vec_B_norm;
	Matrix x3x3_normalize_A(4,4),x3x3_normalize_B(4,4);
	vec_A_norm=vec_A;	vec_B_norm=vec_B;
	q_normalize_points_2D(vec_A,vec_A_norm,x3x3_normalize_A);
	q_normalize_points_2D(vec_B,vec_B_norm,x3x3_normalize_B);

	//	 fill matrix A
	//
	//	  | h1, h2, h3 |    |x1| |x2|
	//	  | h4, h5, h6 | *  |y1|=|y2| <=>
	//	  | 0 ,  0,  1 |    |1 | |1 |
	//
	//	  |x1, y1, 1,  0,  0,  0, -x2 |
	//	  | 0,  0,  0, x1, y1, 1, -y2 | * |h1,h2,h3,h4,h5,h6,h7|=0
	//
	int n_point=vec_A.size();
	Matrix A(2*n_point,7);
	int row=1;
	for(int i=0;i<n_point;i++)
	{
		A(row,1)=vec_A_norm[i].x;	A(row,2)=vec_A_norm[i].y;	A(row,3)=1.0;
		A(row,4)=0.0;				A(row,5)=0.0;				A(row,6)=0.0;
		A(row,7)=-vec_B_norm[i].x;

		A(row+1,1)=0.0;				A(row+1,2)=0.0;				A(row+1,3)=0.0;
		A(row+1,4)=vec_A_norm[i].x;	A(row+1,5)=vec_A_norm[i].y;	A(row+1,6)=1.0;
		A(row+1,7)=-vec_B_norm[i].y;

		row+=2;
	}

	//compute T  --> bug? SVD in newmat need row>=col?
	DiagonalMatrix D;
	Matrix U,V;
	SVD(A,D,U,V);	//A = U * D * V.t()

	Matrix h=V.column(7);	//A*h=0
	if(D(6,6)==0)			//degenerate case
	{
		x3x3_affinematrix=0.0;	//check with A.is_zero()
		printf("Degenerate singular values in SVD! \n");
		//		return false;
	}

	//de-homo
	for(int i=1;i<=7;i++)
	{
		h(i,1) /= h(7,1);
	}

	//reshape h:7*1 to 3*3 matrix
	x3x3_affinematrix(1,1)=h(1,1);	x3x3_affinematrix(1,2)=h(2,1);	x3x3_affinematrix(1,3)=h(3,1);
	x3x3_affinematrix(2,1)=h(4,1);	x3x3_affinematrix(2,2)=h(5,1);	x3x3_affinematrix(2,3)=h(6,1);
	x3x3_affinematrix(3,1)=0.0;		x3x3_affinematrix(3,2)=0.0;		x3x3_affinematrix(3,3)=1.0;

	//denormalize
	x3x3_affinematrix=x3x3_normalize_B.i()*x3x3_affinematrix*x3x3_normalize_A;

	return true;
}

bool q_affine_compute_affinmatrix_3D(const vector<Coord3D_PCM> &vec_A,const vector<Coord3D_PCM> &vec_B,Matrix &x4x4_affinematrix)
{
	//check parameters
	if(vec_A.size()<4 || vec_A.size()!=vec_B.size())
	{
		fprintf(stderr,"ERROR: Invalid input parameters! \n");
		return false;
	}
	if(x4x4_affinematrix.nrows()!=4 || x4x4_affinematrix.ncols()!=4)
	{
		x4x4_affinematrix.ReSize(4,4);
	}

	//normalize point set
	vector<Coord3D_PCM> vec_A_norm,vec_B_norm;
	Matrix x4x4_normalize_A(4,4),x4x4_normalize_B(4,4);
	vec_A_norm=vec_A;	vec_B_norm=vec_B;
	q_normalize_points_3D(vec_A,vec_A_norm,x4x4_normalize_A);
	q_normalize_points_3D(vec_B,vec_B_norm,x4x4_normalize_B);

	//fill matrix A
	//
	//	  | h1, h2, h3, h4 |    |x1| |x2|
	//	  | h5, h6, h7, h8 | *  |y1|=|y2| <=>
	//	  | h9, h10,h11,h12|    |z1| |z2|
	//	  | 0 ,  0,  0,  1 |    |1 | |1 |
	//
	//	  |x1, y1, z1, 1,  0,  0,  0,  0,  0,  0,  0,  0, -x2 |
	//	  | 0,  0,  0, 0, x1, y1, z1,  1,  0,  0,  0,  0, -y2 | * |h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13|=0
	//	  | 0,  0,  0, 0, 0, 0, 0, 0,  0, x1, y1, z1,  1, -z2 |
	int n_point=vec_A.size();
	Matrix A(3*n_point,13);
	int row=1;
	for(int i=0;i<n_point;i++)
	{
		A(row,1)=vec_A_norm[i].x;	A(row,2)=vec_A_norm[i].y;	A(row,3)=vec_A_norm[i].z;	A(row,4)=1.0;
		A(row,5)=0.0;				A(row,6)=0.0;				A(row,7)=0.0;				A(row,8)=0.0;
		A(row,9)=0.0;				A(row,10)=0.0;				A(row,11)=0.0;				A(row,12)=0.0;
		A(row,13)=-vec_B_norm[i].x;

		A(row+1,1)=0.0;				A(row+1,2)=0.0;				A(row+1,3)=0.0;				A(row+1,4)=0.0;
		A(row+1,5)=vec_A_norm[i].x;	A(row+1,6)=vec_A_norm[i].y;	A(row+1,7)=vec_A_norm[i].z;	A(row+1,8)=1.0;
		A(row+1,9)=0.0;				A(row+1,10)=0.0;			A(row+1,11)=0.0;			A(row+1,12)=0.0;
		A(row+1,13)=-vec_B_norm[i].y;

		A(row+2,1)=0.0;				A(row+2,2)=0.0;				A(row+2,3)=0.0;				A(row+2,4)=0.0;
		A(row+2,5)=0.0;				A(row+2,6)=0.0;				A(row+2,7)=0.0;				A(row+2,8)=0.0;
		A(row+2,9)=vec_A_norm[i].x;	A(row+2,10)=vec_A_norm[i].y;A(row+2,11)=vec_A_norm[i].z;A(row+2,12)=1.0;
		A(row+2,13)=-vec_B_norm[i].z;

		row+=3;
	}

	//compute T  --> bug? SVD in newmat need row>=col?
	DiagonalMatrix D;
	Matrix U,V;
	SVD(A,D,U,V);	//A = U * D * V.t()

	Matrix h=V.column(13);	//A*h=0
	if(D(12,12)==0)			//degenerate case
	{
		x4x4_affinematrix=0.0;	//check with A.is_zero()
		printf("Degenerate singular values in SVD! \n");
		//		return false;
	}

	//de-homo
	for(int i=1;i<=13;i++)
	{
		h(i,1) /= h(13,1);
	}

	//reshape h:13*1 to 4*4 matrix
	x4x4_affinematrix(1,1)=h(1,1);	x4x4_affinematrix(1,2)=h(2,1);	x4x4_affinematrix(1,3)=h(3,1);	x4x4_affinematrix(1,4)=h(4,1);
	x4x4_affinematrix(2,1)=h(5,1);	x4x4_affinematrix(2,2)=h(6,1);	x4x4_affinematrix(2,3)=h(7,1);	x4x4_affinematrix(2,4)=h(8,1);
	x4x4_affinematrix(3,1)=h(9,1);	x4x4_affinematrix(3,2)=h(10,1);	x4x4_affinematrix(3,3)=h(11,1);	x4x4_affinematrix(3,4)=h(12,1);
	x4x4_affinematrix(4,1)=0.0;		x4x4_affinematrix(4,2)=0.0;		x4x4_affinematrix(4,3)=0.0;		x4x4_affinematrix(4,4)=1.0;

	//denormalize
	x4x4_affinematrix=x4x4_normalize_B.i()*x4x4_affinematrix*x4x4_normalize_A;

	return true;
}


//affine image warp ����任����@CYF2022.6.15
bool q_imagewarp_affine_chazhi_CYF(int type_gpu, const long long *sz_img_output,const Matrix &x4x4_affinematrix, long long x_offset, long long y_offset,long long z_offset, const long long *sz_img_sub,
	unsigned char ****&p_img_sub_4d, unsigned char ****&p_img_sub2tar_4d, const long long *sz_img_read_block, const long long *sz_img_read_size,const long long *sz_img_sort,const int cal_mode,
	const unsigned char *p_img_sub, unsigned char *p_img_affine)
{
	if (type_gpu == 0)
	{
     #pragma omp parallel for
		for (long long x = 0; x < sz_img_output[0]; x++)
		{
			//printf("affine: [%d/%d]\n", sz_img_output[0], x);
			for (long long y = 0; y < sz_img_output[1]; y++)
				for (long long z = 0; z < sz_img_output[2]; z++)
				{
					Matrix x_pt_sub2tar_homo(4, 1), x_pt_sub_homo(4, 1);
					//compute the inverse affine projected coordinate in subject image���������ͶӰ����������ͼ��
					x_pt_sub2tar_homo(1, 1) = x + x_offset;
					x_pt_sub2tar_homo(2, 1) = y + y_offset;
					x_pt_sub2tar_homo(3, 1) = z + z_offset;
					x_pt_sub2tar_homo(4, 1) = 1.0;
					x_pt_sub_homo = x4x4_affinematrix*x_pt_sub2tar_homo;
					//printf("mark31;(%d,%d,%d)\n",x,y,z);
					//------------------------------------------------------------------
					//linear interpolate���Բ���
					//coordinate in subject image����ͼ���е�����
					double cur_pos[3];//x,y,z
					double cur_pos_1[3];//x,y,z
					cur_pos[0] = x_pt_sub_homo(1, 1) - sz_img_read_block[0];
					cur_pos[1] = x_pt_sub_homo(2, 1) - sz_img_read_block[1];
					cur_pos[2] = x_pt_sub_homo(3, 1) - sz_img_read_block[2];

					cur_pos_1[0] = x_pt_sub_homo(1, 1) ;
					cur_pos_1[1] = x_pt_sub_homo(2, 1);
					cur_pos_1[2] = x_pt_sub_homo(3, 1) ;
					//printf("marked1");
					//if interpolate pixel is out of subject image region, set to -inf�����ֵ���س�������ͼ����������Ϊ-inf
					if (cur_pos_1[0]<0 || cur_pos_1[0]>sz_img_sub[0] - 1 ||
						cur_pos_1[1]<0 || cur_pos_1[1]>sz_img_sub[1] - 1 ||
						cur_pos_1[2]<0 || cur_pos_1[2]>sz_img_sub[2] - 1)
					{
						p_img_sub2tar_4d[0][z][y][x] = 0.0;
						continue;
					}
					if (cur_pos[0]<0 || cur_pos[0]>sz_img_read_size[0] - 1 ||
						cur_pos[1]<0 || cur_pos[1]>sz_img_read_size[1] - 1 ||
						cur_pos[2]<0 || cur_pos[2]>sz_img_read_size[2] - 1)
					{
						p_img_sub2tar_4d[0][z][y][x] = 0.0;
						continue;
					}
					//printf("mark32;(%d,%d,%d)\n", x, y, z);
					//find 8 neighbor pixels boundary�ҳ�8���������صı߽�
					long long x_s, x_b, y_s, y_b, z_s, z_b;
					x_s = floor(cur_pos[0]);		x_b = ceil(cur_pos[0]);
					y_s = floor(cur_pos[1]);		y_b = ceil(cur_pos[1]);
					z_s = floor(cur_pos[2]);		z_b = ceil(cur_pos[2]);
					//printf("(%d,%d,%d,%d,%d,%d)\n", x_s, y_s, z_s, x_b, y_b, z_b);
					//compute weight for left and right, top and bottom -- 4 neighbor pixel's weight in a slice�������ҡ��ϡ��µ�Ȩ�ء���һ����Ƭ��4���ھ����ص�Ȩ��
					double l_w, r_w, t_w, b_w;
					l_w = 1.0 - (cur_pos[0] - x_s);	r_w = 1.0 - l_w;//x_d
					t_w = 1.0 - (cur_pos[1] - y_s);	b_w = 1.0 - t_w;//y_d
					//compute weight for higher slice and lower slice����ϸ���Ƭ�ͽϵ���Ƭ��Ȩ��
					double u_w, d_w;
					u_w = 1.0 - (cur_pos[2] - z_s);	d_w = 1.0 - u_w;//z_d
					//printf("mark33;(%d,%d,%d)\n", x, y, z);
					//linear interpolate each channel
					//printf("marked2");
					for (long long c = 0; c < sz_img_output[3]; c++)
					{
						//linear interpolate in higher slice [t_w*(l_w*lt+r_w*rt)+b_w*(l_w*lb+r_w*rb)]
						double higher_slice;
						higher_slice = t_w*(l_w*p_img_sub_4d[c][z_s][y_s][x_s] + r_w*p_img_sub_4d[c][z_s][y_s][x_b]) +
							b_w*(l_w*p_img_sub_4d[c][z_s][y_b][x_s] + r_w*p_img_sub_4d[c][z_s][y_b][x_b]);
						//linear interpolate in lower slice [t_w*(l_w*lt+r_w*rt)+b_w*(l_w*lb+r_w*rb)]
						double lower_slice;
						lower_slice = t_w*(l_w*p_img_sub_4d[c][z_b][y_s][x_s] + r_w*p_img_sub_4d[c][z_b][y_s][x_b]) +
							b_w*(l_w*p_img_sub_4d[c][z_b][y_b][x_s] + r_w*p_img_sub_4d[c][z_b][y_b][x_b]);
						//linear interpolate the current position [u_w*higher_slice+d_w*lower_slice]
						p_img_sub2tar_4d[c][z][y][x] = u_w*higher_slice + d_w*lower_slice;
					}
					//printf("marked3");
					//printf("mark34;(%d,%d,%d)\n", x, y, z);
				}
		}
	}
	else if (type_gpu == 1)
	{
		//GPU����
		long long gsA2 = sz_img_read_size[2];
		long long gsA1 = sz_img_read_size[1];
		long long gsA0 = sz_img_read_size[0];//����Ƕ�ȡ�Ŀ�Ĵ�Сread
		/*
		long long gsz2 = sz_img_output[2];
		long long gsz1 = sz_img_output[1];
		long long gsz0 = sz_img_output[0];//��������ɵı����Ŀ�Ĵ�С������ֵ��
		*/
		long long gsz2 = sz_img_sort[2];
		long long gsz1 = sz_img_sort[1];
		long long gsz0 = sz_img_sort[0];//��������ɵı����Ŀ�Ĵ�С����������

		long long x_read_offset = sz_img_read_block[0];
		long long y_read_offset = sz_img_read_block[1];
		long long z_read_offset = sz_img_read_block[2];//����Ƕ�ȡ�Ŀ��ƫ��

		long long gs_ori2 = sz_img_sub[2];
		long long gs_ori1 = sz_img_sub[1];
		long long gs_ori0 = sz_img_sub[0];//�����ori

		gpu_interpolation_new(cal_mode, gsz2, gsz1, gsz0, x4x4_affinematrix, p_img_sub_4d, p_img_sub2tar_4d, sz_img_output, gsA2, gsA1, gsA0, x_offset, y_offset, z_offset, x_read_offset, y_read_offset, z_read_offset, gs_ori2, gs_ori1, gs_ori0, p_img_sub, p_img_affine);
	}

	return true;
}


//�ֿ�warp+ƴ�ӣ���Ҫ������
//bool q_imagewarp_affine_cyf(const vector<Coord3D_PCM> &vec_ctlpt_tar, const vector<Coord3D_PCM>  &vec_ctlpt_sub, char* &imgSrcFile, int type_gpu)
bool q_imagewarp_affine_cyf(const QList<ImageMarker> &ql_marker_tar, const QList<ImageMarker> &ql_marker_sub, string *path_TIF, long long sz_img_ori[4], string TeraWarp_save_base, int type_gpu)
{
	//check parameters
	
	if (ql_marker_tar.size() == 0 || ql_marker_sub.size() == 0 || ql_marker_tar.size() != ql_marker_sub.size())
	{
		printf("ERROR: target or subject control points is invalid!\n");
		return false;
	}
	//------------------------------------------------------------------------------------------------------------------------------------
	//re-formate to vector
	vector<Coord3D_PCM> vec_tar, vec_sub;
	long l_minlength = min(ql_marker_tar.size(), ql_marker_sub.size());

	// test
	int x_power, y_power, z_power;
	x_power = 4;// 4;
	y_power = 4;//4;
	z_power = 4;//8;
	for (long i = 0; i < l_minlength; i++)
	{
		vec_tar.push_back(Coord3D_PCM(x_power * ql_marker_tar[i].x, y_power * ql_marker_tar[i].y, z_power * ql_marker_tar[i].z));//22.9.6�޸�;Ϊ����Ӧ����������ͼ��
		vec_sub.push_back(Coord3D_PCM(x_power * ql_marker_sub[i].x, y_power * ql_marker_sub[i].y, z_power * ql_marker_sub[i].z));
	}

	// inference
//    for (long i = 0; i < l_minlength; i++)
//    {
//        vec_tar.push_back(Coord3D_PCM(ql_marker_tar[i].x, ql_marker_tar[i].y, ql_marker_tar[i].z));
//        vec_sub.push_back(Coord3D_PCM(ql_marker_sub[i].x, ql_marker_sub[i].y, ql_marker_sub[i].z));
//    }


	//------------------------------------------------------------------------------------------------------------------------------------
	//estimate the affine matrix
	Matrix x4x4_affinematrix;

	if (!q_affine_compute_affinmatrix_3D(vec_tar, vec_sub, x4x4_affinematrix))	//B=T*A
	{
		printf("ERROR: q_affine_compute_affinmatrix_2D() return false.\n");
		return false;
	}

	//------------------------------------------------------------�ֿ�warp---------------------------------------------------------------
	//��������ѭ����ȡCYF2022.6.15
	long long sz_img_block_size[4] = { 100, 100, 100, 1 };
	long long sz_img_block_size_sort[4] = { 100, 100, 100, 1 };
	long long *sz_img_read_block = 0;
	long long *sz_img_read_end_block = 0;

	int X_SIZE, Y_SIZE, Z_SIZE;
	long long per_block_siza =  1000;//����һ warp���ɵĿ�Ĵ�С
	long long per_sz_block = 2500;//������ ��ȡ�Ŀ�Ĵ�С

	int x_th = 0;
	int read_offset = 10;
	long long Z_Slice_control = 500;//����ÿһ��ƴ��ʱ�Ĳ���
	//control
	int warp_type = 1;
	int type_PJ = 1;

    std::stringstream base_block_path;
    std::stringstream slice_path;
    slice_path << TeraWarp_save_base << "_affine_slice";
    std::stringstream block_path;
    block_path << TeraWarp_save_base << "_affine_block";
	std::stringstream abs_pos_z;
	std::stringstream abs_pos_y;
	std::stringstream abs_pos_x;

//    if (0 != access(block_path.str().c_str(), 0))
//    {
//        mkdir(block_path.str().c_str());
//    }
//    if (0 != access(slice_path.str().c_str(), 0))
//    {
//        mkdir(slice_path.str().c_str());
//    }
    MkdirWithPath(block_path.str().c_str());
    MkdirWithPath(slice_path.str().c_str());


	Z_SIZE = sz_img_ori[2] % per_block_siza == 0 ? sz_img_ori[2] / per_block_siza : sz_img_ori[2] / per_block_siza + 1;
	Y_SIZE = sz_img_ori[1] % per_block_siza == 0 ? sz_img_ori[1] / per_block_siza : sz_img_ori[1] / per_block_siza + 1;
	X_SIZE = sz_img_ori[0] % per_block_siza == 0 ? sz_img_ori[0] / per_block_siza : sz_img_ori[0] / per_block_siza + 1;
	int sum_SIZE = X_SIZE*Y_SIZE*Z_SIZE;
	long long *sz_offset_x = 0;//ƴ�Ӻ�������Ҫ������,2022.8.29�Ѿ����Թ���������ʹ��
	long long *sz_offset_y = 0;
	sz_offset_x = new long long[X_SIZE*Y_SIZE*Z_SIZE];
	sz_offset_y = new long long[X_SIZE*Y_SIZE*Z_SIZE];

	//int real_type = 0;
	int cal_mode=0;
	
	clock_t aff_warp;
	clock_t read_time;
	clock_t save_time;
	clock_t PJ_time;

	float read_sum = 0;
	float aff_sum = 0;
	float save_sum = 0;
	vector<long long> sort_point_x;
	vector<long long> sort_point_y;
	vector<long long> sort_point_z;

	long long z_offset = 0;
	long long y_offset = 0;
	long long x_offset = 0;
	//string block_path = "E:/STPS/Outpu_t/2d_tif/2d_tif_200";
	//std::string block_path = "D:/BIG_Warp_File/OUTPUT/2D_raw_18052_1000_0920";
	//std::string block_path = "E:/STPS/Outpu_t/2d_raw_800_18052_18_0914";
	
	if (warp_type)
	{
		for (z_offset = 0; z_offset < sz_img_ori[2]; z_offset = z_offset + per_block_siza)
			for (y_offset = 0; y_offset < sz_img_ori[1]; y_offset = y_offset + per_block_siza)
			{
				for (x_offset = 0; x_offset < sz_img_ori[0]; x_offset = x_offset + per_block_siza)
					//x_th,y_thƴ�Ӻ�����Ҫ�����������22.8.29
				{
					sz_img_block_size[0] = (x_offset + per_block_siza) >= sz_img_ori[0] ? (sz_img_ori[0] - x_offset) : per_block_siza;
					sz_img_block_size[1] = (y_offset + per_block_siza) >= sz_img_ori[1] ? (sz_img_ori[1] - y_offset) : per_block_siza;
					sz_img_block_size[2] = (z_offset + per_block_siza) >= sz_img_ori[2] ? (sz_img_ori[2] - z_offset) : per_block_siza;

					//���ĵ�
					long long start_point[3];
					long long low_right_point[3];
					long long low_left_point[3];
					long long low_end_point[3];
					//���ĵ�
					long long high_right_point[3];
					long long high_left_point[3];
					long long high_start_point[3];
					long long end_point[3];

					//���ĵ�
					long long *read_start_point = 0;
					long long *read_low_right_point = 0;
					long long *read_low_left_point = 0;
					long long *read_low_end_point = 0;
					//���ĵ�
					long long *read_high_right_point = 0;
					long long *read_high_left_point = 0;
					long long *read_high_start_point = 0;
					long long *read_end_point = 0;

					//���ĵ�
					read_start_point = new long long[3];
					read_low_right_point = new long long[3];
					read_low_left_point = new long long[3];
					read_low_end_point = new long long[3];
					//���ĵ�
					read_high_right_point = new long long[3];
					read_high_left_point = new long long[3];
					read_high_start_point = new long long[3];
					read_end_point = new long long[3];


					//-------------------------------------����ת��--------------------------------------
					//��ʼ��
					start_point[0] = x_offset;
					start_point[1] = y_offset;
					start_point[2] = z_offset;
					if (!conversion_point(start_point, read_start_point, x4x4_affinematrix, 1))
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}

					//low_right
					low_right_point[0] = sz_img_block_size[0] + x_offset;
					low_right_point[1] = y_offset;
					low_right_point[2] = z_offset;
					if (!conversion_point(low_right_point, read_low_right_point, x4x4_affinematrix, 0))
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}

					//low_left
					low_left_point[0] = x_offset;
					low_left_point[1] = sz_img_block_size[1] + y_offset;
					low_left_point[2] = z_offset;
					if (!conversion_point(low_left_point, read_low_left_point, x4x4_affinematrix, 0))
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}

					//low_end
					low_end_point[0] = sz_img_block_size[0] + x_offset;
					low_end_point[1] = sz_img_block_size[1] + y_offset;
					low_end_point[2] = z_offset;
					if (!conversion_point(low_end_point, read_low_end_point, x4x4_affinematrix, 0))
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}

					//high_start
					high_start_point[0] = x_offset;
					high_start_point[1] = y_offset;
					high_start_point[2] = sz_img_block_size[2] + z_offset;
					if (!conversion_point(high_start_point, read_high_start_point, x4x4_affinematrix, 0))
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}

					//high_right
					high_right_point[0] = sz_img_block_size[0] + x_offset;
					high_right_point[1] = y_offset;
					high_right_point[2] = sz_img_block_size[2] + z_offset;
					if (!conversion_point(high_right_point, read_high_right_point, x4x4_affinematrix, 0))
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}

					//high_left
					high_left_point[0] = x_offset;
					high_left_point[1] = sz_img_block_size[1] + y_offset;
					high_left_point[2] = sz_img_block_size[2] + z_offset;
					if (!conversion_point(high_left_point, read_high_left_point, x4x4_affinematrix, 0))
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}

					//�յ�
					end_point[0] = sz_img_block_size[0] + x_offset;
					end_point[1] = sz_img_block_size[1] + y_offset;
					end_point[2] = sz_img_block_size[2] + z_offset;
					if (!conversion_point(end_point, read_end_point, x4x4_affinematrix, 0))
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}
					//-----------------------------------------------------------------------------------

					printf("warp����ʼ�����꣺%d,%d,%d\n", start_point[0], start_point[1], start_point[2]);
					printf("warp�Ľ��������꣺%d,%d,%d\n", end_point[0], end_point[1], end_point[2]);

					printf("��ȡ����ʼ�����꣺%d,%d,%d\n", read_start_point[0], read_start_point[1], read_start_point[2]);
					printf("��ȡ��low_right�����꣺%d,%d,%d\n", read_low_right_point[0], read_low_right_point[1], read_low_right_point[2]);
					printf("��ȡ��low_left�����꣺%d,%d,%d\n", read_low_left_point[0], read_low_left_point[1], read_low_left_point[2]);
					printf("��ȡ��low_end�����꣺%d,%d,%d\n", read_low_end_point[0], read_low_end_point[1], read_low_end_point[2]);

					printf("��ȡ��high_start�����꣺%d,%d,%d\n", read_high_start_point[0], read_high_start_point[1], read_high_start_point[2]);
					printf("��ȡ��high_right�����꣺%d,%d,%d\n", read_high_right_point[0], read_high_right_point[1], read_high_right_point[2]);
					printf("��ȡ��high_left�����꣺%d,%d,%d\n", read_high_left_point[0], read_high_left_point[1], read_high_left_point[2]);
					printf("��ȡ�Ľ��������꣺%d,%d,%d\n", read_end_point[0], read_end_point[1], read_end_point[2]);

					//--------------------------------------Sort-------------------------------------------
					sort_point_x.push_back(read_start_point[0]);
					sort_point_x.push_back(read_low_right_point[0]);
					sort_point_x.push_back(read_low_left_point[0]);
					sort_point_x.push_back(read_low_end_point[0]);
					sort_point_x.push_back(read_high_start_point[0]);
					sort_point_x.push_back(read_high_right_point[0]);
					sort_point_x.push_back(read_high_left_point[0]);
					sort_point_x.push_back(read_end_point[0]);

					sort(sort_point_x.begin(), sort_point_x.end());
//					printf("x����С��%d\n", sort_point_x.front());
//					printf("x�����%d\n", sort_point_x.back());
					//bubbleSort(sort_point_x);

					sort_point_y.push_back(read_start_point[1]);
					sort_point_y.push_back(read_low_right_point[1]);
					sort_point_y.push_back(read_low_left_point[1]);
					sort_point_y.push_back(read_low_end_point[1]);
					sort_point_y.push_back(read_high_start_point[1]);
					sort_point_y.push_back(read_high_right_point[1]);
					sort_point_y.push_back(read_high_left_point[1]);
					sort_point_y.push_back(read_end_point[1]);

					sort(sort_point_y.begin(), sort_point_y.end());
//					printf("y����С��%d\n", sort_point_y.front());
//					printf("y�����%d\n", sort_point_y.back());
					//bubbleSort(sort_point_y);

					sort_point_z.push_back(read_start_point[2]);
					sort_point_z.push_back(read_low_right_point[2]);
					sort_point_z.push_back(read_low_left_point[2]);
					sort_point_z.push_back(read_low_end_point[2]);
					sort_point_z.push_back(read_high_start_point[2]);
					sort_point_z.push_back(read_high_right_point[2]);
					sort_point_z.push_back(read_high_left_point[2]);
					sort_point_z.push_back(read_end_point[2]);

					sort(sort_point_z.begin(), sort_point_z.end());
//					printf("z����С��%d\n", sort_point_z.front());
//					printf("z�����%d\n", sort_point_z.back());
					//------------------------------------------------------------------------------------

				
					sz_img_read_block = new long long[4];
					//�����С�ĵ�
					sz_img_read_block[0] = (sort_point_x.front() - read_offset) < 0 ? 0 : (sort_point_x.front() - read_offset);
					sz_img_read_block[1] = (sort_point_y.front() - read_offset) < 0 ? 0 : (sort_point_y.front() - read_offset);
					sz_img_read_block[2] = (sort_point_z.front() - read_offset) < 0 ? 0 : (sort_point_z.front() - read_offset);

					sz_img_read_block[3] = 1;
					//������ĵ�
					sz_img_read_end_block = new long long[4];
					sz_img_read_end_block[0] = (sort_point_x.back() + read_offset) > sz_img_ori[0] ? sz_img_ori[0] : (sort_point_x.back() + read_offset);
					sz_img_read_end_block[1] = (sort_point_y.back() + read_offset) > sz_img_ori[1] ? sz_img_ori[1] : (sort_point_y.back() + read_offset);
					sz_img_read_end_block[2] = (sort_point_z.back() + read_offset) > sz_img_ori[2] ? sz_img_ori[2] : (sort_point_z.back() + read_offset);

					sz_img_read_end_block[3] = 1;
					//��Ҫ����ľ�������ڴ�����
					unsigned char *p_img_sub2tar_affine = 0;
					unsigned char *p_img_read_block = 0;
					unsigned char ****p_img_read_4d_block = 0;
					long long *sz_img_read_size = 0;
					unsigned char ****p_img_sub2tar_4d = 0;
					sz_img_read_size = new long long[4];

					p_img_sub2tar_affine = new unsigned char[sz_img_block_size[0] * sz_img_block_size[1] * sz_img_block_size[2] * sz_img_block_size[3]]();


					if (!new4dpointer(p_img_sub2tar_4d, sz_img_block_size[0], sz_img_block_size[1], sz_img_block_size[2], sz_img_block_size[3], p_img_sub2tar_affine))
					{
						printf("ERROR: Fail to allocate memory for the 4d pointer of image.\n");
						if (p_img_sub2tar_affine) 		{ delete[]p_img_sub2tar_affine;		p_img_sub2tar_affine = 0; }
						//if (p_img_sub_4d) 		{ delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
						if (p_img_sub2tar_4d) 	{ delete4dpointer(p_img_sub2tar_4d, sz_img_block_size[0], sz_img_block_size[1], sz_img_block_size[2], sz_img_block_size[3]); }
						return false;
					}
					//printf("mark1");
					read_time = clock();
					if (!read_2Dtif_BLOCK(path_TIF, p_img_read_block, p_img_read_4d_block, sz_img_read_block, sz_img_ori, sz_img_read_size, per_sz_block, sz_img_read_end_block))//sz_img_read_size
					{
						printf("ERROR: Fail to read_2Dtif_BLOCK.\n");
						//if (p_img_read_4d_block) 	{ delete4dpointer(p_img_read_4d_block, sz0, sz1, sz2, sz3); }
						return false;
					}
					read_sum += (float)(clock() - read_time) / CLOCKS_PER_SEC;//readʱ�����
					printf("\t>>read_time time consume %.2f s\n", (float)(clock() - read_time) / CLOCKS_PER_SEC);
					//printf("mark2\n");
					printf("warp�Ŀ�����꣺%d;%d;%d\n", x_offset, y_offset, z_offset);//warp�Ŀ������
					printf("warp�Ŀ����ʵ�ߴ磺%d;%d;%d\n", sz_img_block_size[0], sz_img_block_size[1], sz_img_block_size[2]);//warp�Ŀ�ĳߴ�
					//printf("warp�Ŀ������ߴ磺%d;%d;%d\n", sz_img_block_size_fake[0], sz_img_block_size_fake[1], sz_img_block_size_fake[2]);//warp�Ŀ�ĳߴ�
					printf("��ȡ�Ŀ�����꣺%d;%d;%d\n", sz_img_read_block[0], sz_img_read_block[1], sz_img_read_block[2]);//��ȡ�Ŀ������
					printf("��ȡ�Ŀ�ĳߴ磺%d;%d;%d\n", sz_img_read_size[0], sz_img_read_size[1], sz_img_read_size[2]);//��ȡ�Ŀ�ĳߴ�


                    sz_img_block_size_sort[0] = sz_img_block_size[0];
                    sz_img_block_size_sort[1] = sz_img_block_size[1];
                    sz_img_block_size_sort[2] = sz_img_block_size[2];

					aff_warp = clock();
					if (!q_imagewarp_affine_chazhi_CYF(type_gpu, sz_img_block_size, x4x4_affinematrix, x_offset, y_offset, z_offset, sz_img_ori, p_img_read_4d_block, p_img_sub2tar_4d,
						sz_img_read_block, sz_img_read_size, sz_img_block_size_sort, cal_mode, p_img_read_block, p_img_sub2tar_affine))
					{
						printf("ERROR: chazhi is wrong.\n");
						//if (p_img_read_4d_block) 	{ delete4dpointer(p_img_read_4d_block, sz0, sz1, sz2, sz3); }
						return false;
					}
					aff_sum += (float)(clock() - aff_warp) / CLOCKS_PER_SEC;
					printf("\t>>aff_tps time consume %.2f s\n", (float)(clock() - aff_warp) / CLOCKS_PER_SEC);
					//printf("mark3\n");
					abs_pos_z.width(6);
					abs_pos_y.width(6);
					abs_pos_x.width(6);

					abs_pos_z.fill('0');
					abs_pos_y.fill('0');
					abs_pos_x.fill('0');
					abs_pos_z << (int)(z_offset * 10);
					abs_pos_y << (int)(y_offset * 10);
					abs_pos_x << (int)(x_offset * 10);
					base_block_path << TeraWarp_save_base << "_affine_block" << "/" << abs_pos_z.str() << "_" << abs_pos_y.str() << "_" << abs_pos_x.str() << ".v3draw";
					sz_offset_x[x_th] = x_offset;
					sz_offset_y[x_th] = y_offset;

					x_th++;
					save_time = clock();

                    V3DLONG sz_img_block_size_[4] = {sz_img_block_size[0], sz_img_block_size[1], sz_img_block_size[2], sz_img_block_size[3]};
					saveImage(base_block_path.str().c_str(), p_img_sub2tar_affine, sz_img_block_size_, 1);
					save_sum += (float)(clock() - save_time) / CLOCKS_PER_SEC;
					printf("\t>>save time consume %.2f s\n", (float)(clock() - save_time) / CLOCKS_PER_SEC);
					//free memory
					string numStr;
					base_block_path.clear();
					base_block_path.str("");

					abs_pos_z.clear();
					abs_pos_z.str("");
					abs_pos_y.clear();
					abs_pos_y.str("");
					abs_pos_x.clear();
					abs_pos_x.str("");

					sort_point_x.clear();
					sort_point_y.clear();
					sort_point_z.clear();
					if (p_img_read_4d_block) 	{ delete4dpointer(p_img_read_4d_block, sz_img_read_size[0], sz_img_read_size[1], sz_img_read_size[2], sz_img_read_size[3]); }
					if (p_img_read_block) 				{ delete[]p_img_read_block;			p_img_read_block = 0; }
					if (p_img_sub2tar_4d) 	{ delete4dpointer(p_img_sub2tar_4d, sz_img_block_size[0], sz_img_block_size[1], sz_img_block_size[2], sz_img_block_size[3]); }
					if (p_img_sub2tar_affine) 				{ delete[]p_img_sub2tar_affine;			p_img_sub2tar_affine = 0; }
					if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
					if (read_low_right_point){ delete[]read_low_right_point; read_low_right_point = 0; }
					if (read_low_left_point){ delete[]read_low_left_point; read_low_left_point = 0; }
					if (read_low_end_point){ delete[]read_low_end_point; read_low_end_point = 0; }
					if (read_high_start_point){ delete[]read_high_start_point; read_high_start_point = 0; }
					if (read_high_right_point){ delete[]read_high_right_point; read_high_right_point = 0; }
					if (read_high_left_point){ delete[]read_high_left_point; read_high_left_point = 0; }
					if (read_end_point){ delete[]read_end_point; read_end_point = 0; }
					if (sz_img_read_block){ delete[]sz_img_read_block; sz_img_read_block = 0; }
					if (sz_img_read_size){ delete[]sz_img_read_size; sz_img_read_size = 0; }

				}
			}
	}

	


	PJ_time = clock();
	if (type_PJ)
	{
		q_imagewarp_stitch_2DRAW_new(block_path.str().c_str(), slice_path.str().c_str(), sz_img_ori, per_block_siza, sz_offset_x, sz_offset_y, Z_Slice_control);
	}
	printf("\t>>read all time consume %.2f s\n", read_sum);
	printf("\t>>stps all time consume %.2f s\n", aff_sum);
	printf("\t>>save all time consume %.2f s\n", save_sum);
	printf("\t>>PJ_time time consume %.2f s\n", (float)(clock() - PJ_time) / CLOCKS_PER_SEC);
	
//	if (path_TIF){ delete[]path_TIF; }
	if (sz_offset_x){ delete[]sz_offset_x; }
	if (sz_offset_y){ delete[]sz_offset_y; }

	

	//------------------------------------------------------------------------------------------------------------------------------------
	printf("6. free memory. \n");
    block_path.clear();
    block_path.str("");
    slice_path.clear();
    slice_path.str("");

	return true;
}


bool q_imagewarp_affine_image(const QList<ImageMarker> &ql_marker_tar, const QList<ImageMarker> &ql_marker_sub, unsigned char *p_img_sub, long long *sz_img_sub, string TeraWarp_save_base, int type_gpu)
{
    if (ql_marker_tar.size() == 0 || ql_marker_sub.size() == 0 || ql_marker_tar.size() != ql_marker_sub.size())
    {
        printf("ERROR: target or subject control points is invalid!\n");
        return false;
    }

    vector<Coord3D_PCM> vec_tar, vec_sub;
    long l_minlength = min(ql_marker_tar.size(), ql_marker_sub.size());

    for (long i = 0; i < l_minlength; i++)
    {
        vec_tar.push_back(Coord3D_PCM(ql_marker_tar[i].x, ql_marker_tar[i].y, ql_marker_tar[i].z));
        vec_sub.push_back(Coord3D_PCM(ql_marker_sub[i].x, ql_marker_sub[i].y, ql_marker_sub[i].z));
    }

    Matrix x4x4_affinematrix;

    if (!q_affine_compute_affinmatrix_3D(vec_tar, vec_sub, x4x4_affinematrix))	//B=T*A
    {
        printf("ERROR: q_affine_compute_affinmatrix_2D() return false.\n");
        return false;
    }

    long long sz_img_output[4] = { sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3] };
    unsigned char *p_img_affine = 0;
    unsigned char ****p_img_sub2tar_4d = 0;
    unsigned char ****p_img_sub_4d = 0;
    p_img_affine = new unsigned char[sz_img_output[0] * sz_img_output[1] * sz_img_output[2] * sz_img_output[3]]();

    if (!new4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3], p_img_sub) ||
        !new4dpointer(p_img_sub2tar_4d, sz_img_output[0], sz_img_output[1], sz_img_output[2], sz_img_output[3], p_img_affine))
    {
        printf("ERROR: Fail to allocate memory for the 4d pointer of image.\n");
        if (p_img_affine) 		{ delete[]p_img_affine;		p_img_affine = 0; }
        if (p_img_sub_4d) 		{ delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
        if (p_img_sub2tar_4d) 	{ delete4dpointer(p_img_sub2tar_4d, sz_img_output[0], sz_img_output[1], sz_img_output[2], sz_img_output[3]); }
        return false;
    }

    if (type_gpu == 0)
    {
        #pragma omp parallel for
        for (long long x = 0; x < sz_img_output[0]; x++)
        {
//            printf("affine: [%d/%d]\n", sz_img_output[0], x);
            for (long long y = 0; y < sz_img_output[1]; y++)
                for (long long z = 0; z < sz_img_output[2]; z++)
                {
                    Matrix x_pt_sub2tar_homo(4, 1), x_pt_sub_homo(4, 1);
                    //compute the inverse affine projected coordinate in subject image
                    x_pt_sub2tar_homo(1, 1) = x;
                    x_pt_sub2tar_homo(2, 1) = y;
                    x_pt_sub2tar_homo(3, 1) = z;
                    x_pt_sub2tar_homo(4, 1) = 1.0;
                    x_pt_sub_homo = x4x4_affinematrix*x_pt_sub2tar_homo;

                    //------------------------------------------------------------------
                    //linear interpolate
                    //coordinate in subject image
                    double cur_pos[3];//x,y,z
                    cur_pos[0] = x_pt_sub_homo(1, 1);
                    cur_pos[1] = x_pt_sub_homo(2, 1);
                    cur_pos[2] = x_pt_sub_homo(3, 1);

                    //if interpolate pixel is out of subject image region, set to -inf
                    if (cur_pos[0]<0 || cur_pos[0]>sz_img_sub[0] - 1 ||
                        cur_pos[1]<0 || cur_pos[1]>sz_img_sub[1] - 1 ||
                        cur_pos[2]<0 || cur_pos[2]>sz_img_sub[2] - 1)
                    {
                        p_img_sub2tar_4d[0][z][y][x] = 0.0;
                        continue;
                    }

                    //find 8 neighbor pixels boundary
                    long long x_s, x_b, y_s, y_b, z_s, z_b;
                    x_s = floor(cur_pos[0]);		x_b = ceil(cur_pos[0]);
                    y_s = floor(cur_pos[1]);		y_b = ceil(cur_pos[1]);
                    z_s = floor(cur_pos[2]);		z_b = ceil(cur_pos[2]);

                    //compute weight for left and right, top and bottom -- 4 neighbor pixel's weight in a slice
                    double l_w, r_w, t_w, b_w;
                    l_w = 1.0 - (cur_pos[0] - x_s);	r_w = 1.0 - l_w;
                    t_w = 1.0 - (cur_pos[1] - y_s);	b_w = 1.0 - t_w;
                    //compute weight for higher slice and lower slice
                    double u_w, d_w;
                    u_w = 1.0 - (cur_pos[2] - z_s);	d_w = 1.0 - u_w;

                    //linear interpolate each channel
                    for (long long c = 0; c < sz_img_output[3]; c++)
                    {
                        //linear interpolate in higher slice [t_w*(l_w*lt+r_w*rt)+b_w*(l_w*lb+r_w*rb)]
                        double higher_slice;
                        higher_slice = t_w*(l_w*p_img_sub_4d[c][z_s][y_s][x_s] + r_w*p_img_sub_4d[c][z_s][y_s][x_b]) +
                                       b_w*(l_w*p_img_sub_4d[c][z_s][y_b][x_s] + r_w*p_img_sub_4d[c][z_s][y_b][x_b]);
                        //linear interpolate in lower slice [t_w*(l_w*lt+r_w*rt)+b_w*(l_w*lb+r_w*rb)]
                        double lower_slice;
                        lower_slice = t_w*(l_w*p_img_sub_4d[c][z_b][y_s][x_s] + r_w*p_img_sub_4d[c][z_b][y_s][x_b]) +
                                      b_w*(l_w*p_img_sub_4d[c][z_b][y_b][x_s] + r_w*p_img_sub_4d[c][z_b][y_b][x_b]);
                        //linear interpolate the current position [u_w*higher_slice+d_w*lower_slice]
                        p_img_sub2tar_4d[c][z][y][x] = u_w*higher_slice + d_w*lower_slice;
                    }
                }
        }
    }
    else
    {
        long long gsA2 = sz_img_sub[2];
        long long gsA1 = sz_img_sub[1];
        long long gsA0 = sz_img_sub[0];

        long long gsz2 = sz_img_output[2];
        long long gsz1 = sz_img_output[1];
        long long gsz0 = sz_img_output[0];
//        printf("gsz2:%d,gsz1:%d,gsz0:%d\n", gsz2, gsz1, gsz0);
//        printf("gsA2:%d,gsA1:%d,gsA0:%d\n", gsA2, gsA1, gsA0);
//        printf("********************GPU***********************");
        //GPU???????
        clock_t aff_gpu_interpolation;
        aff_gpu_interpolation = clock();
        gpu_interpolation_affine(gsz2, gsz1, gsz0, x4x4_affinematrix, p_img_sub_4d, p_img_sub2tar_4d, sz_img_output, gsA2, gsA1, gsA0);//???:sz_img_sub??sz_img_output?????,?????????
        printf("\t>>aff_gpu_interpolation time consume %.2f s\n", (float)(clock() - aff_gpu_interpolation) / CLOCKS_PER_SEC);
    }

    V3DLONG sz_img_output_[4] = {sz_img_output[0], sz_img_output[1], sz_img_output[2], sz_img_output[3]};
    saveImage(TeraWarp_save_base.c_str(), p_img_affine, sz_img_output_, 1);

//------------------------------------------------------------------------------------------------------------------------------------
    printf("6. free memory. \n");
    if (p_img_sub_4d) 		{ delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
    if (p_img_sub2tar_4d) 	{ delete4dpointer(p_img_sub2tar_4d, sz_img_output[0], sz_img_output[1], sz_img_output[2], sz_img_output[3]); }
    if (p_img_affine) {delete[]p_img_affine; p_img_affine=0;}

    return true;


}




bool q_imagewarp_affine_cyf_downsample(const QList<ImageMarker> &ql_marker_tar, const QList<ImageMarker> &ql_marker_sub, string *path_TIF, long long sz_img_ori[4], string TeraWarp_save_base,int type_gpu)
{
	//check parameters

    long long factor = sz_img_ori[0] > sz_img_ori[1] ? max(sz_img_ori[0], sz_img_ori[2]) : max(sz_img_ori[1], sz_img_ori[2]);
    int downsample_factor = factor / 600 + (factor % 600 >= 300 ? 1 : 0);

	if (ql_marker_tar.size() == 0 || ql_marker_sub.size() == 0 || ql_marker_tar.size() != ql_marker_sub.size())
	{
		printf("ERROR: target or subject control points is invalid!\n");
		return false;
	}
	//------------------------------------------------------------------------------------------------------------------------------------
	//re-formate to vector
	vector<Coord3D_PCM> vec_tar, vec_sub;
	long l_minlength = min(ql_marker_tar.size(), ql_marker_sub.size());
	int x_power, y_power, z_power;
	int x_power_sub, y_power_sub, z_power_sub;
	x_power = 1;// 4;
	y_power = 1;//4;
	z_power = 1;//8;
    x_power_sub = downsample_factor;// 4;
    y_power_sub = downsample_factor;//4;
    z_power_sub = downsample_factor;//8;
	for (long i = 0; i < l_minlength; i++)
	{
		vec_tar.push_back(Coord3D_PCM(x_power * ql_marker_tar[i].x, y_power * ql_marker_tar[i].y, z_power * ql_marker_tar[i].z));//22.11.14�޸�;Ϊ����Ӧ������,�������Ļ�tarӦΪ���,subӦΪС��
		vec_sub.push_back(Coord3D_PCM(x_power_sub*ql_marker_sub[i].x, y_power_sub*ql_marker_sub[i].y, z_power_sub*ql_marker_sub[i].z));
	}//2022.11.14�����ǽ����˽����ģ�ʵ���ϴ�ĵ�Ӧ����tar��С�ĵ���sub��������ʱ�����˽�����������ⲿ����Ļ�����Ӧ�ð������Ǹ�������ǰ��������λ�õ���
	//------------------------------------------------------------------------------------------------------------------------------------
	//estimate the affine matrix
	Matrix x4x4_affinematrix;
	//clock_t affine_mat;
	//affine_mat = clock();
	if (!q_affine_compute_affinmatrix_3D(vec_tar, vec_sub, x4x4_affinematrix))	//B=T*A
	{
		printf("ERROR: q_affine_compute_affinmatrix_2D() return false.\n");
		return false;
	}
	//------------------------------------------------------------�ֿ�warp---------------------------------------------------------------
	//��������ѭ����ȡCYF2022.6.15




    long long sz_img_warp_ori[4] = { sz_img_ori[0] / downsample_factor, sz_img_ori[1] / downsample_factor, sz_img_ori[2] / downsample_factor, 1 };
	//long long sz_img_warp_ori[4] = { 547, 442, 682, 1 };//18052
	long long sz_img_block_size[4] = { 500, 500, 500, 1 };
	//long long sz_img_block_size_fake[4] = { 100, 100, 100, 1 };
	long long sz_img_block_size_sort[4] = { 100, 100, 100, 1 };
	long long *sz_img_read_block = 0;
	long long *sz_img_read_end_block = 0;

	int X_SIZE, Y_SIZE, Z_SIZE;



	long long per_block_siza = 1000;//����һ warp���ɵĿ�Ĵ�С
	long long per_sz_block = 2500;//������ ��ȡ�Ŀ�Ĵ�С

	int x_th = 0;
	int read_offset = 10;
	long long Z_Slice_control = 500;//����ÿһ��ƴ��ʱ�Ĳ���
	//control
	int warp_type = 1;
	int type_PJ = 1;



    std::stringstream base_block_path;
    std::stringstream slice_path;
    slice_path << TeraWarp_save_base << "_downsample_slice";
    std::stringstream block_path;
    block_path << TeraWarp_save_base << "_downsample_block";
    std::stringstream abs_pos_z;
    std::stringstream abs_pos_y;
    std::stringstream abs_pos_x;

//    if (0 != access(block_path.str().c_str(), 0))
//    {
//        mkdir(block_path.str().c_str());
//    }
//    if (0 != access(slice_path.str().c_str(), 0))
//    {
//        mkdir(slice_path.str().c_str());
//    }
    MkdirWithPath(block_path.str().c_str());
    MkdirWithPath(slice_path.str().c_str());


	Z_SIZE = sz_img_warp_ori[2] % per_block_siza == 0 ? sz_img_warp_ori[2] / per_block_siza : sz_img_warp_ori[2] / per_block_siza + 1;
	Y_SIZE = sz_img_warp_ori[1] % per_block_siza == 0 ? sz_img_warp_ori[1] / per_block_siza : sz_img_warp_ori[1] / per_block_siza + 1;
	X_SIZE = sz_img_warp_ori[0] % per_block_siza == 0 ? sz_img_warp_ori[0] / per_block_siza : sz_img_warp_ori[0] / per_block_siza + 1;
	int sum_SIZE = X_SIZE*Y_SIZE*Z_SIZE;
	long long *sz_offset_x = 0;//ƴ�Ӻ�������Ҫ������,2022.8.29�Ѿ����Թ���������ʹ��
	long long *sz_offset_y = 0;
	sz_offset_x = new long long[X_SIZE*Y_SIZE*Z_SIZE];
	sz_offset_y = new long long[X_SIZE*Y_SIZE*Z_SIZE];

	//int real_type = 0;
	int cal_mode;

	clock_t aff_warp;
	clock_t read_time;
	clock_t save_time;
	clock_t PJ_time;

	float read_sum = 0;
	float aff_sum = 0;
	float save_sum = 0;
	vector<long long> sort_point_x;
	vector<long long> sort_point_y;
	vector<long long> sort_point_z;

	long long z_offset = 0;
	long long y_offset = 0;
	long long x_offset = 0;

	//-------------------------------------------------------------------------------------------------------------------------------
	if (warp_type)
	{
		for (z_offset = 0; z_offset < sz_img_warp_ori[2]; z_offset = z_offset + per_block_siza)
			for (y_offset = 0; y_offset < sz_img_warp_ori[1]; y_offset = y_offset + per_block_siza)
			{
				for (x_offset = 0; x_offset < sz_img_warp_ori[0]; x_offset = x_offset + per_block_siza)
					//x_th,y_thƴ�Ӻ�����Ҫ�����������22.8.29
				{
					sz_img_block_size[0] = (x_offset + per_block_siza) >= sz_img_warp_ori[0] ? (sz_img_warp_ori[0] - x_offset) : per_block_siza;
					sz_img_block_size[1] = (y_offset + per_block_siza) >= sz_img_warp_ori[1] ? (sz_img_warp_ori[1] - y_offset) : per_block_siza;
					sz_img_block_size[2] = (z_offset + per_block_siza) >= sz_img_warp_ori[2] ? (sz_img_warp_ori[2] - z_offset) : per_block_siza;

					/*long long centre_point[3];
					long long start_point[3];
					long long end_point[3];
					long long *read_centre_point = 0;
					long long *read_start_point = 0;
					long long *read_end_point = 0;*/

					//���ĵ�
					long long start_point[3];
					long long low_right_point[3];
					long long low_left_point[3];
					long long low_end_point[3];
					//���ĵ�
					long long high_right_point[3];
					long long high_left_point[3];
					long long high_start_point[3];
					long long end_point[3];

					//���ĵ�
					long long *read_start_point = 0;
					long long *read_low_right_point = 0;
					long long *read_low_left_point = 0;
					long long *read_low_end_point = 0;
					//���ĵ�
					long long *read_high_right_point = 0;
					long long *read_high_left_point = 0;
					long long *read_high_start_point = 0;
					long long *read_end_point = 0;

					//���ĵ�
					read_start_point = new long long[3];
					read_low_right_point = new long long[3];
					read_low_left_point = new long long[3];
					read_low_end_point = new long long[3];
					//���ĵ�
					read_high_right_point = new long long[3];
					read_high_left_point = new long long[3];
					read_high_start_point = new long long[3];
					read_end_point = new long long[3];

					/*read_centre_point = new long long[3];
					read_start_point = new long long[3];
					read_end_point = new long long[3];*/
					//-------------------------------------����ת��--------------------------------------
					//��ʼ��
					start_point[0] = x_offset;
					start_point[1] = y_offset;
					start_point[2] = z_offset;
					if (!conversion_point(start_point, read_start_point, x4x4_affinematrix, 1))
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}

					//low_right
					low_right_point[0] = sz_img_block_size[0] + x_offset;
					low_right_point[1] = y_offset;
					low_right_point[2] = z_offset;
					if (!conversion_point(low_right_point, read_low_right_point, x4x4_affinematrix, 0))
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}

					//low_left
					low_left_point[0] = x_offset;
					low_left_point[1] = sz_img_block_size[1] + y_offset;
					low_left_point[2] = z_offset;
					if (!conversion_point(low_left_point, read_low_left_point, x4x4_affinematrix, 0))
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}

					//low_end
					low_end_point[0] = sz_img_block_size[0] + x_offset;
					low_end_point[1] = sz_img_block_size[1] + y_offset;
					low_end_point[2] = z_offset;
					if (!conversion_point(low_end_point, read_low_end_point, x4x4_affinematrix, 0))
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}

					//high_start
					high_start_point[0] = x_offset;
					high_start_point[1] = y_offset;
					high_start_point[2] = sz_img_block_size[2] + z_offset;
					if (!conversion_point(high_start_point, read_high_start_point, x4x4_affinematrix, 0))
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}

					//high_right
					high_right_point[0] = sz_img_block_size[0] + x_offset;
					high_right_point[1] = y_offset;
					high_right_point[2] = sz_img_block_size[2] + z_offset;
					if (!conversion_point(high_right_point, read_high_right_point, x4x4_affinematrix, 0))
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}

					//high_left
					high_left_point[0] = x_offset;
					high_left_point[1] = sz_img_block_size[1] + y_offset;
					high_left_point[2] = sz_img_block_size[2] + z_offset;
					if (!conversion_point(high_left_point, read_high_left_point, x4x4_affinematrix, 0))
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}

					//�յ�
					end_point[0] = sz_img_block_size[0] + x_offset;
					end_point[1] = sz_img_block_size[1] + y_offset;
					end_point[2] = sz_img_block_size[2] + z_offset;
					if (!conversion_point(end_point, read_end_point, x4x4_affinematrix, 0))
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}
					//-----------------------------------------------------------------------------------

//					printf("warp����ʼ�����꣺%d,%d,%d\n", start_point[0], start_point[1], start_point[2]);
//					printf("warp�Ľ��������꣺%d,%d,%d\n", end_point[0], end_point[1], end_point[2]);
//
//					printf("��ȡ����ʼ�����꣺%d,%d,%d\n", read_start_point[0], read_start_point[1], read_start_point[2]);
//					printf("��ȡ��low_right�����꣺%d,%d,%d\n", read_low_right_point[0], read_low_right_point[1], read_low_right_point[2]);
//					printf("��ȡ��low_left�����꣺%d,%d,%d\n", read_low_left_point[0], read_low_left_point[1], read_low_left_point[2]);
//					printf("��ȡ��low_end�����꣺%d,%d,%d\n", read_low_end_point[0], read_low_end_point[1], read_low_end_point[2]);
//
//					printf("��ȡ��high_start�����꣺%d,%d,%d\n", read_high_start_point[0], read_high_start_point[1], read_high_start_point[2]);
//					printf("��ȡ��high_right�����꣺%d,%d,%d\n", read_high_right_point[0], read_high_right_point[1], read_high_right_point[2]);
//					printf("��ȡ��high_left�����꣺%d,%d,%d\n", read_high_left_point[0], read_high_left_point[1], read_high_left_point[2]);
//					printf("��ȡ�Ľ��������꣺%d,%d,%d\n", read_end_point[0], read_end_point[1], read_end_point[2]);

					//--------------------------------------Sort-------------------------------------------
					sort_point_x.push_back(read_start_point[0]);
					sort_point_x.push_back(read_low_right_point[0]);
					sort_point_x.push_back(read_low_left_point[0]);
					sort_point_x.push_back(read_low_end_point[0]);
					sort_point_x.push_back(read_high_start_point[0]);
					sort_point_x.push_back(read_high_right_point[0]);
					sort_point_x.push_back(read_high_left_point[0]);
					sort_point_x.push_back(read_end_point[0]);

					sort(sort_point_x.begin(), sort_point_x.end());
//					printf("x����С��%d\n", sort_point_x.front());
//					printf("x�����%d\n", sort_point_x.back());
					//bubbleSort(sort_point_x);

					sort_point_y.push_back(read_start_point[1]);
					sort_point_y.push_back(read_low_right_point[1]);
					sort_point_y.push_back(read_low_left_point[1]);
					sort_point_y.push_back(read_low_end_point[1]);
					sort_point_y.push_back(read_high_start_point[1]);
					sort_point_y.push_back(read_high_right_point[1]);
					sort_point_y.push_back(read_high_left_point[1]);
					sort_point_y.push_back(read_end_point[1]);

					sort(sort_point_y.begin(), sort_point_y.end());
//					printf("y����С��%d\n", sort_point_y.front());
//					printf("y�����%d\n", sort_point_y.back());
					//bubbleSort(sort_point_y);

					sort_point_z.push_back(read_start_point[2]);
					sort_point_z.push_back(read_low_right_point[2]);
					sort_point_z.push_back(read_low_left_point[2]);
					sort_point_z.push_back(read_low_end_point[2]);
					sort_point_z.push_back(read_high_start_point[2]);
					sort_point_z.push_back(read_high_right_point[2]);
					sort_point_z.push_back(read_high_left_point[2]);
					sort_point_z.push_back(read_end_point[2]);

					sort(sort_point_z.begin(), sort_point_z.end());
//					printf("z����С��%d\n", sort_point_z.front());
//					printf("z�����%d\n", sort_point_z.back());


					sz_img_read_block = new long long[4];
					//�����С�ĵ�
					sz_img_read_block[0] = (sort_point_x.front() - read_offset) < 0 ? 0 : (sort_point_x.front() - read_offset);
					sz_img_read_block[1] = (sort_point_y.front() - read_offset) < 0 ? 0 : (sort_point_y.front() - read_offset);
					sz_img_read_block[2] = (sort_point_z.front() - read_offset) < 0 ? 0 : (sort_point_z.front() - read_offset);

					sz_img_read_block[3] = 1;
					//������ĵ�
					sz_img_read_end_block = new long long[4];
					sz_img_read_end_block[0] = (sort_point_x.back() + read_offset) > sz_img_ori[0] ? sz_img_ori[0] : (sort_point_x.back() + read_offset);
					sz_img_read_end_block[1] = (sort_point_y.back() + read_offset) > sz_img_ori[1] ? sz_img_ori[1] : (sort_point_y.back() + read_offset);
					sz_img_read_end_block[2] = (sort_point_z.back() + read_offset) > sz_img_ori[2] ? sz_img_ori[2] : (sort_point_z.back() + read_offset);

					sz_img_read_end_block[3] = 1;
					//��Ҫ����ľ�������ڴ�����
					unsigned char *p_img_sub2tar_affine = 0;
					unsigned char *p_img_read_block = 0;
					unsigned char ****p_img_read_4d_block = 0;
					long long *sz_img_read_size = 0;
					unsigned char ****p_img_sub2tar_4d = 0;
					sz_img_read_size = new long long[4];

					p_img_sub2tar_affine = new unsigned char[sz_img_block_size[0] * sz_img_block_size[1] * sz_img_block_size[2] * sz_img_block_size[3]]();


					if (!new4dpointer(p_img_sub2tar_4d, sz_img_block_size[0], sz_img_block_size[1], sz_img_block_size[2], sz_img_block_size[3], p_img_sub2tar_affine))
					{
						printf("ERROR: Fail to allocate memory for the 4d pointer of image.\n");
						if (p_img_sub2tar_affine) 		{ delete[]p_img_sub2tar_affine;		p_img_sub2tar_affine = 0; }
						//if (p_img_sub_4d) 		{ delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
						if (p_img_sub2tar_4d) 	{ delete4dpointer(p_img_sub2tar_4d, sz_img_block_size[0], sz_img_block_size[1], sz_img_block_size[2], sz_img_block_size[3]); }
						return false;
					}
					//printf("mark1");
					read_time = clock();
					if (!read_2Dtif_BLOCK(path_TIF, p_img_read_block, p_img_read_4d_block, sz_img_read_block, sz_img_ori, sz_img_read_size, per_sz_block, sz_img_read_end_block))//sz_img_read_size
					{
						printf("ERROR: Fail to read_2Dtif_BLOCK.\n");
						//if (p_img_read_4d_block) 	{ delete4dpointer(p_img_read_4d_block, sz0, sz1, sz2, sz3); }
						return false;
					}
					read_sum += (float)(clock() - read_time) / CLOCKS_PER_SEC;//readʱ�����
					printf("\t>>read_time time consume %.2f s\n", (float)(clock() - read_time) / CLOCKS_PER_SEC);
					//printf("mark2\n");
					printf("warp�Ŀ�����꣺%d;%d;%d\n", x_offset, y_offset, z_offset);//warp�Ŀ������
					printf("warp�Ŀ����ʵ�ߴ磺%d;%d;%d\n", sz_img_block_size[0], sz_img_block_size[1], sz_img_block_size[2]);//warp�Ŀ�ĳߴ�
					//printf("warp�Ŀ������ߴ磺%d;%d;%d\n", sz_img_block_size_fake[0], sz_img_block_size_fake[1], sz_img_block_size_fake[2]);//warp�Ŀ�ĳߴ�
					printf("��ȡ�Ŀ�����꣺%d;%d;%d\n", sz_img_read_block[0], sz_img_read_block[1], sz_img_read_block[2]);//��ȡ�Ŀ������
					printf("��ȡ�Ŀ�ĳߴ磺%d;%d;%d\n", sz_img_read_size[0], sz_img_read_size[1], sz_img_read_size[2]);//��ȡ�Ŀ�ĳߴ�

					if (sz_img_block_size[0] >= GPU_LIMIT && sz_img_block_size[1] >= GPU_LIMIT && sz_img_block_size[2] >= GPU_LIMIT)//X,Y,Z>0
					{
						type_gpu = 1;
						cal_mode = 0;
						sz_img_block_size_sort[0] = sz_img_block_size[0];
						sz_img_block_size_sort[1] = sz_img_block_size[1];
						sz_img_block_size_sort[2] = sz_img_block_size[2];

					}
					else if (sz_img_block_size[0] >= GPU_LIMIT && sz_img_block_size[1] >= GPU_LIMIT && sz_img_block_size[2] < GPU_LIMIT)//Z<0
					{
						type_gpu = 1;
						cal_mode = 0;
						sz_img_block_size_sort[0] = sz_img_block_size[0];
						sz_img_block_size_sort[1] = sz_img_block_size[1];
						sz_img_block_size_sort[2] = sz_img_block_size[2];
					}
					else if (sz_img_block_size[0] >= GPU_LIMIT && sz_img_block_size[1] < GPU_LIMIT && sz_img_block_size[2] >= GPU_LIMIT)//Y<0
					{
						type_gpu = 1;
						cal_mode = 1;
						sz_img_block_size_sort[0] = sz_img_block_size[0];
						sz_img_block_size_sort[1] = sz_img_block_size[2];
						sz_img_block_size_sort[2] = sz_img_block_size[1];
					}
					else if (sz_img_block_size[0] < GPU_LIMIT && sz_img_block_size[1] >= GPU_LIMIT && sz_img_block_size[2] >= GPU_LIMIT)//X<0
					{
						type_gpu = 1;
						cal_mode = 2;
						sz_img_block_size_sort[0] = sz_img_block_size[1];
						sz_img_block_size_sort[1] = sz_img_block_size[2];
						sz_img_block_size_sort[2] = sz_img_block_size[0];
					}
					else
					{
						type_gpu = 0;
						cal_mode = 0;
						sz_img_block_size_sort[0] = sz_img_block_size[0];
						sz_img_block_size_sort[1] = sz_img_block_size[1];
						sz_img_block_size_sort[2] = sz_img_block_size[2];
					}

					aff_warp = clock();
					if (!q_imagewarp_affine_chazhi_CYF(type_gpu, sz_img_block_size, x4x4_affinematrix, x_offset, y_offset, z_offset, sz_img_ori, p_img_read_4d_block, p_img_sub2tar_4d,
						sz_img_read_block, sz_img_read_size, sz_img_block_size_sort, cal_mode, p_img_read_block, p_img_sub2tar_affine))
					{
						printf("ERROR: chazhi is wrong.\n");
						//if (p_img_read_4d_block) 	{ delete4dpointer(p_img_read_4d_block, sz0, sz1, sz2, sz3); }
						return false;
					}
					aff_sum += (float)(clock() - aff_warp) / CLOCKS_PER_SEC;
					printf("\t>>aff_tps time consume %.2f s\n", (float)(clock() - aff_warp) / CLOCKS_PER_SEC);
					//printf("mark3\n");
					abs_pos_z.width(6);
					abs_pos_y.width(6);
					abs_pos_x.width(6);

					abs_pos_z.fill('0');
					abs_pos_y.fill('0');
					abs_pos_x.fill('0');
					abs_pos_z << (int)(z_offset * 10);
					abs_pos_y << (int)(y_offset * 10);
					abs_pos_x << (int)(x_offset * 10);
					//base_block_path << block_path << "/" << (int)z_offset << "x" << (int)y_offset << "x" << (int)x_offset << ".v3draw";
					//base_block_path << block_path << "/" << (int)x_th  << ".v3draw";
					base_block_path << TeraWarp_save_base << "_downsample_block" << "/" << abs_pos_z.str() << "_" << abs_pos_y.str() << "_" << abs_pos_x.str() << ".v3draw";
					//pinjie_need
					sz_offset_x[x_th] = x_offset;
					sz_offset_y[x_th] = y_offset;

					x_th++;
					save_time = clock();
					//printf("%s\n", base_block_path.str().c_str());
                    V3DLONG sz_img_block_size_[4] = {sz_img_block_size[0], sz_img_block_size[1], sz_img_block_size[2], sz_img_block_size[3]};
					saveImage(base_block_path.str().c_str(), p_img_sub2tar_affine, sz_img_block_size_, 1);
					save_sum += (float)(clock() - save_time) / CLOCKS_PER_SEC;
					printf("\t>>save time consume %.2f s\n", (float)(clock() - save_time) / CLOCKS_PER_SEC);
					//free memory
					string numStr;
					base_block_path.clear();
					base_block_path.str("");

					abs_pos_z.clear();
					abs_pos_z.str("");
					abs_pos_y.clear();
					abs_pos_y.str("");
					abs_pos_x.clear();
					abs_pos_x.str("");

					sort_point_x.clear();
					sort_point_y.clear();
					sort_point_z.clear();
					if (p_img_read_4d_block) 	{ delete4dpointer(p_img_read_4d_block, sz_img_read_size[0], sz_img_read_size[1], sz_img_read_size[2], sz_img_read_size[3]); }
					if (p_img_read_block) 				{ delete[]p_img_read_block;			p_img_read_block = 0; }
					if (p_img_sub2tar_4d) 	{ delete4dpointer(p_img_sub2tar_4d, sz_img_block_size[0], sz_img_block_size[1], sz_img_block_size[2], sz_img_block_size[3]); }
					if (p_img_sub2tar_affine) 				{ delete[]p_img_sub2tar_affine;			p_img_sub2tar_affine = 0; }
					if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
					if (read_low_right_point){ delete[]read_low_right_point; read_low_right_point = 0; }
					if (read_low_left_point){ delete[]read_low_left_point; read_low_left_point = 0; }
					if (read_low_end_point){ delete[]read_low_end_point; read_low_end_point = 0; }
					if (read_high_start_point){ delete[]read_high_start_point; read_high_start_point = 0; }
					if (read_high_right_point){ delete[]read_high_right_point; read_high_right_point = 0; }
					if (read_high_left_point){ delete[]read_high_left_point; read_high_left_point = 0; }
					if (read_end_point){ delete[]read_end_point; read_end_point = 0; }
					if (sz_img_read_block){ delete[]sz_img_read_block; sz_img_read_block = 0; }
					if (sz_img_read_size){ delete[]sz_img_read_size; sz_img_read_size = 0; }

				}
			}
	}



	//********************************
	//ƴ�Ӳ�������
	/*
	x_th = 0;
	for (long long z_offset = 0; z_offset < sz_img_ori[2]; z_offset = z_offset + per_block_siza)
	{
	for (long long y_offset = 0; y_offset < sz_img_ori[1]; y_offset = y_offset + per_block_siza)
	{
	for (long long x_offset = 0; x_offset < sz_img_ori[0]; x_offset = x_offset + per_block_siza)


	//x_th,y_thƴ�Ӻ�����Ҫ�����������22.8.29
	{
	sz_offset_x[x_th] = x_offset;
	sz_offset_y[x_th] = y_offset;
	x_th++;
	}
	}
	}
	for (int i = 0; i < X_SIZE*Y_SIZE*Z_SIZE; i++)
	{
	printf("X_OFFSET:%d\n", sz_offset_x[i]);
	printf("Y_OFFSET:%d\n", sz_offset_y[i]);
	}
	*/
	//********************************
	PJ_time = clock();
	if (type_PJ)
	{
		q_imagewarp_stitch_2DRAW_new(block_path.str().c_str(), slice_path.str().c_str(), sz_img_warp_ori, per_block_siza, sz_offset_x, sz_offset_y, Z_Slice_control);
	}
	printf("\t>>read all time consume %.2f s\n", read_sum);
	printf("\t>>stps all time consume %.2f s\n", aff_sum);
	printf("\t>>save all time consume %.2f s\n", save_sum);
	printf("\t>>PJ_time time consume %.2f s\n", (float)(clock() - PJ_time) / CLOCKS_PER_SEC);


	if (sz_offset_x){ delete[]sz_offset_x; }
	if (sz_offset_y){ delete[]sz_offset_y; }



	//------------------------------------------------------------------------------------------------------------------------------------
	printf("6. free memory. \n");
	//if (p_img_sub_4d) 		{ delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
	//if (p_img_sub2tar_4d) 	{ delete4dpointer(p_img_sub2tar_4d, sz_img_output[0], sz_img_output[1], sz_img_output[2], sz_img_output[3]); }

	return true;
}
