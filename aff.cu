#include "ele.cuh"
#include "cuda_runtime_api.h"
#include <malloc.h>
#include <time.h>
#include"newmat.h"
#include <math.h>
//#include "q_warp_affine_tps.h"
#include "cuda_runtime.h"
#include "cublas_v2.h"
#include "cusolverDn.h"
#include <assert.h>
#include <stdlib.h>
//#include <QtGui>

extern "C" Matrix matrixMultiply(const int m, const int n, const int k, Matrix &A, Matrix &B);

#define number 5000   //Ua_size in TPS
#define EPS 0.0001

__host__ __device__ void MUL_aff(float *A, float *B, float *C)//A*B=C
{
	for (int i = 0; i < 4; i++)
	{


		float sum = 0;
		for (int k = 0; k < 4; k++)
		{
			sum += A[i * 4 + k] * B[k];
		}
		C[i] = sum;

	}
}


extern "C" Matrix matrixMultiply(const int m, const int n, const int k, Matrix &A, Matrix &B)
{
	//A*B=C
	//m:A.row;n:B.col;k:A.col
	Matrix C(m, n);
	cudaError_t cudaStat;
	cublasStatus_t stat;
	float *H_A, *H_B, *H_C;

	float *D_A, *D_B, *D_C;

	H_A = (float*)malloc(m * k * sizeof(float));
	H_B = (float*)malloc(k * n * sizeof(float));
	H_C = (float*)malloc(m * n * sizeof(float));

	cudaStat = cudaMalloc((void**)&D_A, m * k * sizeof(float));
	cudaStat = cudaMalloc((void**)&D_B, k * n * sizeof(float));
	cudaStat = cudaMalloc((void**)&D_C, m * n * sizeof(float));

	/*	cudaStat = cudaMalloc((void**)&D_A, r_size * r_size * sizeof(float));
	cudaStat = cudaMalloc((void**)&D_B, r_size * r_size * sizeof(float));
	cudaStat = cudaMalloc((void**)&D_C, r_size * r_size * sizeof(float));*/
	printf("cudaStat %d\n", cudaStat);
	for (int i = 0; i < A.nrows(); i++)
	{
		for (int j = 0; j < A.ncols(); j++)
		{

			H_A[i * A.ncols() + j] = A(i + 1, j + 1);
			//if (H_A[i * A.ncols() + j] > EPS)printf("%f", H_A[i * A.ncols() + j]);
		}

	}

	for (int i = 0; i < B.nrows(); i++)
	{
		for (int j = 0; j < B.ncols(); j++)
		{

			H_B[i * B.ncols() + j] = B(i + 1, j + 1);
			//if (H_B[i * B.ncols() + j] > EPS)printf("%.2f\n", H_B[i * B.ncols() + j]);
		}

	}

	HANDLE_ERROR(cudaMemcpy(D_A, H_A, m * k * sizeof(float), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(D_B, H_B, k * n * sizeof(float), cudaMemcpyHostToDevice));

	//	stat = cublasSetMatrix(r_size, r_size, sizeof(*H_A), H_A, r_size, D_A, r_size);
	//	stat = cublasSetMatrix(r_size, r_size, sizeof(*H_B), H_B, r_size, D_B, r_size);
	//	stat = cublasSetMatrix(r_size, r_size, sizeof(*H_C), H_C, r_size, D_C, r_size);

	const float alpha = 1.0f;
	const float beta = 0.0f;

	cublasHandle_t handle;
	cublasCreate(&handle);

	stat = cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, n, m, k, &alpha, D_B, n, D_A, k, &beta, D_C, n);

	cublasDestroy(handle);

	printf("cublas %d\n", stat);
	HANDLE_ERROR(cudaMemcpy(H_C, D_C, m * n * sizeof(float), cudaMemcpyDeviceToHost));


	for (int i = 0; i < C.nrows(); i++)
	{
		for (int j = 0; j < C.ncols(); j++)
		{

			C(i + 1, j + 1) = H_C[i * C.ncols() + j];
			//if (C(i + 1, j + 1) > EPS)printf("%f\n", C(i + 1, j + 1));
		}

	}

	/*	for (int i = 0; i < C.nrows(); i++)
	{
	for (int j = 0; j < C.ncols(); j++)
	{

	C(i + 1, j + 1) = H_C[j * C.ncols() + i];

	}

	}*/

	free(H_A);
	free(H_B);
	free(H_C);
	cudaFree(D_A);
	cudaFree(D_B);
	cudaFree(D_C);
	HANDLE_ERROR(cudaDeviceReset());
	return C;
}


__global__ void get_Displacement_affine_new_Z(const int k, const long long gsz1, const long long gsz0, float *D_x4x4_affinematrix, const int D_sz_img_sub0,
	const int D_sz_img_sub1, const int D_sz_img_sub2, const int D_sz_img_sub3, unsigned char *D_p_img_sub_4d, unsigned char *D_p_img_sub2tar_4d, const int D_sz_img_A_sub0,
	const int D_sz_img_A_sub1, const int D_sz_img_A_sub2, const long long start_block_x, const long long start_block_y, const long long start_block_z,
	const long long x_read_offset, const long long y_read_offset, const long long z_read_offset, const long long gs_ori2, const long long gs_ori1, const long long gs_ori0)
{
	const int row = blockIdx.y * blockDim.y + threadIdx.y;
	const int col = blockIdx.x * blockDim.x + threadIdx.x;

	if (row >= gsz1 || col >= gsz0)return;


	float x_pt_sub2tar_homo[4];
	float x_pt_sub_homo[4];


	x_pt_sub2tar_homo[0] = col + start_block_x;
	x_pt_sub2tar_homo[1] = row + start_block_y;
	x_pt_sub2tar_homo[2] = k + start_block_z;
	x_pt_sub2tar_homo[3] = 1.0;

	MUL_aff(D_x4x4_affinematrix, x_pt_sub2tar_homo, x_pt_sub_homo);

	double cur_pos[3];//x,y,z
	double cur_pos_1[3];//x,y,z
	cur_pos[0] = x_pt_sub_homo[0] - x_read_offset;
	cur_pos[1] = x_pt_sub_homo[1] - y_read_offset;
	cur_pos[2] = x_pt_sub_homo[2] - z_read_offset;

	cur_pos_1[0] = x_pt_sub_homo[0];
	cur_pos_1[1] = x_pt_sub_homo[1];
	cur_pos_1[2] = x_pt_sub_homo[2];

	if (cur_pos_1[0]<0 || cur_pos_1[0]>gs_ori0 - 1 ||
		cur_pos_1[1]<0 || cur_pos_1[1]>gs_ori1 - 1 ||
		cur_pos_1[2]<0 || cur_pos_1[2]>gs_ori2 - 1 ||
		cur_pos[0]<0 || cur_pos[0]>D_sz_img_A_sub0 - 1 ||
		cur_pos[1]<0 || cur_pos[1]>D_sz_img_A_sub1 - 1 ||
		cur_pos[2]<0 || cur_pos[2]>D_sz_img_A_sub2 - 1)
	{
		for (long long c = 0; c < D_sz_img_sub3; c++)
		{
			D_p_img_sub2tar_4d[c * D_sz_img_sub2 * D_sz_img_sub1 * D_sz_img_sub0 + k * D_sz_img_sub1 * D_sz_img_sub0 + row * D_sz_img_sub0 + col] = 0.0;

		}

	}
	else{


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


		float a = t_w*(l_w*D_p_img_sub_4d[z_s*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_s*D_sz_img_A_sub0 + x_s] + r_w*D_p_img_sub_4d[z_s*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_s*D_sz_img_A_sub0 + x_b]) + b_w*(l_w*D_p_img_sub_4d[z_s*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_b*D_sz_img_A_sub0 + x_s] + r_w*D_p_img_sub_4d[z_s*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_b*D_sz_img_A_sub0 + x_b]);

		float b = t_w*(l_w*D_p_img_sub_4d[z_b*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_s*D_sz_img_A_sub0 + x_s] + r_w*D_p_img_sub_4d[z_b*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_s*D_sz_img_A_sub0 + x_b]) + b_w*(l_w*D_p_img_sub_4d[z_b*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_b*D_sz_img_A_sub0 + x_s] + r_w*D_p_img_sub_4d[z_b*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_b*D_sz_img_A_sub0 + x_b]);

		long long d = k * D_sz_img_sub1 * D_sz_img_sub0 + row * D_sz_img_sub0 + col;

		D_p_img_sub2tar_4d[d] = u_w*a + d_w*b;



	}

}


__global__ void get_Displacement_affine_new_Y(const int k, const long long gsz1, const long long gsz0, float *D_x4x4_affinematrix, const int D_sz_img_sub0,
	const int D_sz_img_sub1, const int D_sz_img_sub2, const int D_sz_img_sub3, unsigned char *D_p_img_sub_4d, unsigned char *D_p_img_sub2tar_4d, const int D_sz_img_A_sub0,
	const int D_sz_img_A_sub1, const int D_sz_img_A_sub2, const long long start_block_x, const long long start_block_y, const long long start_block_z,
	const long long x_read_offset, const long long y_read_offset, const long long z_read_offset, const long long gs_ori2, const long long gs_ori1, const long long gs_ori0)
{
	const int row = blockIdx.y * blockDim.y + threadIdx.y;
	const int col = blockIdx.x * blockDim.x + threadIdx.x;

	if (row >= gsz1 || col >= gsz0)return;


	float x_pt_sub2tar_homo[4];
	float x_pt_sub_homo[4];


	x_pt_sub2tar_homo[0] = col + start_block_x;
	x_pt_sub2tar_homo[1] = k + start_block_y;
	x_pt_sub2tar_homo[2] = row + start_block_z;
	x_pt_sub2tar_homo[3] = 1.0;

	MUL_aff(D_x4x4_affinematrix, x_pt_sub2tar_homo, x_pt_sub_homo);

	double cur_pos[3];//x,y,z
	double cur_pos_1[3];//x,y,z
	cur_pos[0] = x_pt_sub_homo[0] - x_read_offset;
	cur_pos[1] = x_pt_sub_homo[1] - y_read_offset;
	cur_pos[2] = x_pt_sub_homo[2] - z_read_offset;

	cur_pos_1[0] = x_pt_sub_homo[0];
	cur_pos_1[1] = x_pt_sub_homo[1];
	cur_pos_1[2] = x_pt_sub_homo[2];

	if (cur_pos_1[0]<0 || cur_pos_1[0]>gs_ori0 - 1 ||
		cur_pos_1[1]<0 || cur_pos_1[1]>gs_ori1 - 1 ||
		cur_pos_1[2]<0 || cur_pos_1[2]>gs_ori2 - 1 ||
		cur_pos[0]<0 || cur_pos[0]>D_sz_img_A_sub0 - 1 ||
		cur_pos[1]<0 || cur_pos[1]>D_sz_img_A_sub1 - 1 ||
		cur_pos[2]<0 || cur_pos[2]>D_sz_img_A_sub2 - 1)
	{
		for (long long c = 0; c < D_sz_img_sub3; c++)
		{
			D_p_img_sub2tar_4d[c * D_sz_img_sub2 * D_sz_img_sub1 * D_sz_img_sub0 + row * D_sz_img_sub1 * D_sz_img_sub0 + k * D_sz_img_sub0 + col] = 0.0;

		}

	}
	else{


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


		float a = t_w*(l_w*D_p_img_sub_4d[z_s*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_s*D_sz_img_A_sub0 + x_s] + r_w*D_p_img_sub_4d[z_s*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_s*D_sz_img_A_sub0 + x_b]) + b_w*(l_w*D_p_img_sub_4d[z_s*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_b*D_sz_img_A_sub0 + x_s] + r_w*D_p_img_sub_4d[z_s*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_b*D_sz_img_A_sub0 + x_b]);

		float b = t_w*(l_w*D_p_img_sub_4d[z_b*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_s*D_sz_img_A_sub0 + x_s] + r_w*D_p_img_sub_4d[z_b*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_s*D_sz_img_A_sub0 + x_b]) + b_w*(l_w*D_p_img_sub_4d[z_b*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_b*D_sz_img_A_sub0 + x_s] + r_w*D_p_img_sub_4d[z_b*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_b*D_sz_img_A_sub0 + x_b]);

		long long d = row * D_sz_img_sub1 * D_sz_img_sub0 + k * D_sz_img_sub0 + col;

		D_p_img_sub2tar_4d[d] = u_w*a + d_w*b;



	}

}

__global__ void get_Displacement_affine_new_X(const int k, const long long gsz1, const long long gsz0, float *D_x4x4_affinematrix, const int D_sz_img_sub0,
	const int D_sz_img_sub1, const int D_sz_img_sub2, const int D_sz_img_sub3, unsigned char *D_p_img_sub_4d, unsigned char *D_p_img_sub2tar_4d, const int D_sz_img_A_sub0,
	const int D_sz_img_A_sub1, const int D_sz_img_A_sub2, const long long start_block_x, const long long start_block_y, const long long start_block_z,
	const long long x_read_offset, const long long y_read_offset, const long long z_read_offset, const long long gs_ori2, const long long gs_ori1, const long long gs_ori0)
{
	const int row = blockIdx.y * blockDim.y + threadIdx.y;
	const int col = blockIdx.x * blockDim.x + threadIdx.x;

	if (row >= gsz1 || col >= gsz0)return;


	float x_pt_sub2tar_homo[4];
	float x_pt_sub_homo[4];


	x_pt_sub2tar_homo[0] = k + start_block_x;
	x_pt_sub2tar_homo[1] = col + start_block_y;
	x_pt_sub2tar_homo[2] = row + start_block_z;
	x_pt_sub2tar_homo[3] = 1.0;

	MUL_aff(D_x4x4_affinematrix, x_pt_sub2tar_homo, x_pt_sub_homo);

	double cur_pos[3];//x,y,z
	double cur_pos_1[3];//x,y,z
	cur_pos[0] = x_pt_sub_homo[0] - x_read_offset;
	cur_pos[1] = x_pt_sub_homo[1] - y_read_offset;
	cur_pos[2] = x_pt_sub_homo[2] - z_read_offset;

	cur_pos_1[0] = x_pt_sub_homo[0];
	cur_pos_1[1] = x_pt_sub_homo[1];
	cur_pos_1[2] = x_pt_sub_homo[2];

	if (cur_pos_1[0]<0 || cur_pos_1[0]>gs_ori0 - 1 ||
		cur_pos_1[1]<0 || cur_pos_1[1]>gs_ori1 - 1 ||
		cur_pos_1[2]<0 || cur_pos_1[2]>gs_ori2 - 1 ||
		cur_pos[0]<0 || cur_pos[0]>D_sz_img_A_sub0 - 1 ||
		cur_pos[1]<0 || cur_pos[1]>D_sz_img_A_sub1 - 1 ||
		cur_pos[2]<0 || cur_pos[2]>D_sz_img_A_sub2 - 1)
	{
		for (long long c = 0; c < D_sz_img_sub3; c++)
		{
			D_p_img_sub2tar_4d[c * D_sz_img_sub2 * D_sz_img_sub1 * D_sz_img_sub0 + row * D_sz_img_sub1 * D_sz_img_sub0 + col * D_sz_img_sub0 + k] = 0.0;

		}

	}
	else{


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


		float a = t_w*(l_w*D_p_img_sub_4d[z_s*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_s*D_sz_img_A_sub0 + x_s] + r_w*D_p_img_sub_4d[z_s*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_s*D_sz_img_A_sub0 + x_b]) + b_w*(l_w*D_p_img_sub_4d[z_s*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_b*D_sz_img_A_sub0 + x_s] + r_w*D_p_img_sub_4d[z_s*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_b*D_sz_img_A_sub0 + x_b]);

		float b = t_w*(l_w*D_p_img_sub_4d[z_b*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_s*D_sz_img_A_sub0 + x_s] + r_w*D_p_img_sub_4d[z_b*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_s*D_sz_img_A_sub0 + x_b]) + b_w*(l_w*D_p_img_sub_4d[z_b*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_b*D_sz_img_A_sub0 + x_s] + r_w*D_p_img_sub_4d[z_b*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_b*D_sz_img_A_sub0 + x_b]);

		long long d = row * D_sz_img_sub1 * D_sz_img_sub0 + col * D_sz_img_sub0 + k;

		D_p_img_sub2tar_4d[d] = u_w*a + d_w*b;



	}

}



__global__ void get_Displacement_affine(const int k, const long long gsz1, const long long gsz0, float *D_x4x4_affinematrix, const int D_sz_img_sub0,
	const int D_sz_img_sub1, const int D_sz_img_sub2, const int D_sz_img_sub3, float *D_p_img_sub_4d, float *D_p_img_sub2tar_4d, const int D_sz_img_A_sub0,
	const int D_sz_img_A_sub1, const int D_sz_img_A_sub2)
{
	const int row = blockIdx.y * blockDim.y + threadIdx.y;
	const int col = blockIdx.x * blockDim.x + threadIdx.x;

	if (row >= gsz1 || col >= gsz0)return;

	float x_pt_sub2tar_homo[4];
	float x_pt_sub_homo[4];


	x_pt_sub2tar_homo[0] = col;
	x_pt_sub2tar_homo[1] = row;
	x_pt_sub2tar_homo[2] = k;
	x_pt_sub2tar_homo[3] = 1.0;

	MUL_aff(D_x4x4_affinematrix, x_pt_sub2tar_homo, x_pt_sub_homo);

	double cur_pos[3];//x,y,z
	cur_pos[0] = x_pt_sub_homo[0];
	cur_pos[1] = x_pt_sub_homo[1];
	cur_pos[2] = x_pt_sub_homo[2];

	if (cur_pos[0]<0 || cur_pos[0]>D_sz_img_A_sub0 - 1 ||
		cur_pos[1]<0 || cur_pos[1]>D_sz_img_A_sub1 - 1 ||
		cur_pos[2]<0 || cur_pos[2]>D_sz_img_A_sub2 - 1)
	{
		for (long long c = 0; c < D_sz_img_sub3; c++)
		{
			D_p_img_sub2tar_4d[c * D_sz_img_sub2 * D_sz_img_sub1 * D_sz_img_sub0 + k * D_sz_img_sub1 * D_sz_img_sub0 + row * D_sz_img_sub0 + col] = 0.0;

		}

	}
	else{


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


		float a = t_w*(l_w*D_p_img_sub_4d[z_s*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_s*D_sz_img_A_sub0 + x_s] + r_w*D_p_img_sub_4d[z_s*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_s*D_sz_img_A_sub0 + x_b]) + b_w*(l_w*D_p_img_sub_4d[z_s*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_b*D_sz_img_A_sub0 + x_s] + r_w*D_p_img_sub_4d[z_s*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_b*D_sz_img_A_sub0 + x_b]);

		float b = t_w*(l_w*D_p_img_sub_4d[z_b*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_s*D_sz_img_A_sub0 + x_s] + r_w*D_p_img_sub_4d[z_b*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_s*D_sz_img_A_sub0 + x_b]) + b_w*(l_w*D_p_img_sub_4d[z_b*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_b*D_sz_img_A_sub0 + x_s] + r_w*D_p_img_sub_4d[z_b*D_sz_img_A_sub1 * D_sz_img_A_sub0 + y_b*D_sz_img_A_sub0 + x_b]);

		int d = k * D_sz_img_sub1 * D_sz_img_sub0 + row * D_sz_img_sub0 + col;

		D_p_img_sub2tar_4d[d] = u_w*a + d_w*b;



	}

}

extern "C" bool gpu_interpolation_new(const int mode, const long long gsz2, const long long gsz1, const long long gsz0, const Matrix &x4x4_affinematrix, unsigned char ****&p_img_sub_4d,
	unsigned char ****&p_img_sub2tar_4d, const long long *sz_img_sub, const long long gsA2, const long long gsA1, const long long gsA0, const long long start_block_x, const long long start_block_y,
	const long long start_block_z, const long long x_read_offset, const long long y_read_offset, const long long z_read_offset, const long long gs_ori2, const long long gs_ori1, const long long gs_ori0,
	const unsigned char *p_img_sub, unsigned char *p_img_affine)
{


	float *H_x4x4_affinematrix, *H_sz_img_sub;

	float *D_x4x4_affinematrix, *D_H_sz_img_sub;

	unsigned char *D_p_img_sub_4d, *D_p_img_sub2tar_4d;

	H_x4x4_affinematrix = (float*)malloc(x4x4_affinematrix.nrows() * x4x4_affinematrix.ncols() * sizeof(float));



	cudaMalloc((void**)&D_x4x4_affinematrix, x4x4_affinematrix.nrows() * x4x4_affinematrix.ncols() * sizeof(float));
	cudaMalloc((void**)&D_p_img_sub_4d, gsA0 * gsA1 * gsA2 * sz_img_sub[3] * sizeof(unsigned char));
	cudaMalloc((void**)&D_p_img_sub2tar_4d, sz_img_sub[0] * sz_img_sub[1] * sz_img_sub[2] * sz_img_sub[3] * sizeof(unsigned char));


	int aa = sz_img_sub[0]; int bb = sz_img_sub[1]; int cc = sz_img_sub[2]; int dd = sz_img_sub[3];


	for (int i = 0; i < x4x4_affinematrix.nrows(); i++)
	{
		for (int j = 0; j < x4x4_affinematrix.ncols(); j++)
		{

			H_x4x4_affinematrix[i * x4x4_affinematrix.ncols() + j] = x4x4_affinematrix(i + 1, j + 1);

		}

	}

	HANDLE_ERROR(cudaMemcpy(D_x4x4_affinematrix, H_x4x4_affinematrix, x4x4_affinematrix.nrows() * x4x4_affinematrix.ncols() * sizeof(float), cudaMemcpyHostToDevice));

	HANDLE_ERROR(cudaMemcpy(D_p_img_sub_4d, p_img_sub, gsA0 * gsA1 * gsA2 * sz_img_sub[3] * sizeof(unsigned char), cudaMemcpyHostToDevice));
	/*
	*mode:0; z:gsz2; y:gsz1; x:gsz0; z=k;y=row;x=col;
	*mode:1; y:gsz2; z:gsz1; x:gsz0; y=k;z=row;x=col;
	*mode:2; x:gsz2; z:gsz1; y:gsz0; x=k;z=row;y=col;
	*/
	switch (mode)
	{
	case 0:
		for (long long k = 0; k < gsz2; k++)
		{
			dim3 grid((gsz0 + threads_num - 1) / threads_num, (gsz1 + threads_num - 1) / threads_num);
			dim3 block(threads_num, threads_num);
			//get_Displacement_affine_new << <grid, block >> >(k, bb_real, aa_real, D_x4x4_affinematrix, aa, bb, cc, dd, D_p_img_sub_4d, D_p_img_sub2tar_4d, gsA0, gsA1, gsA2, start_block_x, start_block_y, start_block_z, x_read_offset, y_read_offset, z_read_offset, gs_ori2, gs_ori1, gs_ori0, aa_real, bb_real, cc_real);
			get_Displacement_affine_new_Z << <grid, block >> >(k, gsz1, gsz0, D_x4x4_affinematrix, aa, bb, cc, dd, D_p_img_sub_4d, D_p_img_sub2tar_4d, gsA0, gsA1, gsA2, start_block_x, start_block_y, start_block_z, x_read_offset, y_read_offset, z_read_offset, gs_ori2, gs_ori1, gs_ori0);
			cudaDeviceSynchronize();
		}
		break;
	case 1:
		for (long long k = 0; k < gsz2; k++)
		{
			dim3 grid((gsz0 + threads_num - 1) / threads_num, (gsz1 + threads_num - 1) / threads_num);
			dim3 block(threads_num, threads_num);
			//get_Displacement_affine_new << <grid, block >> >(k, bb_real, aa_real, D_x4x4_affinematrix, aa, bb, cc, dd, D_p_img_sub_4d, D_p_img_sub2tar_4d, gsA0, gsA1, gsA2, start_block_x, start_block_y, start_block_z, x_read_offset, y_read_offset, z_read_offset, gs_ori2, gs_ori1, gs_ori0, aa_real, bb_real, cc_real);
			get_Displacement_affine_new_Y << <grid, block >> >(k, gsz1, gsz0, D_x4x4_affinematrix, aa, bb, cc, dd, D_p_img_sub_4d, D_p_img_sub2tar_4d, gsA0, gsA1, gsA2, start_block_x, start_block_y, start_block_z, x_read_offset, y_read_offset, z_read_offset, gs_ori2, gs_ori1, gs_ori0);
			cudaDeviceSynchronize();
		}
		break;
	case 2:
		for (long long k = 0; k < gsz2; k++)
		{
			dim3 grid((gsz0 + threads_num - 1) / threads_num, (gsz1 + threads_num - 1) / threads_num);
			dim3 block(threads_num, threads_num);
			//get_Displacement_affine_new << <grid, block >> >(k, bb_real, aa_real, D_x4x4_affinematrix, aa, bb, cc, dd, D_p_img_sub_4d, D_p_img_sub2tar_4d, gsA0, gsA1, gsA2, start_block_x, start_block_y, start_block_z, x_read_offset, y_read_offset, z_read_offset, gs_ori2, gs_ori1, gs_ori0, aa_real, bb_real, cc_real);
			get_Displacement_affine_new_X << <grid, block >> >(k, gsz1, gsz0, D_x4x4_affinematrix, aa, bb, cc, dd, D_p_img_sub_4d, D_p_img_sub2tar_4d, gsA0, gsA1, gsA2, start_block_x, start_block_y, start_block_z, x_read_offset, y_read_offset, z_read_offset, gs_ori2, gs_ori1, gs_ori0);
			cudaDeviceSynchronize();
		}
		break;
	default:
		printf("Error:calculate mode is wrong!\n"); return false;

	}

	HANDLE_ERROR(cudaMemcpy(p_img_affine, D_p_img_sub2tar_4d, sz_img_sub[0] * sz_img_sub[1] * sz_img_sub[2] * sz_img_sub[3] * sizeof(unsigned char), cudaMemcpyDeviceToHost));

	free(H_x4x4_affinematrix);
	cudaFree(D_x4x4_affinematrix); cudaFree(D_p_img_sub_4d); cudaFree(D_p_img_sub2tar_4d);
	HANDLE_ERROR(cudaDeviceReset());
	return true;
}


extern "C" bool gpu_interpolation_affine(const long long gsz2, const long long gsz1, const long long gsz0, const Matrix &x4x4_affinematrix, unsigned char ****&p_img_sub_4d, unsigned char ****&p_img_sub2tar_4d, const long long *sz_img_sub,
	const long long gsA2, const long long gsA1, const long long gsA0)
{


	float *H_x4x4_affinematrix, *H_p_img_sub_4d, *H_sz_img_sub, *H_p_img_sub2tar_4d;

	float *D_x4x4_affinematrix, *D_p_img_sub_4d, *D_H_sz_img_sub, *D_p_img_sub2tar_4d;


	H_x4x4_affinematrix = (float*)malloc(x4x4_affinematrix.nrows() * x4x4_affinematrix.ncols() * sizeof(float));
	H_p_img_sub_4d = (float*)malloc(2 * gsA0 * gsA1 * gsA2 * sz_img_sub[3] * sizeof(float));
	H_p_img_sub2tar_4d = (float*)malloc(2 * sz_img_sub[0] * sz_img_sub[1] * sz_img_sub[2] * sz_img_sub[3] * sizeof(float));


	cudaMalloc((void**)&D_x4x4_affinematrix, x4x4_affinematrix.nrows() * x4x4_affinematrix.ncols() * sizeof(float));
	cudaMalloc((void**)&D_p_img_sub_4d, 2 * gsA0 * gsA1 * gsA2 * sz_img_sub[3] * sizeof(float));
	cudaMalloc((void**)&D_p_img_sub2tar_4d, 2 * sz_img_sub[0] * sz_img_sub[1] * sz_img_sub[2] * sz_img_sub[3] * sizeof(float));


	int aa = sz_img_sub[0]; int bb = sz_img_sub[1]; int cc = sz_img_sub[2]; int dd = sz_img_sub[3];

	for (int i = 0; i < x4x4_affinematrix.nrows(); i++)
	{
		for (int j = 0; j < x4x4_affinematrix.ncols(); j++)
		{

			H_x4x4_affinematrix[i * x4x4_affinematrix.ncols() + j] = x4x4_affinematrix(i + 1, j + 1);

		}

	}

	for (long long a = 0; a <gsA2; a++)
	{
		for (long long b = 0; b < gsA1; b++)
		{
			for (long long c = 0; c < gsA0; c++)
			{
				H_p_img_sub_4d[a*gsA1 * gsA0 + b*gsA0 + c] = p_img_sub_4d[0][a][b][c];

			}
		}
	}//ע�⣺H_p_img_sub_4d��H_p_img_sub2tar_4d��ÿ��ά�ȴ�С��ͬ




	cudaMemcpy(D_x4x4_affinematrix, H_x4x4_affinematrix, x4x4_affinematrix.nrows() * x4x4_affinematrix.ncols() * sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpy(D_p_img_sub_4d, H_p_img_sub_4d, 2 * gsA0 * gsA1 * gsA2 * sz_img_sub[3] * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(D_p_img_sub2tar_4d, H_p_img_sub2tar_4d, 2 * sz_img_sub[0] * sz_img_sub[1] * sz_img_sub[2] * sz_img_sub[3] * sizeof(float), cudaMemcpyHostToDevice);





	for (long long k = 0; k < gsz2; k++)
	{
		dim3 grid((gsz0 + threads_num - 1) / threads_num, (gsz1 + threads_num - 1) / threads_num);
		dim3 block(threads_num, threads_num);
		get_Displacement_affine << <grid, block >> >(k, gsz1, gsz0, D_x4x4_affinematrix, aa, bb, cc, dd, D_p_img_sub_4d, D_p_img_sub2tar_4d, gsA0, gsA1, gsA2);
		cudaDeviceSynchronize();
	}


	cudaMemcpy(H_p_img_sub2tar_4d, D_p_img_sub2tar_4d, 2 * sz_img_sub[0] * sz_img_sub[1] * sz_img_sub[2] * sz_img_sub[3] * sizeof(float), cudaMemcpyDeviceToHost);

	for (long long a = 0; a < sz_img_sub[2]; a++)
	{
		for (long long b = 0; b < sz_img_sub[1]; b++)
		{
			for (long long c = 0; c < sz_img_sub[0]; c++)
			{
				p_img_sub2tar_4d[0][a][b][c] = H_p_img_sub2tar_4d[a*sz_img_sub[1] * sz_img_sub[0] + b*sz_img_sub[0] + c];
			}
		}
	}

	cudaFree(D_x4x4_affinematrix); cudaFree(D_p_img_sub_4d); cudaFree(D_p_img_sub2tar_4d);
	return true;
}