#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>
#include <malloc.h>
#include<time.h>
#include"newmat.h"
#include "ele.cuh"
#include "cuda_runtime_api.h"
#define nCpt_Max 20000

#if defined(_MSC_VER) && (_WIN64)
//#if defined(_MSC_VER) && defined(_WIN64) //correct?

#define V3DLONG long long

#else

#define V3DLONG long

#endif

typedef unsigned char UINT8_JBA;
typedef float MYFLOAT_JBA;

class DisplaceFieldF3D
{
public:
	UINT8_JBA b_transform;

	MYFLOAT_JBA sx, sy, sz; //shift of x,y,z
	DisplaceFieldF3D() { sx = sy = sz = 0; b_transform = 0; }
	DisplaceFieldF3D(double vv) { sx = sy = sz = vv; b_transform = 0; }
	void scale(double dfactor) { sx *= dfactor; sy *= dfactor; sz *= dfactor; }
	void resetToDefault() //070517
	{
		sx = 0; sy = 0; sz = 0;
		b_transform = 0;
	}
	bool copy(DisplaceFieldF3D *wp)
	{
		if (!wp) return false;
		sx = wp->sx;  sy = wp->sy;  sz = wp->sz;
		b_transform = wp->b_transform;
		return true;
	}
	bool copy(DisplaceFieldF3D &wp)
	{
		sx = wp.sx;  sy = wp.sy;  sz = wp.sz;
		b_transform = wp.b_transform;
		return true;
	}

};


__device__ void MUL(float *A, float *B, float *C, float *D, float E[4], int nCpt)//MUL(x_ori, D_x4x4_d, xmxn_K1, D_xnx4_c, x_stps, nCpt)
{
	int k, a; int p1 = 0; int p2 = 0;

	for (int i = 0; i < 4; i++)
	{
		float sum1 = 0; float sum2 = 0;
		for (k = 0; k < 4; k++)
		{
			sum1 += A[k] * B[k + p1];
		}
		for (a = 0; a < nCpt; a++)
		{
			sum2 += C[a] * D[a + p2];
		}
		E[i] = sum1 + sum2;
		p1 = p1 + 4;
		p2 = p2 + nCpt;
	}


}
__device__ void assignment(V3DLONG k, int row, int col, V3DLONG gsz1, V3DLONG gsz0, V3DLONG gfactor_x, V3DLONG gfactor_y, V3DLONG gfactor_z, float *D_RESULT_X, float *D_RESULT_Y, float *D_RESULT_Z, const float x_stps[4])
{
	//printf("\t>>gfactor_z[%d] [%d] : %f  %f  %f  %f  k=%d:\n", row, col, x_stps[0], x_stps[1], x_stps[2], x_stps[3], k);
	int id = k * 88 + row*gsz0 + col;
	float a = x_stps[1], b = x_stps[2], c = x_stps[3]; __syncthreads();
	//printf("\t>>gfactor_z[%d] [%d] : %f  %f  %f  %f  k=%d:\n", row, col, x_stps[0], a, b, c, k);
	D_RESULT_Y[id] = b - (row - 1)*gfactor_y;// printf("\t>>gfactor_z[%d] [%d] : %f  %f  %f  %f  k=%d:\n", row, col, x_stps[0], x_stps[1], x_stps[2], x_stps[3], k);

	D_RESULT_X[id] = a - (col - 1)*gfactor_x; //__syncthreads();

	D_RESULT_Z[id] = c - (k - 1)*gfactor_z;
	//	printf("\t>>gfactor_z[%d] [%d] :%d  %f  %f  %f  k=%d:\n", row, col, id,  D_RESULT_X[id], D_RESULT_Y[id], D_RESULT_Z[id], k);
}
__global__ void get_cd(int nCpt, V3DLONG k, const V3DLONG gsz1, const V3DLONG gsz0, V3DLONG gfactor_x, V3DLONG gfactor_y, V3DLONG gfactor_z,
	float * D_X, float * D_Y, float * D_Z, float *D_x4x4_d, float *D_xnx4_c, float *D_RESULT_X, float *D_RESULT_Y, float *D_RESULT_Z,
	long long x_offset, long long y_offset, long long z_offset)
{
	const int row = blockIdx.y * blockDim.y + threadIdx.y;
	const int col = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= gsz1 || col >= gsz0)return;
	float x_ori[4]; float xmxn_K1[nCpt_Max]; float x_stps[4];
	x_stps[0] = 0; x_stps[1] = 0; x_stps[2] = 0; x_stps[3] = 0;
	x_ori[0] = 1.0; x_ori[1] = (col - 1)*gfactor_x + x_offset; x_ori[2] = (row - 1)*gfactor_y + y_offset; x_ori[3] = (k - 1)*gfactor_z + z_offset;
	for (int n = 0; n < nCpt; n++)
	{
		xmxn_K1[n] = -sqrt(((col - 1)*gfactor_x + x_offset - D_X[n])*((col - 1)*gfactor_x + x_offset - D_X[n]) + ((row - 1)*gfactor_y + y_offset - D_Y[n])*((row - 1)*gfactor_y + y_offset - D_Y[n]) + ((k - 1)*gfactor_z + z_offset - D_Z[n])*((k - 1)*gfactor_z + z_offset - D_Z[n]));
	}


	MUL(x_ori, D_x4x4_d, xmxn_K1, D_xnx4_c, x_stps, nCpt);


	D_RESULT_X[k*gsz1*gsz0 + row*gsz0 + col] = x_stps[1] - ((col - 1)*gfactor_x + x_offset);
	D_RESULT_Y[k*gsz1*gsz0 + row*gsz0 + col] = x_stps[2] - ((row - 1)*gfactor_y + y_offset);
	D_RESULT_Z[k*gsz1*gsz0 + row*gsz0 + col] = x_stps[3] - ((k - 1)*gfactor_z + z_offset);

}



extern "C" bool gpu_computedistance(int nCpt, const V3DLONG gsz2, const V3DLONG gsz1, const V3DLONG gsz0, V3DLONG gfactor_x, V3DLONG gfactor_y, V3DLONG gfactor_z,
	Matrix &x4x4_d, Matrix &xnx4_c, float * H_X, float * H_Y, float * H_Z, DisplaceFieldF3D *** df_local_3d, long long x_offset, long long y_offset, long long z_offset)
{
	float *D_X, *D_Y, *D_Z, *H_xmxn_K, *D_xmxn_K, *H_ori, *D_ori, *H_stps, *D_stps, *H_x4x4_d, *D_x4x4_d, *H_xnx4_c, *D_xnx4_c, *H_RESULT_X, *H_RESULT_Y, *H_RESULT_Z, *D_RESULT_X, *D_RESULT_Y, *D_RESULT_Z;
	H_xmxn_K = (float*)malloc(nCpt * sizeof(float));
	H_ori = (float*)malloc(4 * sizeof(float));
	H_stps = (float*)malloc(4 * sizeof(float));
	H_x4x4_d = (float*)malloc(4 * 4 * sizeof(float));
	H_xnx4_c = (float*)malloc(nCpt * 4 * sizeof(float));
	H_RESULT_X = (float*)malloc(gsz2 * gsz1 *gsz0* sizeof(float));
	H_RESULT_Y = (float*)malloc(gsz2 * gsz1 *gsz0* sizeof(float));
	H_RESULT_Z = (float*)malloc(gsz2 * gsz1 *gsz0* sizeof(float));

	cudaMalloc((void**)&D_X, nCpt * sizeof(float));
	cudaMalloc((void**)&D_Y, nCpt * sizeof(float));
	cudaMalloc((void**)&D_Z, nCpt * sizeof(float));
	cudaMalloc((void**)&D_xmxn_K, nCpt * sizeof(float));
	cudaMalloc((void**)&D_x4x4_d, 4 * 4 * sizeof(float));
	cudaMalloc((void**)&D_xnx4_c, nCpt * 4 * sizeof(float));
	cudaMalloc((void**)&D_ori, 4 * sizeof(float));
	cudaMalloc((void**)&D_stps, 4 * sizeof(float));
	cudaMalloc((void**)&D_RESULT_X, gsz2 * gsz1 *gsz0* sizeof(float));
	cudaMalloc((void**)&D_RESULT_Y, gsz2 * gsz1 *gsz0* sizeof(float));
	cudaMalloc((void**)&D_RESULT_Z, gsz2 * gsz1 *gsz0* sizeof(float));
	printf("\t>>gsz2 : %d \n", gsz2);
	printf("\t>>gsz1 : %d \n", gsz1);
	printf("\t>>gsz0 : %d \n", gsz0);
	for (int i = 0; i < x4x4_d.nrows(); i++)//x4x4_d����ת��������H_X4X4_d
	{
		for (int j = 0; j < x4x4_d.ncols(); j++)
		{

			H_x4x4_d[j * x4x4_d.nrows() + i] = x4x4_d(i + 1, j + 1);

		}
	}



	for (int i = 0; i < xnx4_c.nrows(); i++)//xnx4_c����ת��������H_xnx4_c
	{
		for (int j = 0; j < xnx4_c.ncols(); j++)
		{

			H_xnx4_c[j * xnx4_c.nrows() + i] = xnx4_c(i + 1, j + 1);

		}
	}


	HANDLE_ERROR(cudaMemcpy(D_X, H_X, nCpt * sizeof(float), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(D_Y, H_Y, nCpt * sizeof(float), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(D_Z, H_Z, nCpt * sizeof(float), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(D_x4x4_d, H_x4x4_d, 4 * 4 * sizeof(float), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(D_xnx4_c, H_xnx4_c, 4 * nCpt* sizeof(float), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(D_RESULT_X, H_RESULT_X, gsz2 * gsz1 *gsz0* sizeof(float), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(D_RESULT_Y, H_RESULT_Y, gsz2 * gsz1 *gsz0* sizeof(float), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(D_RESULT_Z, H_RESULT_Z, gsz2 * gsz1 *gsz0* sizeof(float), cudaMemcpyHostToDevice));

	for (int k = 0; k < gsz2; k++)
	{
		dim3 grid((gsz0 + threads_num - 1) / threads_num, (gsz1 + threads_num - 1) / threads_num);//dim3 ����һ����ά�������������һάΪ1
		dim3 block(threads_num, threads_num);
		get_cd << <grid, block >> >(nCpt, k, gsz1, gsz0, gfactor_x, gfactor_y, gfactor_z, D_X, D_Y, D_Z, D_x4x4_d, D_xnx4_c, D_RESULT_X, D_RESULT_Y, D_RESULT_Z, x_offset, y_offset, z_offset);
		cudaDeviceSynchronize();
		//printf("---------mark[%d]-------\n\n", k);
	}

	HANDLE_ERROR(cudaMemcpy(H_RESULT_X, D_RESULT_X, gsz2 * gsz1 * gsz0 * sizeof(float), cudaMemcpyDeviceToHost));//�豸������ ��������
	HANDLE_ERROR(cudaMemcpy(H_RESULT_Y, D_RESULT_Y, gsz2 * gsz1 * gsz0 * sizeof(float), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(H_RESULT_Z, D_RESULT_Z, gsz2 * gsz1 * gsz0 * sizeof(float), cudaMemcpyDeviceToHost));



	for (V3DLONG a = 0; a < gsz2; a++)
	{
		for (V3DLONG b = 0; b < gsz1; b++)
		{
			for (V3DLONG c = 0; c < gsz0; c++)
			{
				df_local_3d[a][b][c].sx = H_RESULT_X[a*gsz1*gsz0 + b*gsz0 + c]; //printf("\t>>gfactor_z[%d] [%d] :x=%.3f\n ", b, c, df_local_3d[a][b][c].sx);
				df_local_3d[a][b][c].sy = H_RESULT_Y[a*gsz1*gsz0 + b*gsz0 + c]; //printf("y=%.3f ", df_local_3d[a][b][c].sy);
				df_local_3d[a][b][c].sz = H_RESULT_Z[a*gsz1*gsz0 + b*gsz0 + c];// printf("z=%.3f  ", df_local_3d[a][b][c].sz);
			}
		}
	}

	free(H_xmxn_K);
	free(H_ori);
	free(H_stps);
	free(H_xnx4_c);
	free(H_x4x4_d);
	free(H_RESULT_X);
	free(H_RESULT_Y);
	free(H_RESULT_Z);

	cudaFree(&D_X); cudaFree(&D_ori); cudaFree(&D_stps); cudaFree(&D_x4x4_d); cudaFree(&D_xnx4_c);
	cudaFree(&D_Y); cudaFree(&D_RESULT_X); cudaFree(&D_RESULT_Y); cudaFree(&D_RESULT_Z);
	cudaFree(&D_Z); cudaFree(&D_xmxn_K);
	HANDLE_ERROR(cudaDeviceReset());
	return true;
}