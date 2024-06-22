//imgwarp_tps_quicksmallmemory.cpp
//
// by Lei Qu
//2012-07-08

//#include "basic_memory.cpp"
#include "q_imgwarp_tps_quicksmallmemory.h"
#include "Bigwarp.h"
#include<time.h>
#include "RawFmtMngr.h"
#include "omp.h"
#if defined(_MSC_VER) && (_WIN64)
//#if defined(_MSC_VER) && defined(_WIN64) //correct?

#define V3DLONG long long

#else

#define V3DLONG long

#endif

extern "C" int gpu_QR(int ncpt, const Matrix &A, Matrix &Q, Matrix &R);
extern "C" int gpu_extendornormal(int ncpt, int n, Matrix &Q);
extern "C" bool gpu_A(int ncpt, const Matrix &q2_t, const Matrix &xnxn_K, Matrix &A);
extern "C" int gpu_A_i_new(int ncpt, const Matrix &A, Matrix &A_i);
extern "C" bool gpu_xnxn(int ncpt, const Matrix &q2, const Matrix &A_t, Matrix &C);
extern "C" bool gpu_computedistance(int nCpt, const V3DLONG gsz2, const V3DLONG gsz1, const V3DLONG gsz0, V3DLONG gfactor_x, V3DLONG gfactor_y, V3DLONG gfactor_z,
	Matrix &x4x4_d, Matrix &xnx4_c, float * H_X, float * H_Y, float * H_Z, DisplaceFieldF3D *** df_local_3d, long long x_offset, long long y_offset, long long z_offset);
extern "C" bool gpu_interpolation(const int gsz2, const int gsz1, const int gsz0, DisplaceFieldF3D ***&pppSubDF, const Matrix &x_bsplinebasis, const V3DLONG sz_gridwnd, DisplaceFieldF3D ***&pppDFBlock,
	unsigned char ****&p_img_sub_4d, const V3DLONG *sz_img_sub, const V3DLONG szBlock_x, const V3DLONG szBlock_y, const V3DLONG szBlock_z, const int i_interpmethod_img, unsigned char ****&p_img_warp_4d,
	const V3DLONG *sz_img_sub_read, const unsigned char *p_img_sub, unsigned char *p_img_warp, const long long start_block_x, const long long start_block_y, const long long start_block_z,
	const long long x_read_offset, const long long y_read_offset, const long long z_read_offset, const long long gs_ori2, const long long gs_ori1, const long long gs_ori0);
extern "C" bool gpu_interpolation_sort(const int gsz2, const int gsz1, const int gsz0, DisplaceFieldF3D ***&pppSubDF, const Matrix &x_bsplinebasis, const V3DLONG sz_gridwnd, DisplaceFieldF3D ***&pppDFBlock,
	const V3DLONG *sz_img_sub, const V3DLONG szBlock_x, const V3DLONG szBlock_y, const V3DLONG szBlock_z, const int i_interpmethod_img, const long long start_block_x, const long long start_block_y,
	const long long start_block_z, const long long gs_ori2, const long long gs_ori1, const long long gs_ori0, float * &sort_x, float * &sort_y, float * &sort_z);
extern "C" Matrix matrixMultiply(const int m, const int n, const int k, Matrix &A, Matrix &B);
extern "C" bool gpu_interpolation_stps(const int gsz2, const int gsz1, const int gsz0, DisplaceFieldF3D ***&pppSubDF, const Matrix &x_bsplinebasis, const V3DLONG sz_gridwnd, DisplaceFieldF3D ***&pppDFBlock,
                                       unsigned char ****&p_img_sub_4d, const V3DLONG *sz_img_sub, const V3DLONG szBlock_x, const V3DLONG szBlock_y, const V3DLONG szBlock_z, const int i_interpmethod_img, unsigned char ****&p_img_warp_4d,
                                       const V3DLONG *sz_img_sub_ori, const unsigned char *p_img_sub, unsigned char *p_img_warp);


bool STPS_interpolate_CYF(const int GPU_type, Vol3DSimple<DisplaceFieldF3D> *pSubDF, DisplaceFieldF3D ***pppSubDF, Matrix x_bsplinebasis,
                          const V3DLONG sz_gridwnd, DisplaceFieldF3D ***&pppDFBlock, unsigned char ****&p_img_sub_4d, const long long *sz_img_sub_, const V3DLONG szBlock_x,
                          const V3DLONG szBlock_y, const V3DLONG szBlock_z, const int i_interpmethod_img, unsigned char ****&p_img_warp_4d, const long long *sz_img_sub_ori_,
                          const unsigned char *p_img_sub, unsigned char *p_img_warp, long long x_offset, long long y_offset, long long z_offset, const long long *sz_img_ori,
                          const long long *sz_img_read_block)
{
    V3DLONG sz_img_sub[4] = {sz_img_sub_[0], sz_img_sub_[1], sz_img_sub_[2], sz_img_sub_[3]};
    V3DLONG sz_img_sub_ori[4] = {sz_img_sub_ori_[0], sz_img_sub_ori_[1], sz_img_sub_ori_[2], sz_img_sub_ori_[3]};
    if (GPU_type == 1)
    {
        int gsz2 = pSubDF->sz2();
        int gsz1 = pSubDF->sz1();
        int gsz0 = pSubDF->sz0();
        clock_t stps_interpolation;
        stps_interpolation = clock();


        gpu_interpolation(gsz2, gsz1, gsz0, pppSubDF, x_bsplinebasis, sz_gridwnd, pppDFBlock, p_img_sub_4d, sz_img_sub, szBlock_x, szBlock_y, szBlock_z, i_interpmethod_img,
                          p_img_warp_4d, sz_img_sub_ori, p_img_sub, p_img_warp, x_offset, y_offset, z_offset, sz_img_read_block[0], sz_img_read_block[1], sz_img_read_block[2],sz_img_ori[2],
                          sz_img_ori[1],sz_img_ori[0]);


        printf("\t>>interpolation time consume %.2f s\n", (float)(clock() - stps_interpolation) / CLOCKS_PER_SEC);
    }
    else
    {
        for (V3DLONG substart_z = 0; substart_z<pSubDF->sz2() - 1 - 2; substart_z++)
            for (V3DLONG substart_y = 0; substart_y<pSubDF->sz1() - 1 - 2; substart_y++)
                for (V3DLONG substart_x = 0; substart_x<pSubDF->sz0() - 1 - 2; substart_x++)
                {
                    //bspline interpolate the SubDfBlock to DFBlock
                    q_dfblcokinterp_bspline(pppSubDF, x_bsplinebasis, sz_gridwnd, substart_x, substart_y, substart_z, pppDFBlock);
                    //warp image block using DFBlock
                    q_imgblockwarp_CYF(p_img_sub_4d, sz_img_sub, pppDFBlock, szBlock_x, szBlock_y, szBlock_z, i_interpmethod_img, substart_x, substart_y, substart_z, p_img_warp_4d,
                                       sz_img_sub_ori, x_offset, y_offset, z_offset, sz_img_ori, sz_img_read_block);
                }
    }
    return true;
}



bool q_imgblockwarp(unsigned char ****&p_img_sub_4d, const V3DLONG *sz_img_sub, DisplaceFieldF3D ***&pppDFBlock,
                    const V3DLONG szBlock_x,const V3DLONG szBlock_y,const V3DLONG szBlock_z,const int i_interpmethod_img,
                    const V3DLONG substart_x,const V3DLONG substart_y,const V3DLONG substart_z,
                    unsigned char ****&p_img_warp_4d, const V3DLONG *sz_img_sub_ori)
{
    V3DLONG start_x,start_y,start_z;
    start_x=substart_x*szBlock_x;
    start_y=substart_y*szBlock_y;
    start_z=substart_z*szBlock_z;
    for(V3DLONG z=0;z<szBlock_z;z++)
        for(V3DLONG y=0;y<szBlock_y;y++)
            for(V3DLONG x=0;x<szBlock_x;x++)
            {
                V3DLONG pos_warp[3];
                pos_warp[0]=start_x+x;
                pos_warp[1]=start_y+y;
                pos_warp[2]=start_z+z;
                if(pos_warp[0]>=sz_img_sub[0] || pos_warp[1]>=sz_img_sub[1] || pos_warp[2]>=sz_img_sub[2])
                    continue;

                double pos_sub[3];
                pos_sub[0]=pos_warp[0]+pppDFBlock[z][y][x].sx;
                pos_sub[1]=pos_warp[1]+pppDFBlock[z][y][x].sy;
                pos_sub[2]=pos_warp[2]+pppDFBlock[z][y][x].sz;
                if (pos_sub[0]<0 || pos_sub[0]>sz_img_sub_ori[0] - 1 ||
                    pos_sub[1]<0 || pos_sub[1]>sz_img_sub_ori[1] - 1 ||
                    pos_sub[2]<0 || pos_sub[2]>sz_img_sub_ori[2] - 1)
                {
                    for(V3DLONG c=0;c<sz_img_sub[3];c++)
                        p_img_warp_4d[c][pos_warp[2]][pos_warp[1]][pos_warp[0]]=0;
                    continue;
                }

                //nearest neighbor interpolate
                if(i_interpmethod_img==1)
                {
                    V3DLONG pos_sub_nn[3];
                    for(int i=0;i<3;i++)
                    {
                        pos_sub_nn[i]=pos_sub[i]+0.5;
                        pos_sub_nn[i]=pos_sub_nn[i]<sz_img_sub[i]?pos_sub_nn[i]:sz_img_sub[i]-1;
                    }
                    for(V3DLONG c=0;c<sz_img_sub[3];c++)
                        p_img_warp_4d[c][pos_warp[2]][pos_warp[1]][pos_warp[0]]=p_img_sub_4d[c][pos_sub_nn[2]][pos_sub_nn[1]][pos_sub_nn[0]];
                }
                    //linear interpolate
                else if(i_interpmethod_img==0)
                {
                    //find 8 neighor pixels boundary
                    V3DLONG x_s,x_b,y_s,y_b,z_s,z_b;
                    x_s=floor(pos_sub[0]);		x_b=ceil(pos_sub[0]);
                    y_s=floor(pos_sub[1]);		y_b=ceil(pos_sub[1]);
                    z_s=floor(pos_sub[2]);		z_b=ceil(pos_sub[2]);

                    //compute weight for left and right, top and bottom -- 4 neighbor pixel's weight in a slice
                    double l_w,r_w,t_w,b_w;
                    l_w=1.0-(pos_sub[0]-x_s);	r_w=1.0-l_w;
                    t_w=1.0-(pos_sub[1]-y_s);	b_w=1.0-t_w;
                    //compute weight for higer slice and lower slice
                    double u_w,d_w;
                    u_w=1.0-(pos_sub[2]-z_s);	d_w=1.0-u_w;

                    //linear interpolate each channel
                    for(V3DLONG c=0;c<sz_img_sub[3];c++)
                    {
                        //linear interpolate in higher slice [t_w*(l_w*lt+r_w*rt)+b_w*(l_w*lb+r_w*rb)]
                        double higher_slice;
                        higher_slice=t_w*(l_w*p_img_sub_4d[c][z_s][y_s][x_s]+r_w*p_img_sub_4d[c][z_s][y_s][x_b])+
                                     b_w*(l_w*p_img_sub_4d[c][z_s][y_b][x_s]+r_w*p_img_sub_4d[c][z_s][y_b][x_b]);
                        //linear interpolate in lower slice [t_w*(l_w*lt+r_w*rt)+b_w*(l_w*lb+r_w*rb)]
                        double lower_slice;
                        lower_slice =t_w*(l_w*p_img_sub_4d[c][z_b][y_s][x_s]+r_w*p_img_sub_4d[c][z_b][y_s][x_b])+
                                     b_w*(l_w*p_img_sub_4d[c][z_b][y_b][x_s]+r_w*p_img_sub_4d[c][z_b][y_b][x_b]);
                        //linear interpolate the current position [u_w*higher_slice+d_w*lower_slice]
                        double intval = (u_w*higher_slice + d_w*lower_slice + 0.5);
                        p_img_warp_4d[c][pos_warp[2]][pos_warp[1]][pos_warp[0]]=intval;
                    }
                }

            }

    return true;
}

bool q_dfblcokinterp_bspline(DisplaceFieldF3D ***&pppSubDF, const Matrix &x_bsplinebasis,
	const V3DLONG sz_gridwnd, const V3DLONG substart_x, const V3DLONG substart_y, const V3DLONG substart_z,
	DisplaceFieldF3D ***&pppDFBlock)
{
	//vectorize the gridblock's nodes position that use for interpolation
	Matrix x1D_gridblock(4 * 4 * 4, 3);
	long ind = 1;
	for (long dep = substart_z; dep<substart_z + 4; dep++)
		for (long col = substart_x; col<substart_x + 4; col++)
			for (long row = substart_y; row<substart_y + 4; row++)
			{
				x1D_gridblock(ind, 1) = pppSubDF[dep][row][col].sx;
				x1D_gridblock(ind, 2) = pppSubDF[dep][row][col].sy;
				x1D_gridblock(ind, 3) = pppSubDF[dep][row][col].sz;
				ind++;
			}
	//printf("grid[%d %d],%f %f %f\n", substart_y, substart_x, x1D_gridblock(1, 3), x1D_gridblock(1, 2), x1D_gridblock(1, 1));

	//cubic B-spline interpolate the vectorized grid block
	Matrix x1D_gridblock_int = x_bsplinebasis*x1D_gridblock;
	//printf("grid[%d %d],%f %f %f\n", substart_y, substart_x, x1D_gridblock_int(1, 3), x1D_gridblock_int(1, 2), x1D_gridblock_int(1, 1));
	//de-vectorize the interpolated grid block and save back to vec4D_grid_int
	ind = 1;
	for (long zz = 0; zz<sz_gridwnd; zz++)
		for (long xx = 0; xx<sz_gridwnd; xx++)
			for (long yy = 0; yy<sz_gridwnd; yy++)
			{
				pppDFBlock[zz][yy][xx].sx = x1D_gridblock_int(ind, 1);
				pppDFBlock[zz][yy][xx].sy = x1D_gridblock_int(ind, 2);
				pppDFBlock[zz][yy][xx].sz = x1D_gridblock_int(ind, 3);
				ind++;
			}

	return true;
}

Vol3DSimple<DisplaceFieldF3D> * compute_df_tps_subsampled_volume(const vector <Coord3D_JBA> & matchTargetPos, const vector <Coord3D_JBA> & matchSubjectPos, V3DLONG sz0, V3DLONG sz1, V3DLONG sz2,
	V3DLONG gfactor_x, V3DLONG gfactor_y, V3DLONG gfactor_z)
{
	int nCpt = matchTargetPos.size();
	if (nCpt != matchSubjectPos.size() || nCpt <= 0)
	{
		fprintf(stderr, "The input vectors are invalid in compute_tps_df_field().\n");
		return 0;
	}

	Image2DSimple<MYFLOAT_JBA> * cpt_target = new Image2DSimple<MYFLOAT_JBA>(3, nCpt);
	Image2DSimple<MYFLOAT_JBA> * cpt_subject = new Image2DSimple<MYFLOAT_JBA>(3, nCpt);
	if (!cpt_target || !cpt_target->valid() || !cpt_subject || !cpt_subject->valid())
	{
		fprintf(stderr, "Fail to allocate memory.");
		if (cpt_target) { delete cpt_target; cpt_target = 0; }
		if (cpt_subject) { delete cpt_subject; cpt_subject = 0; }
		return 0;
	}

	V3DLONG n;

	MYFLOAT_JBA ** cpt_target_ref = cpt_target->getData2dHandle();
	MYFLOAT_JBA ** cpt_subject_ref = cpt_subject->getData2dHandle();
	//printf("\n---------------------------------\n");
	for (n = 0; n<nCpt; n++)
	{
		cpt_target_ref[n][0] = matchTargetPos.at(n).x;
		cpt_target_ref[n][1] = matchTargetPos.at(n).y;
		cpt_target_ref[n][2] = matchTargetPos.at(n).z;

		cpt_subject_ref[n][0] = matchSubjectPos.at(n).x;
		cpt_subject_ref[n][1] = matchSubjectPos.at(n).y;
		cpt_subject_ref[n][2] = matchSubjectPos.at(n).z;

		//printf("n=%d \tx=[%5.3f -> %5.3f] y=[%5.3f -> %5.3f] z=[%5.3f -> %5.3f] \n",
		//       n, cpt_target_ref[n][0], cpt_subject_ref[n][0], cpt_target_ref[n][1], cpt_subject_ref[n][1], cpt_target_ref[n][2], cpt_subject_ref[n][2]);
	}
	//printf("\n#################################\n");

	Matrix wR(nCpt, nCpt);

	double tmp, s;

	V3DLONG i, j, k;
	for (j = 0; j<nCpt; j++)
	{
		for (i = 0; i<nCpt; i++)
		{
			s = 0.0;
			tmp = cpt_target_ref[i][0] - cpt_target_ref[j][0]; s += tmp*tmp;
			tmp = cpt_target_ref[i][1] - cpt_target_ref[j][1]; s += tmp*tmp;
			tmp = cpt_target_ref[i][2] - cpt_target_ref[j][2]; s += tmp*tmp;
			wR(i + 1, j + 1) = 2 * s*log(s + 1e-20);
		}
	}

	Matrix wP(nCpt, 4);
	for (j = 0; j<nCpt; j++)
	{
		wP(j + 1, 1) = 1;
		wP(j + 1, 2) = cpt_target_ref[j][0];
		wP(j + 1, 3) = cpt_target_ref[j][1];
		wP(j + 1, 4) = cpt_target_ref[j][2];
	}

	Matrix wL(nCpt + 4, nCpt + 4);
	wL.submatrix(1, nCpt, 1, nCpt) = wR;
	wL.submatrix(1, nCpt, nCpt + 1, nCpt + 4) = wP;
	wL.submatrix(nCpt + 1, nCpt + 4, 1, nCpt) = wP.t();
	wL.submatrix(nCpt + 1, nCpt + 4, nCpt + 1, nCpt + 4) = 0;

	Matrix wY(nCpt + 4, 3);
	for (j = 0; j<nCpt; j++)
	{
		wY(j + 1, 1) = cpt_subject_ref[j][0];
		wY(j + 1, 2) = cpt_subject_ref[j][1];
		wY(j + 1, 3) = cpt_subject_ref[j][2];
	}
	wY.submatrix(nCpt + 1, nCpt + 4, 1, 3) = 0;

	Matrix wW;

	Try
	{
		wW = wL.i() * wY;
	}
		CatchAll
	{
		fprintf(stderr, "Fail to find the inverse of the wL matrix.\n");

		if (cpt_target) { delete cpt_target; cpt_target = 0; }
		if (cpt_subject) { delete cpt_subject; cpt_subject = 0; }
		return 0;
	}

	V3DLONG p;

	V3DLONG gsz0 = (V3DLONG)(ceil((double(sz0) / gfactor_x))) + 1, gsz1 = (V3DLONG)(ceil((double(sz1) / gfactor_y))) + 1, gsz2 = (V3DLONG)(ceil((double(sz2) / gfactor_z))) + 1;
	Vol3DSimple<DisplaceFieldF3D> * df_local = new Vol3DSimple<DisplaceFieldF3D>(gsz0, gsz1, gsz2);
	DisplaceFieldF3D *** df_local_3d = df_local->getData3dHandle();

	if (!df_local || !df_local->valid())
	{
		fprintf(stderr, "Fail to allocate memory for the subsampled DF volume memory [%d].\n", __LINE__);

		if (cpt_target) { delete cpt_target; cpt_target = 0; }
		if (cpt_subject) { delete cpt_subject; cpt_subject = 0; }
		if (df_local) { delete df_local; df_local = 0; }
		return 0;
	}

	V3DLONG ndimpt = 3;
	double * dist = new double[nCpt + ndimpt + 1];
	if (!dist)
	{
		fprintf(stderr, "Fail to allocate memory dist for tps warping [%d].\n", __LINE__);

		if (cpt_target) { delete cpt_target; cpt_target = 0; }
		if (cpt_subject) { delete cpt_subject; cpt_subject = 0; }
		if (df_local) { delete df_local; df_local = 0; }
		return 0;
	}

	printf("-------------------- Now compute the distances of pixels to the mapping points. -------\n\n");

	DisplaceFieldF3D * df_local_1d = df_local->getData1dHandle();
	for (k = 0; k<df_local->getTotalElementNumber(); k++)
	{
		df_local_1d[k].sz = df_local_1d[k].sy = df_local_1d[k].sx = 0;
	}
	for (k = 0; k<gsz2; k++)
	{
		for (j = 0; j<gsz1; j++)
		{
			for (i = 0; i<gsz0; i++)
			{
				for (n = 0; n<nCpt; n++)
				{
					s = 0;
					tmp = (i*gfactor_x) - cpt_target_ref[n][0]; s += tmp*tmp;
					tmp = (j*gfactor_y) - cpt_target_ref[n][1]; s += tmp*tmp;
					tmp = (k*gfactor_z) - cpt_target_ref[n][2]; s += tmp*tmp;
					dist[n] = 2 * s*log(s + 1e-20);
				}

				dist[nCpt] = 1;
				dist[nCpt + 1] = i*gfactor_x;
				dist[nCpt + 2] = j*gfactor_y;
				dist[nCpt + 3] = k*gfactor_z;

				s = 0;  for (p = 0; p<nCpt + ndimpt + 1; p++) { s += dist[p] * wW(p + 1, 1); }
				df_local_3d[k][j][i].sx = s - i*gfactor_x;

				s = 0;  for (p = 0; p<nCpt + ndimpt + 1; p++) { s += dist[p] * wW(p + 1, 2); }
				df_local_3d[k][j][i].sy = s - j*gfactor_y;

				s = 0;  for (p = 0; p<nCpt + ndimpt + 1; p++) { s += dist[p] * wW(p + 1, 3); }
				df_local_3d[k][j][i].sz = s - k*gfactor_z;
			}//i
		}//j
		printf("z=%ld ", k); fflush(stdout);
	}//k
	printf("\n");

	if (dist) { delete[]dist; dist = 0; }
	if (cpt_target) { delete cpt_target; cpt_target = 0; }
	if (cpt_subject) { delete cpt_subject; cpt_subject = 0; }

	return df_local;
}

Vol3DSimple<DisplaceFieldF3D> * compute_df_stps_subsampled_volume_4bspline(const vector <Coord3D_JBA> & matchTargetPos, const vector <Coord3D_JBA> & matchSubjectPos, V3DLONG sz0, V3DLONG sz1, V3DLONG sz2,
	V3DLONG gfactor_x, V3DLONG gfactor_y, V3DLONG gfactor_z, int gpu_type)
{
	int nCpt = matchTargetPos.size();
	if (nCpt != matchSubjectPos.size() || nCpt <= 0)
	{
		fprintf(stderr, "The input vectors are invalid in compute_tps_df_field().\n");
		return 0;
	}

	Image2DSimple<MYFLOAT_JBA> * cpt_target = new Image2DSimple<MYFLOAT_JBA>(3, nCpt);
	Image2DSimple<MYFLOAT_JBA> * cpt_subject = new Image2DSimple<MYFLOAT_JBA>(3, nCpt);
	if (!cpt_target || !cpt_target->valid() || !cpt_subject || !cpt_subject->valid())
	{
		fprintf(stderr, "Fail to allocate memory.");
		if (cpt_target) { delete cpt_target; cpt_target = 0; }
		if (cpt_subject) { delete cpt_subject; cpt_subject = 0; }
		return 0;
	}

	V3DLONG n;
	Matrix x4x4_d, xnx4_c, xnxn_K;
	if (xnx4_c.nrows() != nCpt || xnx4_c.ncols() != 4)
		xnx4_c.ReSize(nCpt, 4);
	if (x4x4_d.nrows() != 4 || xnx4_c.ncols() != 4)
		x4x4_d.ReSize(4, 4);
	if (xnxn_K.nrows() != nCpt || xnxn_K.ncols() != nCpt)
		xnxn_K.ReSize(nCpt, nCpt);





	MYFLOAT_JBA ** cpt_target_ref = cpt_target->getData2dHandle();
	MYFLOAT_JBA ** cpt_subject_ref = cpt_subject->getData2dHandle();

	printf("\n---------------------------------\n");
	for (n = 0; n<nCpt; n++)
	{
		cpt_target_ref[n][0] = matchTargetPos.at(n).x;
		cpt_target_ref[n][1] = matchTargetPos.at(n).y;
		cpt_target_ref[n][2] = matchTargetPos.at(n).z;

		cpt_subject_ref[n][0] = matchSubjectPos.at(n).x;
		cpt_subject_ref[n][1] = matchSubjectPos.at(n).y;
		cpt_subject_ref[n][2] = matchSubjectPos.at(n).z;

		printf("n=%d \tx=[%5.3f -> %5.3f] y=[%5.3f -> %5.3f] z=[%5.3f -> %5.3f] \n",
			n, cpt_target_ref[n][0], cpt_subject_ref[n][0], cpt_target_ref[n][1], cpt_subject_ref[n][1], cpt_target_ref[n][2], cpt_subject_ref[n][2]);
	}
	printf("\n#################################\n");
	float *H_X, *H_Y, *H_Z;
	H_X = (float*)malloc(nCpt * sizeof(float)); H_Y = (float*)malloc(nCpt * sizeof(float)); H_Z = (float*)malloc(nCpt * sizeof(float));
	for (unsigned V3DLONG j = 0; j<nCpt; j++)
	{
		H_X[j] = cpt_subject_ref[j][0]; H_Y[j] = cpt_subject_ref[j][1]; H_Z[j] = cpt_subject_ref[j][2];
	}
	//compute K=-r=-|xi-xj|

	double d_x, d_y, d_z;
	for (unsigned V3DLONG i = 0; i<nCpt; i++)
		for (unsigned V3DLONG j = 0; j<nCpt; j++)
		{
			d_x = cpt_subject_ref[i][0] - cpt_subject_ref[j][0];
			d_y = cpt_subject_ref[i][1] - cpt_subject_ref[j][1];
			d_z = cpt_subject_ref[i][2] - cpt_subject_ref[j][2];
			xnxn_K(i + 1, j + 1) = -sqrt(d_x*d_x + d_y*d_y + d_z*d_z);
		}
//	printf("\t>>xnxn_K time consume %.6f ms\n", (float)(clock() - stps_start) / CLOCKS_PER_SEC * 1000);

	clock_t c;
	c = clock();
	if (gpu_type == 1)
	{
		Matrix X(nCpt, 4), Y(nCpt, 4);
		Matrix Q_ori(nCpt, nCpt); Q_ori = 0.0;
		for (V3DLONG i = 0; i < nCpt; i++)
		{
			Q_ori(i + 1, 1) = X(i + 1, 1) = 1;
			Q_ori(i + 1, 2) = X(i + 1, 2) = cpt_subject_ref[i][0];
			Q_ori(i + 1, 3) = X(i + 1, 3) = cpt_subject_ref[i][1];
			Q_ori(i + 1, 4) = X(i + 1, 4) = cpt_subject_ref[i][2];

			Y(i + 1, 1) = 1;
			Y(i + 1, 2) = cpt_target_ref[i][0];
			Y(i + 1, 3) = cpt_target_ref[i][1];
			Y(i + 1, 4) = cpt_target_ref[i][2];
		}
		Matrix Q(nCpt, nCpt), Q_x(nCpt, nCpt);
		clock_t stps_start1;
		stps_start1 = clock();

        // qr_gpu
		Matrix R(4, 4);
		clock_t stps_startqr;
		stps_startqr = clock();
		gpu_QR(nCpt, Q_ori, Q_x, R);
		printf("\t>>QR�ֽ� time consume %.6f s\n", (float)(clock() - stps_startqr) / CLOCKS_PER_SEC);
		printf("\t>>QR�ֽ� time consume %.6f s\n", (float)(clock() - stps_start1) / CLOCKS_PER_SEC);
		clock_t stps_start2;
		stps_start2 = clock();
		UpperTriangularMatrix R1;
		clock_t stps_start;
		stps_start = clock();
		gpu_extendornormal(nCpt, 4, Q_x);
		Q = Q_x;
		printf("\t>>extend_orthonormal���� time consume %.6f s\n", (float)(clock() - stps_start2) / CLOCKS_PER_SEC);
		Matrix q1 = Q.columns(1, 4);
		Matrix q2 = Q.columns(5, nCpt);
		Matrix r = R.submatrix(1, 4, 1, 4);
		printf("\t>>q1��q2��r���� time consume %.6f s\n", (float)(clock() - stps_start2) / CLOCKS_PER_SEC);
		clock_t stps_start3;
		stps_start3 = clock();
		Matrix KQ(nCpt, nCpt - 4);
		KQ = matrixMultiply(nCpt, nCpt - 4, nCpt, xnxn_K, q2);
		Matrix q2t = q2.t();
		Matrix A1 = matrixMultiply(nCpt - 4, nCpt - 4, nCpt, q2t, KQ);
		Matrix A = A1 + IdentityMatrix(nCpt - 4)*0.2;
		Matrix A_i(nCpt - 4, nCpt - 4); A_i = 0.0;
		gpu_A_i_new(nCpt, A, A_i);
		Matrix q2A(nCpt, nCpt - 4);
		q2A = matrixMultiply(nCpt, nCpt - 4, nCpt - 4, q2, A_i);
		Matrix c1 = matrixMultiply(nCpt, nCpt, nCpt - 4, q2A, q2t);
		xnx4_c = c1*Y;
		printf("\t>>total xnx4_c���� time consume %.6f s\n", (float)(clock() - stps_start3) / CLOCKS_PER_SEC);
		clock_t stps_start4;
		stps_start4 = clock();
		x4x4_d = r.i()*q1.t()*(Y - xnxn_K*xnx4_c);
		printf("\t>>x4x4_d���� time consume %.6f ms\n", (float)(clock() - stps_start4) / CLOCKS_PER_SEC * 1000);
		printf("\t>>the total consumption %.6f s\n", (float)(clock() - stps_start) / CLOCKS_PER_SEC);


        // qr_cpu
//        clock_t stps_startqr;
//        stps_startqr = clock();
//        clock_t stps_start2;
//        stps_start2 = clock();
//        UpperTriangularMatrix R1;
//        QRZ(Q_ori, R1);
//        Q_x = Q_ori;
//        printf("\t>>QR�ֽ� time consume %.6f s\n", (float)(clock() - stps_startqr) / CLOCKS_PER_SEC);
//        printf("\t>>QR�ֽ� time consume %.6f s\n", (float)(clock() - stps_start1) / CLOCKS_PER_SEC);
//        clock_t stps_start;
//        stps_start = clock();
//        gpu_extendornormal(nCpt, 4, Q_x);
//        Q = Q_x;
//        printf("\t>>extend_orthonormal���� time consume %.6f s\n", (float)(clock() - stps_start2) / CLOCKS_PER_SEC);
//        Matrix q1 = Q.columns(1, 4);
//        Matrix q2 = Q.columns(5, nCpt);
//        Matrix r = R1.submatrix(1, 4, 1, 4);
//        printf("\t>>q1��q2��r���� time consume %.6f s\n", (float)(clock() - stps_start2) / CLOCKS_PER_SEC);
//        clock_t stps_start3;
//        stps_start3 = clock();
//        Matrix KQ(nCpt, nCpt - 4);
//        KQ = matrixMultiply(nCpt, nCpt - 4, nCpt, xnxn_K, q2);
//        Matrix q2t = q2.t();
//        Matrix A1 = matrixMultiply(nCpt - 4, nCpt - 4, nCpt, q2t, KQ);
//        Matrix A = A1 + IdentityMatrix(nCpt - 4)*0.2;
//
//        Matrix A_i(nCpt - 4, nCpt - 4); A_i = 0.0;
//
//        gpu_A_i_new(nCpt, A, A_i);
//
//        Matrix q2A(nCpt, nCpt - 4);
//        q2A = matrixMultiply(nCpt, nCpt - 4, nCpt - 4, q2, A_i);
//        Matrix c1 = matrixMultiply(nCpt, nCpt, nCpt - 4, q2A, q2t);
//        xnx4_c = c1*Y;
//
//        printf("\t>>total xnx4_c���� time consume %.6f s\n", (float)(clock() - stps_start3) / CLOCKS_PER_SEC);
//        clock_t stps_start4;
//        stps_start4 = clock();
//        x4x4_d = r.i()*q1.t()*(Y - xnxn_K*xnx4_c);
//
//        printf("\t>>x4x4_d���� time consume %.6f ms\n", (float)(clock() - stps_start4) / CLOCKS_PER_SEC * 1000);
//        printf("\t>>the total consumption %.6f s\n", (float)(clock() - stps_start) / CLOCKS_PER_SEC);

	}
	else
	{
		//clock_t c;
		//c = clock();
		Matrix X(nCpt, 4), Y(nCpt, 4);
		Matrix Q(nCpt, nCpt); Q = 0.0;
		for (V3DLONG i = 0; i < nCpt; i++)
		{
			Q(i + 1, 1) = X(i + 1, 1) = 1;
			Q(i + 1, 2) = X(i + 1, 2) = cpt_subject_ref[i][0];
			Q(i + 1, 3) = X(i + 1, 3) = cpt_subject_ref[i][1];
			Q(i + 1, 4) = X(i + 1, 4) = cpt_subject_ref[i][2];
			Y(i + 1, 1) = 1;
			Y(i + 1, 2) = cpt_target_ref[i][0];
			Y(i + 1, 3) = cpt_target_ref[i][1];
			Y(i + 1, 4) = cpt_target_ref[i][2];
		}
		UpperTriangularMatrix R;
		QRZ(Q, R);
		clock_t stps_start;
		stps_start = clock();
		extend_orthonormal(Q, 4);//otherwise q2=0

		Matrix q1 = Q.columns(1, 4);
		Matrix q2 = Q.columns(5, nCpt);
		Matrix r = R.submatrix(1, 4, 1, 4);
		//compute non-affine term c which decomposed from TPS
		Matrix A = q2.t()*xnxn_K*q2 + IdentityMatrix(nCpt - 4)*0.2;
		xnx4_c = q2*(A.i()*q2.t()*Y);
		//compute affine term d (normal)
		x4x4_d = r.i()*q1.t()*(Y - xnxn_K*xnx4_c);
		printf("\t>>xnxn_K time consume %.2f s\n", (float)(clock() - stps_start) / CLOCKS_PER_SEC);
	}


	V3DLONG p;

	//	V3DLONG gsz0 = (V3DLONG)(ceil((double(sz0)/gfactor_x)))+1, gsz1 = (V3DLONG)(ceil((double(sz1)/gfactor_y)))+1, gsz2 = (V3DLONG)(ceil((double(sz2)/gfactor_z)))+1;
	V3DLONG gsz0 = (V3DLONG)(ceil((double(sz0) / gfactor_x))) + 1 + 2, gsz1 = (V3DLONG)(ceil((double(sz1) / gfactor_y))) + 1 + 2, gsz2 = (V3DLONG)(ceil((double(sz2) / gfactor_z))) + 1 + 2;//+2 for bspline
	Vol3DSimple<DisplaceFieldF3D> * df_local = new Vol3DSimple<DisplaceFieldF3D>(gsz0, gsz1, gsz2);
	DisplaceFieldF3D *** df_local_3d = df_local->getData3dHandle();

	if (!df_local || !df_local->valid())
	{
		fprintf(stderr, "Fail to allocate memory for the subsampled DF volume memory [%d].\n", __LINE__);

		if (cpt_target) { delete cpt_target; cpt_target = 0; }
		if (cpt_subject) { delete cpt_subject; cpt_subject = 0; }
		if (df_local) { delete df_local; df_local = 0; }
		return 0;
	}

	printf("-------------------- Now compute the distances of pixels to the mapping points. -------\n\n");

	V3DLONG i, j, k;
	DisplaceFieldF3D * df_local_1d = df_local->getData1dHandle();
	for (k = 0; k<df_local->getTotalElementNumber(); k++)
	{
		df_local_1d[k].sz = df_local_1d[k].sy = df_local_1d[k].sx = 0;
	}
	clock_t stps_computedistance;
	stps_computedistance = clock();
	if (gpu_type == 1)
	    gpu_computedistance(nCpt, gsz2, gsz1, gsz0, gfactor_x, gfactor_y, gfactor_z, x4x4_d, xnx4_c, H_X, H_Y, H_Z, df_local_3d, 0, 0, 0);
	else
    {
        for (k = 0; k<gsz2; k++)
        {
            for (j = 0; j<gsz1; j++)
            {
                for (i = 0; i<gsz0; i++)
                {
                    Matrix x_ori(1, 4);
                    x_ori(1, 1) = 1.0;
                    x_ori(1, 2) = (i - 1)*gfactor_x;
                    x_ori(1, 3) = (j - 1)*gfactor_y;
                    x_ori(1, 4) = (k - 1)*gfactor_z;
                    Matrix x_stps(1, 4);
                    Matrix xmxn_K;
                    xmxn_K.resize(1, nCpt);
                    double d_x, d_y, d_z;
                    for (unsigned V3DLONG n = 0; n<nCpt; n++)
                    {
                        d_x = (i - 1)*gfactor_x - H_X[n];
                        d_y = (j - 1)*gfactor_y - H_Y[n];
                        d_z = (k - 1)*gfactor_z - H_Z[n];
                        xmxn_K(1, n + 1) = -sqrt(d_x*d_x + d_y*d_y + d_z*d_z);
                    }
                    x_stps = x_ori*x4x4_d + xmxn_K*xnx4_c;
                    df_local_3d[k][j][i].sx = x_stps(1, 2) - ((i - 1)*gfactor_x);
                    df_local_3d[k][j][i].sy = x_stps(1, 3) - ((j - 1)*gfactor_y);
                    df_local_3d[k][j][i].sz = x_stps(1, 4) - ((k - 1)*gfactor_z);
                }
            }
        }
    }

	printf("\t>>computedistance time consume %.2f s\n", (float)(clock() - stps_computedistance) / CLOCKS_PER_SEC);
	printf("\n");

	free(H_X);
	free(H_Y);
	free(H_Z);
	if (cpt_target) { delete cpt_target; cpt_target = 0; }
	if (cpt_subject) { delete cpt_subject; cpt_subject = 0; }

	return df_local;
}

bool q_nonrigid_ini_bsplinebasis_3D(const long n, Matrix &BxBxB)
{
	//check paras
	if (n <= 0)
	{
		printf("ERROR: n should > 0!\n");
		return false;
	}

	//cubic B-spline basis matrix
	Matrix B(4, 4);
	B.row(1) << -1 << 3 << -3 << 1;
	B.row(2) << 3 << -6 << 3 << 0;
	B.row(3) << -3 << 0 << 3 << 0;
	B.row(4) << 1 << 4 << 1 << 0;
	B /= 6.0;

	//construct T(i,:)=[t^3 t^2 t^1 1]
	Matrix T(n, 4);
	double t_step = 1.0 / n;
	for (long i = 0; i<n; i++)
	{
		double t = t_step*i;
		for (long j = 0; j <= 3; j++)
			T(i + 1, j + 1) = pow(t, 3 - j);
	}

	//construct B-spline basis/blending functions B=T*B
	Matrix TB = T*B;//n x 4

	//construct B-spline basis/blending functions for 2D interpolation B=BxB
	Matrix BxB = KP(TB, TB);//n^2 x 4^2
	//construct B-spline basis/blending functions for 3D interpolation B=BxBxB
	BxBxB = KP(BxB, TB);//n^3 x 4^3

	return true;
}

Vol3DSimple <MYFLOAT_JBA> * linearinterp_regularmesh_3d(V3DLONG sz0, V3DLONG sz1, V3DLONG sz2, Vol3DSimple <MYFLOAT_JBA> * df_regular_grid)
{
	V3DLONG k, j, i;

	if (!df_regular_grid || !df_regular_grid->valid())
	{
		fprintf(stderr, "The pointer is not correct.\n");
		return 0;
	}
	MYFLOAT_JBA *** df_grid3d = df_regular_grid->getData3dHandle();
	//	V3DLONG n0 = df_regular_grid->sz0()-3, n1 = df_regular_grid->sz1()-3, n2 = df_regular_grid->sz2()-3;//-3 for B-spline?
	V3DLONG n0 = df_regular_grid->sz0() - 1, n1 = df_regular_grid->sz1() - 1, n2 = df_regular_grid->sz2() - 1;//modified by qul @ 120710
	if (n0 <= 0 || n1 <= 0 || n2 <= 0)
	{
		fprintf(stderr, "The size  is not correct.\n");
		return 0;
	}
	if (sz0 <= 0 || sz1 <= 0 || sz2 <= 0)
	{
		fprintf(stderr, "The size of the DF to be computed is not correct.\n");
		return 0;
	}

	Vol3DSimple <MYFLOAT_JBA> * df_field = new Vol3DSimple <MYFLOAT_JBA>(sz0, sz1, sz2);
	if (!df_field || !df_field->valid())
	{
		fprintf(stderr, "Fail to allocate memory.\n");
		if (df_field) { delete df_field; df_field = 0; }
		return 0;
	}

	MYFLOAT_JBA * df_field_ref1d = df_field->getData1dHandle();
	for (i = 0; i<df_field->getTotalElementNumber(); i++)
	{
		df_field_ref1d[i] = 0;
	}

	Coord3D_JBA *c = new Coord3D_JBA[df_field->getTotalElementNumber()];
	double nf0 = (double)n0 / sz0, nf1 = (double)n1 / sz1, nf2 = (double)n2 / sz2;
	V3DLONG cnt = 0;
	for (k = 0; k<sz2; k++)
	{
		double k_tmp = (double)k*nf2;
		for (j = 0; j<sz1; j++)
		{
			double j_tmp = (double)j*nf1;
			for (i = 0; i<sz0; i++)
			{
				c[cnt].x = i*nf0;
				c[cnt].y = j_tmp;
				c[cnt].z = k_tmp;
				cnt++;
			}
		}
	}

	interpolate_coord_linear(df_field_ref1d, c, df_field->getTotalElementNumber(),
		df_grid3d, df_regular_grid->sz0(), df_regular_grid->sz1(), df_regular_grid->sz2(),
		0, n0 - 1 + 1, 0, n1 - 1 + 1, 0, n2 - 1 + 1);

	if (c) { delete[]c; c = 0; }
	return df_field;
}

bool interpolate_coord_linear(MYFLOAT_JBA * interpolatedVal, Coord3D_JBA *c, V3DLONG numCoord,
	MYFLOAT_JBA *** templateVol3d, V3DLONG tsz0, V3DLONG tsz1, V3DLONG tsz2,
	V3DLONG tlow0, V3DLONG tup0, V3DLONG tlow1, V3DLONG tup1, V3DLONG tlow2, V3DLONG tup2)
{
	if (!interpolatedVal || !c || numCoord <= 0 ||
		!templateVol3d || tsz0 <= 0 || tsz1 <= 0 || tsz2 <= 0 ||
		tlow0<0 || tlow0 >= tsz0 || tup0<0 || tup0 >= tsz0 || tlow0>tup0 ||
		tlow1<0 || tlow1 >= tsz1 || tup1<0 || tup1 >= tsz1 || tlow1>tup1 ||
		tlow2<0 || tlow2 >= tsz2 || tup2<0 || tup2 >= tsz2 || tlow2>tup2)
	{
		fprintf(stderr, "Invalid parameters! [%s][%d]\n", __FILE__, __LINE__);
		return false;
	}

	double curpx, curpy, curpz;
	V3DLONG cpx0, cpx1, cpy0, cpy1, cpz0, cpz1;

	for (V3DLONG ipt = 0; ipt<numCoord; ipt++)
	{
		if (c[ipt].x< tlow0 || c[ipt].x> tup0 || c[ipt].y< tlow1 || c[ipt].y> tup1 || c[ipt].z< tlow2 || c[ipt].z> tup2)
		{
			interpolatedVal[ipt] = 0;
			continue;
		}

		curpx = c[ipt].x; curpx = (curpx<tlow0) ? tlow0 : curpx; curpx = (curpx>tup0) ? tup0 : curpx;
#ifndef POSITIVE_Y_COORDINATE
		curpy = tsz1 - 1 - c[ipt].y; curpy = (curpy<tlow1) ? tlow1 : curpy; curpy = (curpy>tup1) ? tup1 : curpy;
#else
		curpy = c[ipt].y; curpy = (curpy<tlow1) ? tlow1 : curpy; curpy = (curpy>tup1) ? tup1 : curpy;
#endif
		curpz = c[ipt].z; curpz = (curpz<tlow2) ? tlow2 : curpz; curpz = (curpz>tup2) ? tup2 : curpz;

		cpx0 = V3DLONG(floor(curpx)); cpx1 = V3DLONG(ceil(curpx));
		cpy0 = V3DLONG(floor(curpy)); cpy1 = V3DLONG(ceil(curpy));
		cpz0 = V3DLONG(floor(curpz)); cpz1 = V3DLONG(ceil(curpz));

		if (cpz0 == cpz1)
		{
			if (cpy0 == cpy1)
			{
				if (cpx0 == cpx1)
				{
					interpolatedVal[ipt] = (MYFLOAT_JBA)(templateVol3d[cpz0][cpy0][cpx0]);
				}
				else
				{
					double w0x0y0z = (cpx1 - curpx);
					double w1x0y0z = (curpx - cpx0);
					interpolatedVal[ipt] = (MYFLOAT_JBA)(w0x0y0z * double(templateVol3d[cpz0][cpy0][cpx0]) +
						w1x0y0z * double(templateVol3d[cpz0][cpy0][cpx1]));
				}
			}
			else
			{
				if (cpx0 == cpx1)
				{
					double w0x0y0z = (cpy1 - curpy);
					double w0x1y0z = (curpy - cpy0);
					interpolatedVal[ipt] = (MYFLOAT_JBA)(w0x0y0z * double(templateVol3d[cpz0][cpy0][cpx0]) +
						w0x1y0z * double(templateVol3d[cpz0][cpy1][cpx0]));
				}
				else
				{
					double w0x0y0z = (cpx1 - curpx)*(cpy1 - curpy);
					double w0x1y0z = (cpx1 - curpx)*(curpy - cpy0);
					double w1x0y0z = (curpx - cpx0)*(cpy1 - curpy);
					double w1x1y0z = (curpx - cpx0)*(curpy - cpy0);
					interpolatedVal[ipt] = (MYFLOAT_JBA)(w0x0y0z * double(templateVol3d[cpz0][cpy0][cpx0]) +
						w0x1y0z * double(templateVol3d[cpz0][cpy1][cpx0]) +
						w1x0y0z * double(templateVol3d[cpz0][cpy0][cpx1]) +
						w1x1y0z * double(templateVol3d[cpz0][cpy1][cpx1]));
				}
			}
		}
		else
		{
			if (cpy0 == cpy1)
			{
				if (cpx0 == cpx1)
				{
					double w0x0y0z = (cpz1 - curpz);
					double w0x0y1z = (curpz - cpz0);

					interpolatedVal[ipt] = (MYFLOAT_JBA)(w0x0y0z * double(templateVol3d[cpz0][cpy0][cpx0]) + w0x0y1z * double(templateVol3d[cpz1][cpy0][cpx0]));
				}
				else
				{
					double w0x0y0z = (cpx1 - curpx)*(cpz1 - curpz);
					double w0x0y1z = (cpx1 - curpx)*(curpz - cpz0);

					double w1x0y0z = (curpx - cpx0)*(cpz1 - curpz);
					double w1x0y1z = (curpx - cpx0)*(curpz - cpz0);

					interpolatedVal[ipt] = (MYFLOAT_JBA)(w0x0y0z * double(templateVol3d[cpz0][cpy0][cpx0]) + w0x0y1z * double(templateVol3d[cpz1][cpy0][cpx0]) +
						w1x0y0z * double(templateVol3d[cpz0][cpy0][cpx1]) + w1x0y1z * double(templateVol3d[cpz1][cpy0][cpx1]));
				}
			}
			else
			{
				if (cpx0 == cpx1)
				{
					double w0x0y0z = (cpy1 - curpy)*(cpz1 - curpz);
					double w0x0y1z = (cpy1 - curpy)*(curpz - cpz0);

					double w0x1y0z = (curpy - cpy0)*(cpz1 - curpz);
					double w0x1y1z = (curpy - cpy0)*(curpz - cpz0);

					interpolatedVal[ipt] = (MYFLOAT_JBA)(w0x0y0z * double(templateVol3d[cpz0][cpy0][cpx0]) + w0x0y1z * double(templateVol3d[cpz1][cpy0][cpx0]) +
						w0x1y0z * double(templateVol3d[cpz0][cpy1][cpx0]) + w0x1y1z * double(templateVol3d[cpz1][cpy1][cpx0]));
				}
				else
				{
					double w0x0y0z = (cpx1 - curpx)*(cpy1 - curpy)*(cpz1 - curpz);
					double w0x0y1z = (cpx1 - curpx)*(cpy1 - curpy)*(curpz - cpz0);

					double w0x1y0z = (cpx1 - curpx)*(curpy - cpy0)*(cpz1 - curpz);
					double w0x1y1z = (cpx1 - curpx)*(curpy - cpy0)*(curpz - cpz0);

					double w1x0y0z = (curpx - cpx0)*(cpy1 - curpy)*(cpz1 - curpz);
					double w1x0y1z = (curpx - cpx0)*(cpy1 - curpy)*(curpz - cpz0);

					double w1x1y0z = (curpx - cpx0)*(curpy - cpy0)*(cpz1 - curpz);
					double w1x1y1z = (curpx - cpx0)*(curpy - cpy0)*(curpz - cpz0);

					interpolatedVal[ipt] = (MYFLOAT_JBA)(w0x0y0z * double(templateVol3d[cpz0][cpy0][cpx0]) + w0x0y1z * double(templateVol3d[cpz1][cpy0][cpx0]) +
						w0x1y0z * double(templateVol3d[cpz0][cpy1][cpx0]) + w0x1y1z * double(templateVol3d[cpz1][cpy1][cpx0]) +
						w1x0y0z * double(templateVol3d[cpz0][cpy0][cpx1]) + w1x0y1z * double(templateVol3d[cpz1][cpy0][cpx1]) +
						w1x1y0z * double(templateVol3d[cpz0][cpy1][cpx1]) + w1x1y1z * double(templateVol3d[cpz1][cpy1][cpx1]));
				}
			}
		}

	}

	return true;
}

//The following functions created by CYF 2022.9.28
//Big warp
//
bool compute_df_stps_subsampled_volume_4bspline_per(const vector <Coord3D_JBA> & matchTargetPos, const vector <Coord3D_JBA> & matchSubjectPos, Matrix &x4x4_d, Matrix &xnx4_c, float * &H_X, float * &H_Y,
	float * &H_Z, int nCpt, Image2DSimple<MYFLOAT_JBA> * &cpt_subject,int gpu_type)
{
	//nCpt = matchTargetPos.size();
	if (nCpt != matchSubjectPos.size() || nCpt <= 0)
	{
		fprintf(stderr, "The input vectors are invalid in compute_tps_df_field().\n");
		return 0;
	}

	Image2DSimple<MYFLOAT_JBA> * cpt_target = new Image2DSimple<MYFLOAT_JBA>(3, nCpt);
	//Image2DSimple<MYFLOAT_JBA> * cpt_subject = new Image2DSimple<MYFLOAT_JBA>(3, nCpt);
	cpt_subject = new Image2DSimple<MYFLOAT_JBA>(3, nCpt);
	if (!cpt_target || !cpt_target->valid() || !cpt_subject || !cpt_subject->valid())
	{
		fprintf(stderr, "Fail to allocate memory.");
		if (cpt_target) { delete cpt_target; cpt_target = 0; }
		if (cpt_subject) { delete cpt_subject; cpt_subject = 0; }
		return 0;
	}

	V3DLONG n;
	//Matrix x4x4_d, xnx4_c, xnxn_K;
	Matrix xnxn_K;
	if (xnx4_c.nrows() != nCpt || xnx4_c.ncols() != 4)
		xnx4_c.ReSize(nCpt, 4);
	if (x4x4_d.nrows() != 4 || xnx4_c.ncols() != 4)
		x4x4_d.ReSize(4, 4);
	if (xnxn_K.nrows() != nCpt || xnxn_K.ncols() != nCpt)
		xnxn_K.ReSize(nCpt, nCpt);





	MYFLOAT_JBA ** cpt_target_ref = cpt_target->getData2dHandle();
	MYFLOAT_JBA ** cpt_subject_ref = cpt_subject->getData2dHandle();

//	printf("\n---------------------------------\n");
	for (n = 0; n<nCpt; n++)
	{
		cpt_target_ref[n][0] = matchTargetPos.at(n).x;
		cpt_target_ref[n][1] = matchTargetPos.at(n).y;
		cpt_target_ref[n][2] = matchTargetPos.at(n).z;

		cpt_subject_ref[n][0] = matchSubjectPos.at(n).x;
		cpt_subject_ref[n][1] = matchSubjectPos.at(n).y;
		cpt_subject_ref[n][2] = matchSubjectPos.at(n).z;

//		printf("n=%d \tx=[%5.3f -> %5.3f] y=[%5.3f -> %5.3f] z=[%5.3f -> %5.3f] \n",
//			n, cpt_target_ref[n][0], cpt_subject_ref[n][0], cpt_target_ref[n][1], cpt_subject_ref[n][1], cpt_target_ref[n][2], cpt_subject_ref[n][2]);
	}
//	printf("\n#################################\n");
	//float *H_X, *H_Y, *H_Z;
	H_X = (float*)malloc(nCpt * sizeof(float)); H_Y = (float*)malloc(nCpt * sizeof(float)); H_Z = (float*)malloc(nCpt * sizeof(float));
	for (unsigned V3DLONG j = 0; j<nCpt; j++)
	{
		H_X[j] = cpt_subject_ref[j][0]; H_Y[j] = cpt_subject_ref[j][1]; H_Z[j] = cpt_subject_ref[j][2];
	}
	//compute K=-r=-|xi-xj|

	double d_x, d_y, d_z;
	for (unsigned V3DLONG i = 0; i<nCpt; i++)
		for (unsigned V3DLONG j = 0; j<nCpt; j++)
		{
			d_x = cpt_subject_ref[i][0] - cpt_subject_ref[j][0];
			d_y = cpt_subject_ref[i][1] - cpt_subject_ref[j][1];
			d_z = cpt_subject_ref[i][2] - cpt_subject_ref[j][2];
			xnxn_K(i + 1, j + 1) = -sqrt(d_x*d_x + d_y*d_y + d_z*d_z);
		}
	//	printf("\t>>xnxn_K time consume %.6f ms\n", (float)(clock() - stps_start) / CLOCKS_PER_SEC * 1000);
	//int cdd = 0;
	clock_t c;
	c = clock();

	if (gpu_type==1)
	{
		Matrix X(nCpt, 4), Y(nCpt, 4);
		Matrix Q_ori(nCpt, nCpt); Q_ori = 0.0;
		for (V3DLONG i = 0; i < nCpt; i++)
		{
			Q_ori(i + 1, 1) = X(i + 1, 1) = 1;
			Q_ori(i + 1, 2) = X(i + 1, 2) = cpt_subject_ref[i][0];
			Q_ori(i + 1, 3) = X(i + 1, 3) = cpt_subject_ref[i][1];
			Q_ori(i + 1, 4) = X(i + 1, 4) = cpt_subject_ref[i][2];

			Y(i + 1, 1) = 1;
			Y(i + 1, 2) = cpt_target_ref[i][0];
			Y(i + 1, 3) = cpt_target_ref[i][1];
			Y(i + 1, 4) = cpt_target_ref[i][2];
		}
		Matrix Q(nCpt, nCpt), Q_x(nCpt, nCpt);
		clock_t stps_start1;
		stps_start1 = clock();

		clock_t stps_startqr;
		stps_startqr = clock();

        // qr_gpu
        Matrix R(4, 4);
		gpu_QR(nCpt, Q_ori, Q_x, R);
        printf("\t>>QR分解 time consume %.6f s\n", (float)(clock() - stps_startqr) / CLOCKS_PER_SEC);
        Matrix r = R.submatrix(1, 4, 1, 4);

        // qr_cpu
//        UpperTriangularMatrix R1;
//		clock_t stps_start2;
//		stps_start2 = clock();
//        QRZ(Q_ori, R1);
//        Q_x = Q_ori;
//        Matrix r = R1.submatrix(1, 4, 1, 4);


		// continue
        clock_t stps_start;
        stps_start = clock();
		gpu_extendornormal(nCpt, 4, Q_x);
		Q = Q_x;
		Matrix q1 = Q.columns(1, 4);
		Matrix q2 = Q.columns(5, nCpt);
		clock_t stps_start3;
		stps_start3 = clock();
		Matrix KQ(nCpt, nCpt - 4);
		KQ = matrixMultiply(nCpt, nCpt - 4, nCpt, xnxn_K, q2);
		Matrix q2t = q2.t();
		Matrix A1 = matrixMultiply(nCpt - 4, nCpt - 4, nCpt, q2t, KQ);
		Matrix A = A1 + IdentityMatrix(nCpt - 4)*0.2;
		Matrix A_i(nCpt - 4, nCpt - 4); A_i = 0.0;
		gpu_A_i_new(nCpt, A, A_i);
		Matrix q2A(nCpt, nCpt - 4);
		q2A = matrixMultiply(nCpt, nCpt - 4, nCpt - 4, q2, A_i);
		Matrix c1 = matrixMultiply(nCpt, nCpt, nCpt - 4, q2A, q2t);
		xnx4_c = c1*Y;
		clock_t stps_start4;
		stps_start4 = clock();
		x4x4_d = r.i()*q1.t()*(Y - xnxn_K*xnx4_c);

	}
	else
	{
		//clock_t c;
		//c = clock();
		Matrix X(nCpt, 4), Y(nCpt, 4);
		Matrix Q(nCpt, nCpt); Q = 0.0;
		for (V3DLONG i = 0; i < nCpt; i++)
		{
			Q(i + 1, 1) = X(i + 1, 1) = 1;
			Q(i + 1, 2) = X(i + 1, 2) = cpt_subject_ref[i][0];
			Q(i + 1, 3) = X(i + 1, 3) = cpt_subject_ref[i][1];
			Q(i + 1, 4) = X(i + 1, 4) = cpt_subject_ref[i][2];
			Y(i + 1, 1) = 1;
			Y(i + 1, 2) = cpt_target_ref[i][0];
			Y(i + 1, 3) = cpt_target_ref[i][1];
			Y(i + 1, 4) = cpt_target_ref[i][2];
		}
		UpperTriangularMatrix R;
		QRZ(Q, R);
		clock_t stps_start;
		stps_start = clock();
		extend_orthonormal(Q, 4);//otherwise q2=0

		Matrix q1 = Q.columns(1, 4);
		Matrix q2 = Q.columns(5, nCpt);
		Matrix r = R.submatrix(1, 4, 1, 4);
		//compute non-affine term c which decomposed from TPS
		Matrix A = q2.t()*xnxn_K*q2 + IdentityMatrix(nCpt - 4)*0.2;
		xnx4_c = q2*(A.i()*q2.t()*Y);
		//compute affine term d (normal)
		x4x4_d = r.i()*q1.t()*(Y - xnxn_K*xnx4_c);
//		printf("\t>>xnxn_K time consume %.2f s\n", (float)(clock() - stps_start) / CLOCKS_PER_SEC);
	}

	if (cpt_target) { delete cpt_target; cpt_target = 0; }


	return true;
}

Vol3DSimple<DisplaceFieldF3D> * compute_df_stps_subsampled_volume_4bspline_block(int nCpt ,Matrix x4x4_d, Matrix xnx4_c, V3DLONG sz0, V3DLONG sz1, V3DLONG sz2, V3DLONG gfactor_x,
	V3DLONG gfactor_y, V3DLONG gfactor_z, float * H_X, float * H_Y, float * H_Z, long long x_offset, long long y_offset, long long z_offset, int gpu_mode)
{
	V3DLONG p;

	//	V3DLONG gsz0 = (V3DLONG)(ceil((double(sz0)/gfactor_x)))+1, gsz1 = (V3DLONG)(ceil((double(sz1)/gfactor_y)))+1, gsz2 = (V3DLONG)(ceil((double(sz2)/gfactor_z)))+1;
	V3DLONG gsz0 = (V3DLONG)(ceil((double(sz0) / gfactor_x))) + 1 + 2, gsz1 = (V3DLONG)(ceil((double(sz1) / gfactor_y))) + 1 + 2, gsz2 = (V3DLONG)(ceil((double(sz2) / gfactor_z))) + 1 + 2;//+2 for bspline
	Vol3DSimple<DisplaceFieldF3D> * df_local = new Vol3DSimple<DisplaceFieldF3D>(gsz0, gsz1, gsz2);
	DisplaceFieldF3D *** df_local_3d = df_local->getData3dHandle();

	if (!df_local || !df_local->valid())
	{
		fprintf(stderr, "Fail to allocate memory for the subsampled DF volume memory [%d].\n", __LINE__);
		if (df_local) { delete df_local; df_local = 0; }
		return 0;
	}



//	printf("-------------------- Now compute the distances of pixels to the mapping points. -------\n\n");

	V3DLONG i, j, k;
	DisplaceFieldF3D * df_local_1d = df_local->getData1dHandle();
	for (k = 0; k<df_local->getTotalElementNumber(); k++)
	{
		df_local_1d[k].sz = df_local_1d[k].sy = df_local_1d[k].sx = 0;
	}
	clock_t stps_computedistance;
	stps_computedistance = clock();
	if (gpu_mode==1)
	    gpu_computedistance(nCpt, gsz2, gsz1, gsz0, gfactor_x, gfactor_y, gfactor_z, x4x4_d, xnx4_c, H_X, H_Y, H_Z, df_local_3d, x_offset, y_offset, z_offset);
	else
    {
#pragma omp parallel for
        for (k = 0; k<gsz2; k++)
        {
        	for (j = 0; j<gsz1; j++)
        	{
        		for (i = 0; i<gsz0; i++)
        		{
        			Matrix x_ori(1, 4);
        			x_ori(1, 1) = 1.0;
        			x_ori(1, 2) = (i - 1)*gfactor_x + x_offset;
        			x_ori(1, 3) = (j - 1)*gfactor_y + y_offset;
        			x_ori(1, 4) = (k - 1)*gfactor_z + z_offset;
        			Matrix x_stps(1, 4);
        			Matrix xmxn_K;
        			xmxn_K.resize(1, nCpt);
        			double d_x, d_y, d_z;
        			for (unsigned V3DLONG n = 0; n<nCpt; n++)
        			{
        				d_x = (i - 1)*gfactor_x + x_offset - H_X[n];
        				d_y = (j - 1)*gfactor_y + y_offset - H_Y[n];
        				d_z = (k - 1)*gfactor_z + z_offset - H_Z[n];
        				xmxn_K(1, n + 1) = -sqrt(d_x*d_x + d_y*d_y + d_z*d_z);
        			}
        			x_stps = x_ori*x4x4_d + xmxn_K*xnx4_c;
        			df_local_3d[k][j][i].sx = x_stps(1, 2) - ((i - 1)*gfactor_x + x_offset);
        			df_local_3d[k][j][i].sy = x_stps(1, 3) - ((j - 1)*gfactor_y + y_offset);
        			df_local_3d[k][j][i].sz = x_stps(1, 4) - ((k - 1)*gfactor_z + z_offset);
        		}
        	}
        }
    }
//	printf("\t>>computedistance time consume %.2f s\n", (float)(clock() - stps_computedistance) / CLOCKS_PER_SEC);
	

	return df_local;
}




// TPS_linear_blockbyblock image warping
//	i_interp_method_df:  0-trilinear, 1-bspline
//	i_interp_method_img: 0-trilinear, 1-nearest neighbor
//template <class T>
bool imgwarp_smallmemory_CYF(const QList<ImageMarker> &ql_marker_tar, const QList<ImageMarker> &ql_marker_sub,
V3DLONG szBlock_x, V3DLONG szBlock_y, V3DLONG szBlock_z, int i_interpmethod_df, int i_interpmethod_img, string *path_TIF, long long sz_img_ori[4], string TeraWarp_save_base, int GPU_type)
{
	//check parameters
	if (ql_marker_tar.size() == 0 || ql_marker_sub.size() == 0 || ql_marker_tar.size() != ql_marker_sub.size())
	{
		printf("ERROR: target or subject control points is invalid!\n");
		return false;
	}
	if (szBlock_x <= 0 || szBlock_y <= 0 || szBlock_z <= 0)
	{
		printf("ERROR: block size is invalid!\n");
		return false;
	}
	if (i_interpmethod_df != 0 && i_interpmethod_df != 1)
	{
		printf("ERROR: DF_interp_method should be 0(linear) or 1(bspline)!\n");
		return false;
	}
	if (i_interpmethod_img != 0 && i_interpmethod_img != 1)
	{
		printf("ERROR: img_interp_method should be 0(linear) or 1(nn)!\n");
		return false;
	}
	if (i_interpmethod_df == 1 && (szBlock_x != szBlock_y || szBlock_x != szBlock_z))
	{
		printf("ERROR: df_interp_method=bspline need szBlock_x=szBlock_y=szBlock_z!\n");
		return false;
	}

	//------------------------------------------------------------------------------------------------------------------------------------
	printf(">>>>compute the subsampled displace field \n");
	vector<Coord3D_JBA> matchTargetPos, matchSubjectPos;
	int x_para = 4;
	int y_para = 4;
	int z_para = 4;
	for (V3DLONG i = 0; i<ql_marker_tar.size(); i++)
	{
		Coord3D_JBA tmpc;
		tmpc.x = x_para * ql_marker_tar.at(i).x;	tmpc.y = y_para * ql_marker_tar.at(i).y;	tmpc.z = z_para * ql_marker_tar.at(i).z;
		matchTargetPos.push_back(tmpc);
		tmpc.x = x_para * ql_marker_sub.at(i).x;	tmpc.y = y_para * ql_marker_sub.at(i).y;	tmpc.z = z_para * ql_marker_sub.at(i).z;
		matchSubjectPos.push_back(tmpc);
	}
	int nCpt = matchTargetPos.size();
	Image2DSimple<MYFLOAT_JBA> * cpt_subject = 0;
	Matrix xnx4_c(nCpt, 4);
	Matrix x4x4_d(4, 4);
	float *H_X = 0;
	float *H_Y = 0;
	float *H_Z = 0;
	clock_t BSP;
	BSP = clock();
	if (!(compute_df_stps_subsampled_volume_4bspline_per(matchTargetPos, matchSubjectPos, x4x4_d, xnx4_c, H_X, H_Y, H_Z, nCpt, cpt_subject,GPU_type)))
	{
		printf("ERROR:compute_df_stps_subsampled_volume_4bspline_per() return false. \n");
		return false;
	}
	printf("\t>>BSP time consume: %.2f s\n", (float)(clock() - BSP) / CLOCKS_PER_SEC);
	float bsp_time = (float)(clock() - BSP) / CLOCKS_PER_SEC;
	V3DLONG sz_gridwnd = szBlock_x;
	Matrix x_bsplinebasis(pow(double(sz_gridwnd), 3.0), pow(4.0, 3.0));
	if (!q_nonrigid_ini_bsplinebasis_3D(sz_gridwnd, x_bsplinebasis))
	{
		printf("ERROR: q_ini_bsplinebasis_3D() return false!\n");
		//if (p_img_warp_4d) 		{ delete4dpointer(p_img_warp_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
		//if (p_img_sub_4d) 		{ delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
		//if (pSubDF)				{ delete pSubDF;					pSubDF = 0; }
		return false;
	}
	printf("\t>>x_bsplinebasis:[%d,%d]\n", x_bsplinebasis.nrows(), x_bsplinebasis.ncols());
	//--------------------------------------------------------------------�ֿ�warp-----------------------------------------------------------
	//��������ѭ����ȡCYF2022.6.15


	//����WARP���ַ
    std::stringstream base_block_path;
    std::stringstream Slice_path;
    Slice_path << TeraWarp_save_base << "_stps_slice";
    std::stringstream block_path;
    block_path << TeraWarp_save_base << "_stps_block";
    std::stringstream abs_pos_z;
    std::stringstream abs_pos_y;
    std::stringstream abs_pos_x;
//    if (0 != access(block_path.str().c_str(), 0))
//    {
//        mkdir(block_path.str().c_str());
//    }
//    if (0 != access(Slice_path.str().c_str(), 0))
//    {
//        mkdir(Slice_path.str().c_str());
//    }
    MkdirWithPath(block_path.str().c_str());
    MkdirWithPath(Slice_path.str().c_str());

	long long sz_img_block_size[4] = { 100, 100, 100, 1 };
	//long long sz_img_block_size_fake[4] = { 100, 100, 100, 1 };
	long long sz_img_block_size_sort[4] = { 100, 100, 100, 1 };
	long long *sz_img_read_block = 0;
	long long *sz_img_read_end_block = 0;

	int X_SIZE, Y_SIZE, Z_SIZE;

	long long per_block_siza = 1000;//����һ warp���ɵĿ�Ĵ�С
	long long per_sz_block = 2500;//������ ��ȡ�Ŀ�Ĵ�С
	long long Z_Slice_control = 500;//����ÿһ��ƴ��ʱ�Ĳ���


	Z_SIZE = sz_img_ori[2] % per_block_siza == 0 ? sz_img_ori[2] / per_block_siza : sz_img_ori[2] / per_block_siza + 1;
	Y_SIZE = sz_img_ori[1] % per_block_siza == 0 ? sz_img_ori[1] / per_block_siza : sz_img_ori[1] / per_block_siza + 1;
	X_SIZE = sz_img_ori[0] % per_block_siza == 0 ? sz_img_ori[0] / per_block_siza : sz_img_ori[0] / per_block_siza + 1;
	int sum_SIZE = X_SIZE*Y_SIZE*Z_SIZE;
	long long *sz_offset_x = 0;//ƴ�Ӻ�������Ҫ������,2022.8.29�Ѿ����Թ���������ʹ��
	long long *sz_offset_y = 0;
	sz_offset_x = new long long[X_SIZE*Y_SIZE*Z_SIZE];
	sz_offset_y = new long long[X_SIZE*Y_SIZE*Z_SIZE];

	int x_th = 0;
	int read_offset = 60;
	//control
	int warp_type = 1;
	int type_PJ = 1;

	//int real_type = 0;
	int cal_mode;

	clock_t aff_warp;
	double read_time;
	clock_t save_time;
    double PJ_time;
	clock_t sort_time;
	clock_t aff_warp1;
	vector<long long> sort_point_x;
	vector<long long> sort_point_y;
	vector<long long> sort_point_z;

	float read_sum = 0;
	float stps_sum = 0;
	float save_sum = 0;
	
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
					//���������
					long long up_point[3];
					long long low_point[3];
					//�ĸ�����
					long long front_left[3];
					long long front_right[3];
					long long back_left[3];
					long long back_right[3];


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
					//���������
					long long *read_up_point = 0;
					long long *read_low_point = 0;
					//�ĸ�����
					long long *read_front_left = 0;
					long long *read_front_right = 0;
					long long *read_back_left = 0;
					long long *read_back_right = 0;
					Vol3DSimple<DisplaceFieldF3D> *pSubDF = 0;
					aff_warp1 = clock();
					pSubDF = compute_df_stps_subsampled_volume_4bspline_block(nCpt, x4x4_d, xnx4_c, sz_img_block_size[0], sz_img_block_size[1], sz_img_block_size[2], szBlock_x,
						szBlock_y, szBlock_z, H_X, H_Y, H_Z, x_offset, y_offset, z_offset, GPU_type);
					stps_sum += (float)(clock() - aff_warp1) / CLOCKS_PER_SEC;
						
					//printf("\t>>BSP time consume: %.2f s\n", (float)(clock() - BSP) / CLOCKS_PER_SEC);
					if (!pSubDF)
					{
						printf("Fail to produce the subsampled DF.\n");
						return false;
					}
					DisplaceFieldF3D ***pppSubDF = pSubDF->getData3dHandle();
					printf("subsampled DF size: [%ld,%ld,%ld]\n", pSubDF->sz0(), pSubDF->sz1(), pSubDF->sz2());

					Vol3DSimple<DisplaceFieldF3D> *pDFBlock = new Vol3DSimple<DisplaceFieldF3D>(szBlock_x, szBlock_y, szBlock_z);
					if (!pDFBlock)
					{
						printf("ERROR: Fail to allocate memory for pDFBlock.\n");
						//if (p_img_warp_4d) 		{ delete4dpointer(p_img_warp_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
						//if (p_img_sub_4d) 		{ delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
						if (pSubDF)				{ delete pSubDF;					pSubDF = 0; }
						return false;
					}
					DisplaceFieldF3D ***pppDFBlock = pDFBlock->getData3dHandle();

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
					//���������
					read_up_point = new long long[3];
					read_low_point = new long long[3];
					//�ĸ�����
					read_front_left = new long long[3];
					read_front_right = new long long[3];
					read_back_left = new long long[3];
					read_back_right = new long long[3];
					//����ת��
					//����
					low_point[0] = (sz_img_block_size[0] + x_offset)/2;
					low_point[1] = (sz_img_block_size[1] + y_offset)/2;
					low_point[2] = z_offset;

					if (!q_dfblcokinterp_bspline(pppSubDF, x_bsplinebasis, sz_gridwnd, (pSubDF->sz0() - 4) / 2, (pSubDF->sz1() - 4) / 2, 0, pppDFBlock))//��һ��С��
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}
					read_low_point[0] = (floor(low_point[0] + pppDFBlock[0][1][1].sx)) > 0 ? floor(low_point[0] + pppDFBlock[0][1][1].sx) : 0;
					read_low_point[1] = (floor(low_point[1] + pppDFBlock[0][1][1].sy)) > 0 ? floor(low_point[1] + pppDFBlock[0][1][1].sy) : 0;
					read_low_point[2] = (floor(low_point[2] + pppDFBlock[0][1][1].sz)) > 0 ? floor(low_point[2] + pppDFBlock[0][1][1].sz) : 0;
					//����
					up_point[0] = (sz_img_block_size[0] + x_offset) / 2;
					up_point[1] = (sz_img_block_size[1] + y_offset) / 2;
					up_point[2] = (sz_img_block_size[2] + z_offset) / 2;
					
					if (!q_dfblcokinterp_bspline(pppSubDF, x_bsplinebasis, sz_gridwnd, (pSubDF->sz0() - 4) / 2, (pSubDF->sz1() - 4) / 2, (pSubDF->sz2() - 4) / 2, pppDFBlock))//��һ��С��
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}
					read_up_point[0] = (floor(up_point[0] + pppDFBlock[3][1][1].sx)) > 0 ? floor(up_point[0] + pppDFBlock[3][1][1].sx) : 0;
					read_up_point[1] = (floor(up_point[1] + pppDFBlock[3][1][1].sy)) > 0 ? floor(up_point[1] + pppDFBlock[3][1][1].sy) : 0;
					read_up_point[2] = (floor(up_point[2] + pppDFBlock[3][1][1].sz)) > 0 ? floor(up_point[2] + pppDFBlock[3][1][1].sz) : 0;
					//ǰ��
					front_left[0] = x_offset;
					front_left[1] = (sz_img_block_size[1] + y_offset) / 2;
					front_left[2] = (sz_img_block_size[2] + z_offset) / 2;;

					if (!q_dfblcokinterp_bspline(pppSubDF, x_bsplinebasis, sz_gridwnd, 0, (pSubDF->sz1() - 4) / 2, (pSubDF->sz2() - 4) / 2, pppDFBlock))//��һ��С��
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}
					read_front_left[0] = (floor(front_left[0] + pppDFBlock[1][1][0].sx)) > 0 ? floor(front_left[0] + pppDFBlock[1][1][0].sx) : 0;
					read_front_left[1] = (floor(front_left[1] + pppDFBlock[1][1][0].sy)) > 0 ? floor(front_left[1] + pppDFBlock[1][1][0].sy) : 0;
					read_front_left[2] = (floor(front_left[2] + pppDFBlock[1][1][0].sz)) > 0 ? floor(front_left[2] + pppDFBlock[1][1][0].sz) : 0;
					//����
					back_right[0] = sz_img_block_size[0] + x_offset;
					back_right[1] = (sz_img_block_size[1] + y_offset) / 2;
					back_right[2] = (sz_img_block_size[2] + z_offset) / 2;;

					if (!q_dfblcokinterp_bspline(pppSubDF, x_bsplinebasis, sz_gridwnd, pSubDF->sz0() - 4, (pSubDF->sz1() - 4) / 2, (pSubDF->sz2() - 4) / 2, pppDFBlock))//��һ��С��
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}
					read_back_right[0] = (floor(back_right[0] + pppDFBlock[1][1][3].sx)) > 0 ? floor(back_right[0] + pppDFBlock[1][1][3].sx) : 0;
					read_back_right[1] = (floor(back_right[1] + pppDFBlock[1][1][3].sy)) > 0 ? floor(back_right[1] + pppDFBlock[1][1][3].sy) : 0;
					read_back_right[2] = (floor(back_right[2] + pppDFBlock[1][1][3].sz)) > 0 ? floor(back_right[2] + pppDFBlock[1][1][3].sz) : 0;
					//ǰ��
					front_right[0] = (sz_img_block_size[0] + x_offset) / 2;
					front_right[1] = y_offset;
					front_right[2] = (sz_img_block_size[2] + z_offset) / 2;;

					if (!q_dfblcokinterp_bspline(pppSubDF, x_bsplinebasis, sz_gridwnd, (pSubDF->sz0() - 4) / 2, 0, (pSubDF->sz2() - 4) / 2, pppDFBlock))//��һ��С��
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}
					read_front_right[0] = (floor(front_right[0] + pppDFBlock[1][0][1].sx)) > 0 ? floor(front_right[0] + pppDFBlock[1][0][1].sx) : 0;
					read_front_right[1] = (floor(front_right[1] + pppDFBlock[1][0][1].sy)) > 0 ? floor(front_right[1] + pppDFBlock[1][0][1].sy) : 0;
					read_front_right[2] = (floor(front_right[2] + pppDFBlock[1][0][1].sz)) > 0 ? floor(front_right[2] + pppDFBlock[1][0][1].sz) : 0;
					//����
					back_left[0] = (sz_img_block_size[0] + x_offset) / 2;
					back_left[1] = sz_img_block_size[1] + y_offset;
					back_left[2] = (sz_img_block_size[2] + z_offset) / 2;;

					if (!q_dfblcokinterp_bspline(pppSubDF, x_bsplinebasis, sz_gridwnd, (pSubDF->sz0() - 4) / 2, pSubDF->sz1() - 4, (pSubDF->sz2() - 4) / 2, pppDFBlock))//��һ��С��
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}
					read_back_left[0] = (floor(back_left[0] + pppDFBlock[1][3][1].sx)) > 0 ? floor(back_left[0] + pppDFBlock[1][3][1].sx) : 0;
					read_back_left[1] = (floor(back_left[1] + pppDFBlock[1][3][1].sy)) > 0 ? floor(back_left[1] + pppDFBlock[1][3][1].sy) : 0;
					read_back_left[2] = (floor(back_left[2] + pppDFBlock[1][3][1].sz)) > 0 ? floor(back_left[2] + pppDFBlock[1][3][1].sz) : 0;




					//��ʼ�㡪��������������������������������������������������������������������������������������������������������������
					start_point[0] = x_offset;
					start_point[1] = y_offset;
					start_point[2] = z_offset;

					if (!q_dfblcokinterp_bspline(pppSubDF, x_bsplinebasis, sz_gridwnd, 0, 0, 0, pppDFBlock))//��һ��С��
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}
					read_start_point[0] = (floor(start_point[0] + pppDFBlock[0][0][0].sx)) > 0 ? floor(start_point[0] + pppDFBlock[0][0][0].sx) : 0;
					read_start_point[1] = (floor(start_point[1] + pppDFBlock[0][0][0].sy)) > 0 ? floor(start_point[1] + pppDFBlock[0][0][0].sy) : 0;
					read_start_point[2] = (floor(start_point[2] + pppDFBlock[0][0][0].sz)) > 0 ? floor(start_point[2] + pppDFBlock[0][0][0].sz) : 0;

					//low_right
					low_right_point[0] = sz_img_block_size[0] + x_offset;
					low_right_point[1] = y_offset;
					low_right_point[2] = z_offset;

					if (!q_dfblcokinterp_bspline(pppSubDF, x_bsplinebasis, sz_gridwnd, pSubDF->sz0() - 4, 0, 0, pppDFBlock))//���һ��С��
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}

					read_low_right_point[0] = (floor(low_right_point[0] + pppDFBlock[0][0][3].sx)) > 0 ? floor(low_right_point[0] + pppDFBlock[0][0][3].sx) : 0;
					read_low_right_point[1] = (floor(low_right_point[1] + pppDFBlock[0][0][3].sy)) > 0 ? floor(low_right_point[1] + pppDFBlock[0][0][3].sy) : 0;
					read_low_right_point[2] = (floor(low_right_point[2] + pppDFBlock[0][0][3].sz)) > 0 ? floor(low_right_point[2] + pppDFBlock[0][0][3].sz) : 0;

					//low_left
					low_left_point[0] = x_offset;
					low_left_point[1] = sz_img_block_size[1] + y_offset;
					low_left_point[2] = z_offset;

					if (!q_dfblcokinterp_bspline(pppSubDF, x_bsplinebasis, sz_gridwnd, 0, pSubDF->sz1() - 4, 0, pppDFBlock))//���һ��С��
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}

					read_low_left_point[0] = (floor(low_left_point[0] + pppDFBlock[0][3][0].sx)) > 0 ? floor(low_left_point[0] + pppDFBlock[0][3][0].sx) : 0;
					read_low_left_point[1] = (floor(low_left_point[1] + pppDFBlock[0][3][0].sy)) > 0 ? floor(low_left_point[1] + pppDFBlock[0][3][0].sy) : 0;
					read_low_left_point[2] = (floor(low_left_point[2] + pppDFBlock[0][3][0].sz)) > 0 ? floor(low_left_point[2] + pppDFBlock[0][3][0].sz) : 0;

					//low_end
					low_end_point[0] = sz_img_block_size[0] + x_offset;
					low_end_point[1] = sz_img_block_size[1] + y_offset;
					low_end_point[2] = z_offset;

					if (!q_dfblcokinterp_bspline(pppSubDF, x_bsplinebasis, sz_gridwnd, pSubDF->sz0() - 4, pSubDF->sz1() - 4, 0, pppDFBlock))//���һ��С��
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}

					read_low_end_point[0] = (floor(low_end_point[0] + pppDFBlock[0][3][3].sx)) > 0 ? floor(low_end_point[0] + pppDFBlock[0][3][3].sx) : 0;
					read_low_end_point[1] = (floor(low_end_point[1] + pppDFBlock[0][3][3].sy)) > 0 ? floor(low_end_point[1] + pppDFBlock[0][3][3].sy) : 0;
					read_low_end_point[2] = (floor(low_end_point[2] + pppDFBlock[0][3][3].sz)) > 0 ? floor(low_end_point[2] + pppDFBlock[0][3][3].sz) : 0;


					//high_start
					high_start_point[0] = x_offset;
					high_start_point[1] = y_offset;
					high_start_point[2] = sz_img_block_size[2] + z_offset;

					if (!q_dfblcokinterp_bspline(pppSubDF, x_bsplinebasis, sz_gridwnd, 0, 0, pSubDF->sz2() - 4, pppDFBlock))//���һ��С��
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}

					read_high_start_point[0] = (ceil(high_start_point[0] + pppDFBlock[3][0][0].sx)) > 0 ? ceil(high_start_point[0] + pppDFBlock[3][0][0].sx) : 0;
					read_high_start_point[1] = (ceil(high_start_point[1] + pppDFBlock[3][0][0].sy)) > 0 ? ceil(high_start_point[1] + pppDFBlock[3][0][0].sy) : 0;
					read_high_start_point[2] = (ceil(high_start_point[2] + pppDFBlock[3][0][0].sz)) > 0 ? ceil(high_start_point[2] + pppDFBlock[3][0][0].sz) : 0;

					//high_right
					high_right_point[0] = sz_img_block_size[0] + x_offset;
					high_right_point[1] = y_offset;
					high_right_point[2] = sz_img_block_size[2] + z_offset;

					if (!q_dfblcokinterp_bspline(pppSubDF, x_bsplinebasis, sz_gridwnd, pSubDF->sz0() - 4, 0, pSubDF->sz2() - 4, pppDFBlock))//���һ��С��
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}

					read_high_right_point[0] = (ceil(high_right_point[0] + pppDFBlock[3][0][3].sx)) > 0 ? ceil(high_right_point[0] + pppDFBlock[3][0][3].sx) : 0;
					read_high_right_point[1] = (ceil(high_right_point[1] + pppDFBlock[3][0][3].sy)) > 0 ? ceil(high_right_point[1] + pppDFBlock[3][0][3].sy) : 0;
					read_high_right_point[2] = (ceil(high_right_point[2] + pppDFBlock[3][0][3].sz)) > 0 ? ceil(high_right_point[2] + pppDFBlock[3][0][3].sz) : 0;

					//high_left
					high_left_point[0] = x_offset;
					high_left_point[1] = sz_img_block_size[1] + y_offset;
					high_left_point[2] = sz_img_block_size[2] + z_offset;

					if (!q_dfblcokinterp_bspline(pppSubDF, x_bsplinebasis, sz_gridwnd, 0, pSubDF->sz1() - 4, pSubDF->sz2() - 4, pppDFBlock))//���һ��С��
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}

					read_high_left_point[0] = (ceil(high_left_point[0] + pppDFBlock[3][3][0].sx)) > 0 ? ceil(high_left_point[0] + pppDFBlock[3][3][0].sx) : 0;
					read_high_left_point[1] = (ceil(high_left_point[1] + pppDFBlock[3][3][0].sy)) > 0 ? ceil(high_left_point[1] + pppDFBlock[3][3][0].sy) : 0;
					read_high_left_point[2] = (ceil(high_left_point[2] + pppDFBlock[3][3][0].sz)) > 0 ? ceil(high_left_point[2] + pppDFBlock[3][3][0].sz) : 0;

					//�յ㡪��������������������������������������������������������������������������������������������������������������
					end_point[0] = sz_img_block_size[0] + x_offset;
					end_point[1] = sz_img_block_size[1] + y_offset;
					end_point[2] = sz_img_block_size[2] + z_offset;

					if (!q_dfblcokinterp_bspline(pppSubDF, x_bsplinebasis, sz_gridwnd, pSubDF->sz0() - 4, pSubDF->sz1() - 4, pSubDF->sz2() - 4, pppDFBlock))//���һ��С��
					{
						printf("ERROR: Fail to ues conversion_point.\n");
						//if (read_centre_point){ delete[]read_centre_point; read_centre_point = 0; }
						if (read_start_point){ delete[]read_start_point; read_start_point = 0; }
						if (read_end_point){ delete[]read_end_point; read_end_point = 0; }

						return false;
					}

					read_end_point[0] = (ceil(end_point[0] + pppDFBlock[3][3][3].sx)) > 0 ? ceil(end_point[0] + pppDFBlock[3][3][3].sx) : 0;
					read_end_point[1] = (ceil(end_point[1] + pppDFBlock[3][3][3].sy)) > 0 ? ceil(end_point[1] + pppDFBlock[3][3][3].sy) : 0;
					read_end_point[2] = (ceil(end_point[2] + pppDFBlock[3][3][3].sz)) > 0 ? ceil(end_point[2] + pppDFBlock[3][3][3].sz) : 0;

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
//
//					printf("��ȡ��low�����꣺%d,%d,%d\n", read_low_point[0], read_low_point[1], read_low_point[2]);
//					printf("��ȡ��up�����꣺%d,%d,%d\n", read_up_point[0], read_up_point[1], read_up_point[2]);
//					printf("��ȡ��front_left�����꣺%d,%d,%d\n", read_front_left[0], read_front_left[1], read_front_left[2]);
//					printf("��ȡ��front_right�����꣺%d,%d,%d\n", read_front_right[0], read_front_right[1], read_front_right[2]);
//					printf("��ȡ��back_left�����꣺%d,%d,%d\n", read_back_left[0], read_back_left[1], read_back_left[2]);
//					printf("��ȡ��back_right�����꣺%d,%d,%d\n", read_back_right[0], read_back_right[1], read_back_right[2]);
					//����
					sort_point_x.push_back(read_start_point[0]);
					sort_point_x.push_back(read_low_right_point[0]);
					sort_point_x.push_back(read_low_left_point[0]);
					sort_point_x.push_back(read_low_end_point[0]);
					sort_point_x.push_back(read_high_start_point[0]);
					sort_point_x.push_back(read_high_right_point[0]);
					sort_point_x.push_back(read_high_left_point[0]);
					sort_point_x.push_back(read_end_point[0]);
					sort_point_x.push_back(read_low_point[0]);
					sort_point_x.push_back(read_up_point[0]);
					sort_point_x.push_back(read_front_left[0]);
					sort_point_x.push_back(read_front_right[0]);
					sort_point_x.push_back(read_back_left[0]);
					sort_point_x.push_back(read_back_right[0]);

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
					sort_point_y.push_back(read_low_point[1]);
					sort_point_y.push_back(read_up_point[1]);
					sort_point_y.push_back(read_front_left[1]);
					sort_point_y.push_back(read_front_right[1]);
					sort_point_y.push_back(read_back_left[1]);
					sort_point_y.push_back(read_back_right[1]);

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
					sort_point_z.push_back(read_low_point[2]);
					sort_point_z.push_back(read_up_point[2]);
					sort_point_z.push_back(read_front_left[2]);
					sort_point_z.push_back(read_front_right[2]);
					sort_point_z.push_back(read_back_left[2]);
					sort_point_z.push_back(read_back_right[2]);

					sort(sort_point_z.begin(), sort_point_z.end());
//					printf("z����С��%d\n", sort_point_z.front());
//					printf("z�����%d\n", sort_point_z.back());

					//��ȡ��ƫ�Ƽ���

					sz_img_read_block = new long long[4];
					//sz_img_read_block[0] = (x_offset - read_offset) < 0 ? 0 : (x_offset - read_offset);
					//sz_img_read_block[1] = (y_offset - read_offset) < 0 ? 0 : (y_offset - read_offset);
					//sz_img_read_block[2] = (z_offset - read_offset) < 0 ? 0 : (z_offset - read_offset);//50Ϊ������
					
					//int off_para = 1;//����С�ĵ��ƫ��������һ��
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
						if (p_img_sub2tar_4d) 	{ delete4dpointer(p_img_sub2tar_4d, sz_img_block_size[0], sz_img_block_size[1], sz_img_block_size[2], sz_img_block_size[3]); }
						if (p_img_sub2tar_affine) 		{ delete[]p_img_sub2tar_affine;		p_img_sub2tar_affine = 0; }
						if (pDFBlock)			{ delete pDFBlock;			pDFBlock = 0; }
						if (pSubDF)				{ delete pSubDF;				pSubDF = 0; }
						
						return false;
					}
					//printf("mark1");
					read_time = omp_get_wtime();
					if (!read_2Dtif_BLOCK(path_TIF, p_img_read_block, p_img_read_4d_block, sz_img_read_block, sz_img_ori, sz_img_read_size, per_sz_block, sz_img_read_end_block))//sz_img_read_size�Ƕ�ȡ��ĳߴ�
					{
						printf("ERROR: Fail to read_2Dtif_BLOCK.\n");
						if (p_img_sub2tar_4d) 	{ delete4dpointer(p_img_sub2tar_4d, sz_img_block_size[0], sz_img_block_size[1], sz_img_block_size[2], sz_img_block_size[3]); }
						if (p_img_sub2tar_affine) 		{ delete[]p_img_sub2tar_affine;		p_img_sub2tar_affine = 0; }
						if (p_img_read_4d_block) 	{ delete4dpointer(p_img_read_4d_block, sz_img_read_size[0], sz_img_read_size[1], sz_img_read_size[2], sz_img_read_size[3]); }
						if (p_img_read_block) 		{ delete[]p_img_read_block;		p_img_read_block = 0; }
						if (pDFBlock)			{ delete pDFBlock;			pDFBlock = 0; }
						if (pSubDF)				{ delete pSubDF;				pSubDF = 0; }
						return false;
					}//sz_img_read_block�Ƕ�ȡ���ƫ��
					read_sum += (float)(omp_get_wtime() - read_time);//readʱ�����
					printf("\t>>read_time time consume %.2f s\n", (float)(omp_get_wtime() - read_time));
					//printf("mark2\n");
					printf("warp�Ŀ�����꣺%d;%d;%d\n", x_offset, y_offset, z_offset);//warp�Ŀ������
					printf("warp�Ŀ����ʵ�ߴ磺%d;%d;%d\n", sz_img_block_size[0], sz_img_block_size[1], sz_img_block_size[2]);//warp�Ŀ�ĳߴ�
					//printf("warp�Ŀ������ߴ磺%d;%d;%d\n", sz_img_block_size_fake[0], sz_img_block_size_fake[1], sz_img_block_size_fake[2]);//warp�Ŀ�ĳߴ�
					printf("��ȡ�Ŀ�����꣺%d;%d;%d\n", sz_img_read_block[0], sz_img_read_block[1], sz_img_read_block[2]);//��ȡ�Ŀ������
					printf("��ȡ�Ŀ�ĳߴ磺%d;%d;%d\n", sz_img_read_size[0], sz_img_read_size[1], sz_img_read_size[2]);//��ȡ�Ŀ�ĳߴ�
					//ĿǰSTPS�Ƿ���warp�߶�GPU���Ʋ�������Ȱ�����ע��;
					/*
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
					*/
					aff_warp = clock();

					if (!STPS_interpolate_CYF(GPU_type, pSubDF, pppSubDF, x_bsplinebasis, sz_gridwnd, pppDFBlock, p_img_read_4d_block, sz_img_block_size, szBlock_x, szBlock_y,
						szBlock_z, i_interpmethod_img, p_img_sub2tar_4d, sz_img_read_size, p_img_read_block, p_img_sub2tar_affine, x_offset, y_offset, z_offset, sz_img_ori, sz_img_read_block))
					{
						printf("ERROR: chazhi is wrong.\n");
						if (p_img_sub2tar_4d) 	{ delete4dpointer(p_img_sub2tar_4d, sz_img_block_size[0], sz_img_block_size[1], sz_img_block_size[2], sz_img_block_size[3]); }
						if (p_img_sub2tar_affine) 		{ delete[]p_img_sub2tar_affine;		p_img_sub2tar_affine = 0; }
						if (p_img_read_4d_block) 	{ delete4dpointer(p_img_read_4d_block, sz_img_read_size[0], sz_img_read_size[1], sz_img_read_size[2], sz_img_read_size[3]); }
						if (p_img_read_block) 		{ delete[]p_img_read_block;		p_img_read_block = 0; }
						if (pDFBlock)			{ delete pDFBlock;			pDFBlock = 0; }
						if (pSubDF)				{ delete pSubDF;				pSubDF = 0; }
						return false;
					}
					stps_sum += (float)(clock() - aff_warp) / CLOCKS_PER_SEC;
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
					base_block_path << TeraWarp_save_base << "_stps_block" << "/" << abs_pos_z.str() << "_" << abs_pos_y.str() << "_" << abs_pos_x.str() << ".v3draw";
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
					if (read_up_point){ delete[]read_up_point; read_up_point = 0; }
					if (read_low_point){ delete[]read_low_point; read_low_point = 0; }
					if (read_front_left){ delete[]read_front_left; read_front_left = 0; }
					if (read_front_right){ delete[]read_front_right; read_front_right = 0; }
					if (read_back_left){ delete[]read_back_left; read_back_left = 0; }
					if (read_back_right){ delete[]read_back_right; read_back_right = 0; }
					if (sz_img_read_block){ delete[]sz_img_read_block; sz_img_read_block = 0; }
					if (sz_img_read_size){ delete[]sz_img_read_size; sz_img_read_size = 0; }
					if (pDFBlock)			{ delete pDFBlock;			pDFBlock = 0; }
					if (pSubDF)				{ delete pSubDF;				pSubDF = 0; }
				}
			}
	}
	//long long sz_img_resize[4] = { 568, 320, 456, 1 };
	//std::string imgSrcFile_save = "D:/BIG_Warp_File/OUTPUT/2D_Slice-18052_1000_0920";
	//std::string imgSrcFile_save = "E:/STPS/Outpu_t/2d_raw_slice_400_18052_x2_0906";
	
	//std::string imgSrcFile_save = "F:/BIG_WARP_FILE/OUTPUT/18052_XYx4_Zx8_2D_1000_Slice_0914";
	//std::string imgSrcFile_save = "E:/STPS/Outpu_t/2d_raw_Slice_800_18052_18_0914";
	//char  *imgSrcFile_ori = "E:/STPS/Outpu_t/2d_raw_400_18052_0907";
	/*//********************************
	//ƴ�Ӳ�������
	
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
	
	//********************************
	*/
	
	PJ_time = omp_get_wtime();
	if (type_PJ)
	{
		//q_imagewarp_stitch_2DTIF_new(block_path.c_str(), imgSrcFile_save, sz_img_ori, per_block_siza, sz_offset_x, sz_offset_y);
		q_imagewarp_stitch_2DRAW_new(block_path.str().c_str(), Slice_path.str().c_str(), sz_img_ori, per_block_siza, sz_offset_x, sz_offset_y, Z_Slice_control);
	}
	printf("\t>>read all time consume %.2f s\n", read_sum);
	printf("\t>>BSP time consume: %.2f s\n", bsp_time);
	printf("\t>>stps all time consume %.2f s\n", stps_sum);
	printf("\t>>save all time consume %.2f s\n", save_sum);
	printf("\t>>PJ_time time consume %.2f s\n", (float)(omp_get_wtime() - PJ_time));


//	if (path_TIF){ delete[]path_TIF; }
    printf("6. free memory. \n");
	if (sz_offset_x){ delete[]sz_offset_x; }
	if (sz_offset_y){ delete[]sz_offset_y; }
    if (cpt_subject) { delete cpt_subject; cpt_subject = 0; }


	free(H_X);
	free(H_Y);
	free(H_Z);



	block_path.clear();
	block_path.str("");
	Slice_path.clear();
	Slice_path.str("");


	//------------------------------------------------------�ֿ�warp------------------------------------------------------------------------------
	return true;
}


bool imgwarp_smallmemory_image(const QList<ImageMarker> &ql_marker_tar, const QList<ImageMarker> &ql_marker_sub,
                               V3DLONG szBlock_x, V3DLONG szBlock_y, V3DLONG szBlock_z, int i_interpmethod_df, int i_interpmethod_img,
                               unsigned char *p_img_sub, V3DLONG *sz_img_sub, string TeraWarp_save_base, int GPU_type)
{
    if (p_img_sub == 0 || sz_img_sub == 0)
    {
        printf("ERROR: p_img_sub or sz_img_sub is invalid.\n");
        return false;
    }
    if (ql_marker_tar.size() == 0 || ql_marker_sub.size() == 0 || ql_marker_tar.size() != ql_marker_sub.size())
    {
        printf("ERROR: target or subject control points is invalid!\n");
        return false;
    }
    if (szBlock_x <= 0 || szBlock_y <= 0 || szBlock_z <= 0)
    {
        printf("ERROR: block size is invalid!\n");
        return false;
    }
    if (szBlock_x >= sz_img_sub[0] || szBlock_y >= sz_img_sub[1] || szBlock_z >= sz_img_sub[2])
    {
        printf("ERROR: block size should smaller than the image size!\n");
        return false;
    }
    if (i_interpmethod_df != 0 && i_interpmethod_df != 1)
    {
        printf("ERROR: DF_interp_method should be 0(linear) or 1(bspline)!\n");
        return false;
    }
    if (i_interpmethod_img != 0 && i_interpmethod_img != 1)
    {
        printf("ERROR: img_interp_method should be 0(linear) or 1(nn)!\n");
        return false;
    }
    if (i_interpmethod_df == 1 && (szBlock_x != szBlock_y || szBlock_x != szBlock_z))
    {
        printf("ERROR: df_interp_method=bspline need szBlock_x=szBlock_y=szBlock_z!\n");
        return false;
    }

    vector<Coord3D_JBA> matchTargetPos, matchSubjectPos;

    for (V3DLONG i = 0; i<ql_marker_tar.size(); i++)
    {
        Coord3D_JBA tmpc;
        tmpc.x = 1 * ql_marker_tar.at(i).x;	tmpc.y = 1 * ql_marker_tar.at(i).y;	tmpc.z = 1 * ql_marker_tar.at(i).z;
        matchTargetPos.push_back(tmpc);
        tmpc.x = 1 * ql_marker_sub.at(i).x;	tmpc.y = 1 * ql_marker_sub.at(i).y;	tmpc.z = 1 * ql_marker_sub.at(i).z;
        matchSubjectPos.push_back(tmpc);
    }

    int nCpt = matchTargetPos.size();
    Image2DSimple<MYFLOAT_JBA> * cpt_subject = 0;
    Matrix xnx4_c(nCpt, 4);
    Matrix x4x4_d(4, 4);
    float *H_X = 0;
    float *H_Y = 0;
    float *H_Z = 0;
    clock_t BSP;
    BSP = clock();
//    if (!(compute_df_stps_subsampled_volume_4bspline_per(matchTargetPos, matchSubjectPos, x4x4_d, xnx4_c, H_X, H_Y, H_Z, nCpt, cpt_subject, GPU_type)))
//    {
//        printf("ERROR:compute_df_stps_subsampled_volume_4bspline_per() return false. \n");
//        return false;
//    }

    Vol3DSimple<DisplaceFieldF3D> *pSubDF=0;


    pSubDF = compute_df_stps_subsampled_volume_4bspline(matchTargetPos, matchSubjectPos, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], szBlock_x, szBlock_y, szBlock_z, GPU_type);


    if(!pSubDF)
    {
        printf("Fail to produce the subsampled DF.\n");
        return false;
    }
    DisplaceFieldF3D ***pppSubDF=pSubDF->getData3dHandle();

    unsigned char *p_img_warp = new unsigned char[sz_img_sub[0] * sz_img_sub[1] * sz_img_sub[2] * sz_img_sub[3]]();
    if(!p_img_warp)
    {
        printf("ERROR: Fail to allocate memory for p_img_warp.\n");
        if(pSubDF)				{delete pSubDF;					pSubDF=0;}
        return false;
    }
    unsigned char ****p_img_warp_4d = 0, ****p_img_sub_4d = 0;
    if(!new4dpointer(p_img_warp_4d,sz_img_sub[0],sz_img_sub[1],sz_img_sub[2],sz_img_sub[3],p_img_warp) ||
       !new4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3], p_img_sub))
    {
        printf("ERROR: Fail to allocate memory for the 4d pointer of image.\n");
        if(p_img_warp_4d) 		{delete4dpointer(p_img_warp_4d,sz_img_sub[0],sz_img_sub[1],sz_img_sub[2],sz_img_sub[3]);}
        if (p_img_sub_4d) 		{ delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
        if(p_img_warp) 			{delete []p_img_warp;			p_img_warp=0;}
        if(pSubDF)				{delete pSubDF;					pSubDF=0;}
        return false;
    }

    Vol3DSimple<DisplaceFieldF3D> *pDFBlock=new Vol3DSimple<DisplaceFieldF3D> (szBlock_x,szBlock_y,szBlock_z);
    if(!pDFBlock)
    {
        printf("ERROR: Fail to allocate memory for pDFBlock.\n");
        if(p_img_warp_4d) 		{delete4dpointer(p_img_warp_4d,sz_img_sub[0],sz_img_sub[1],sz_img_sub[2],sz_img_sub[3]);}
        if (p_img_sub_4d) 		{ delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
        if(p_img_warp) 			{delete []p_img_warp;			p_img_warp=0;}
        if(pSubDF)				{delete pSubDF;					pSubDF=0;}
        return false;
    }
    DisplaceFieldF3D ***pppDFBlock=pDFBlock->getData3dHandle();

    //initialize the bspline basis function
    V3DLONG sz_gridwnd=szBlock_x;
    Matrix x_bsplinebasis(pow(double(sz_gridwnd),3.0),pow(4.0,3.0));
    if(!q_nonrigid_ini_bsplinebasis_3D(sz_gridwnd,x_bsplinebasis))
    {
        printf("ERROR: q_ini_bsplinebasis_3D() return false!\n");
        if(p_img_warp_4d) 		{delete4dpointer(p_img_warp_4d,sz_img_sub[0],sz_img_sub[1],sz_img_sub[2],sz_img_sub[3]);}
        if (p_img_sub_4d) 		{ delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
        if(p_img_warp) 			{delete []p_img_warp;			p_img_warp=0;}
        if(pSubDF)				{delete pSubDF;					pSubDF=0;}
        return false;
    }
    printf("\t>>x_bsplinebasis:[%d,%d]\n",x_bsplinebasis.nrows(),x_bsplinebasis.ncols());

    if(GPU_type==0)	//linear interpolate the SubDfBlock to DFBlock and do warp block by block
    {
        for(V3DLONG substart_z=0;substart_z<pSubDF->sz2()-1;substart_z++)
            for(V3DLONG substart_y=0;substart_y<pSubDF->sz1()-1;substart_y++)
                for(V3DLONG substart_x=0;substart_x<pSubDF->sz0()-1;substart_x++)
                {
                    //bspline interpolate the SubDfBlock to DFBlock
                    q_dfblcokinterp_bspline(pppSubDF, x_bsplinebasis, sz_gridwnd, substart_x, substart_y, substart_z, pppDFBlock);
                    //warp image block using DFBlock
                    q_imgblockwarp(p_img_sub_4d,sz_img_sub,pppDFBlock,szBlock_x,szBlock_y,szBlock_z,i_interpmethod_img,substart_x,substart_y,substart_z,p_img_warp_4d, sz_img_sub);

                }
    }
    else						//bspline interpolate the SubDfBlock to DFBlock and do warp block by block
    {

        int gsz2 = pSubDF->sz2();
        int gsz1 = pSubDF->sz1();
        int gsz0 = pSubDF->sz0();

        clock_t stps_interpolation;
        stps_interpolation = clock();
        gpu_interpolation_stps(gsz2, gsz1, gsz0, pppSubDF, x_bsplinebasis, sz_gridwnd, pppDFBlock, p_img_sub_4d, sz_img_sub, szBlock_x, szBlock_y, szBlock_z, i_interpmethod_img, p_img_warp_4d, sz_img_sub, p_img_sub, p_img_warp);


        printf("\t>>interpolation time consume %.2f s\n", (float)(clock() - stps_interpolation) / CLOCKS_PER_SEC);

    }

    V3DLONG sz_img_sub_[4] = {sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]};
    saveImage(TeraWarp_save_base.c_str(), p_img_warp, sz_img_sub_, 1);

    //free memory
    if(p_img_warp_4d) 		{delete4dpointer(p_img_warp_4d,sz_img_sub[0],sz_img_sub[1],sz_img_sub[2],sz_img_sub[3]);}
    if (p_img_sub_4d) 		{delete4dpointer(p_img_sub_4d, sz_img_sub[0],sz_img_sub[1],sz_img_sub[2],sz_img_sub[3]); }
    if(pDFBlock)			{delete pDFBlock;			pDFBlock=0;}
    if(pSubDF)				{delete pSubDF;				pSubDF=0;}
    if(p_img_warp)          {delete[]p_img_warp; p_img_warp=0;}

    return true;

}



//bool STPS_interpolate_CYF(const int GPU_type, Vol3DSimple<DisplaceFieldF3D> *pSubDF, DisplaceFieldF3D ***pppSubDF, Matrix x_bsplinebasis,
//	const V3DLONG sz_gridwnd, DisplaceFieldF3D ***&pppDFBlock, unsigned char ****&p_img_sub_4d, const V3DLONG *sz_img_sub_, const V3DLONG szBlock_x,
//	const V3DLONG szBlock_y, const V3DLONG szBlock_z, const int i_interpmethod_img, unsigned char ****&p_img_warp_4d, const long long *sz_img_sub_ori_,
//	const unsigned char *p_img_sub, unsigned char *p_img_warp, long long x_offset, long long y_offset, long long z_offset, const long long *sz_img_ori,
//	const long long *sz_img_read_block)
//{
//    V3DLONG sz_img_sub[4] = {sz_img_sub_[0], sz_img_sub_[1], sz_img_sub_[2], sz_img_sub_[3]};
//    V3DLONG sz_img_sub_ori[4] = {sz_img_sub_ori_[0], sz_img_sub_ori_[1], sz_img_sub_ori_[2], sz_img_sub_ori_[3]};
//	if (GPU_type)
//	{
//		int gsz2 = pSubDF->sz2();
//		int gsz1 = pSubDF->sz1();
//		int gsz0 = pSubDF->sz0();
//		clock_t stps_interpolation;
//		stps_interpolation = clock();
//
//		gpu_interpolation(gsz2, gsz1, gsz0, pppSubDF, x_bsplinebasis, sz_gridwnd, pppDFBlock, p_img_sub_4d, sz_img_sub, szBlock_x, szBlock_y, szBlock_z, i_interpmethod_img,
//			p_img_warp_4d, sz_img_sub_ori, p_img_sub, p_img_warp, x_offset, y_offset, z_offset, sz_img_read_block[0], sz_img_read_block[1], sz_img_read_block[2],sz_img_ori[2],
//			sz_img_ori[1],sz_img_ori[0]);
//
//
//		printf("\t>>interpolation time consume %.2f s\n", (float)(clock() - stps_interpolation) / CLOCKS_PER_SEC);
//	}
//	else
//	{
//		for (V3DLONG substart_z = 0; substart_z<pSubDF->sz2() - 1 - 2; substart_z++)
//			for (V3DLONG substart_y = 0; substart_y<pSubDF->sz1() - 1 - 2; substart_y++)
//				for (V3DLONG substart_x = 0; substart_x<pSubDF->sz0() - 1 - 2; substart_x++)
//				{
//					//bspline interpolate the SubDfBlock to DFBlock
//					q_dfblcokinterp_bspline(pppSubDF, x_bsplinebasis, sz_gridwnd, substart_x, substart_y, substart_z, pppDFBlock);
//					//warp image block using DFBlock
//					q_imgblockwarp_CYF(p_img_sub_4d, sz_img_sub, pppDFBlock, szBlock_x, szBlock_y, szBlock_z, i_interpmethod_img, substart_x, substart_y, substart_z, p_img_warp_4d,
//						sz_img_sub_ori, x_offset, y_offset, z_offset, sz_img_ori, sz_img_read_block);
//				}
//	}
//	return true;
//}




bool q_imgblockwarp_CYF(unsigned char ****&p_img_sub_4d, const V3DLONG *sz_img_sub, DisplaceFieldF3D ***&pppDFBlock,
	const V3DLONG szBlock_x, const V3DLONG szBlock_y, const V3DLONG szBlock_z, const int i_interpmethod_img,
	const V3DLONG substart_x, const V3DLONG substart_y, const V3DLONG substart_z,unsigned char ****&p_img_warp_4d, 
	const V3DLONG *sz_img_sub_read, long long x_offset, long long y_offset, long long z_offset, const long long *sz_img_ori,
	const long long *sz_img_read_block)
{
	V3DLONG start_x, start_y, start_z;
	start_x = substart_x*szBlock_x;
	start_y = substart_y*szBlock_y;
	start_z = substart_z*szBlock_z;
	for (V3DLONG z = 0; z<szBlock_z; z++)
		for (V3DLONG y = 0; y<szBlock_y; y++)
			for (V3DLONG x = 0; x<szBlock_x; x++)
			{
				V3DLONG pos_warp[3];
				pos_warp[0] = start_x + x;
				pos_warp[1] = start_y + y;
				pos_warp[2] = start_z + z;
				if (pos_warp[0] >= sz_img_sub[0] || pos_warp[1] >= sz_img_sub[1] || pos_warp[2] >= sz_img_sub[2])//sz_img_sub����warp��Ĵ�С
					continue;

				double pos_sub[3], pos_sub_l[3];
				pos_sub[0] = pos_warp[0] + pppDFBlock[z][y][x].sx + x_offset - sz_img_read_block[0];
				pos_sub[1] = pos_warp[1] + pppDFBlock[z][y][x].sy + y_offset - sz_img_read_block[1];
				pos_sub[2] = pos_warp[2] + pppDFBlock[z][y][x].sz + z_offset - sz_img_read_block[2];//sz_img_read_block���Ƕ�ȡ�Ŀ�Ĵ�С

				pos_sub_l[0] = pos_warp[0] + pppDFBlock[z][y][x].sx + x_offset;
				pos_sub_l[1] = pos_warp[1] + pppDFBlock[z][y][x].sy + y_offset;
				pos_sub_l[2] = pos_warp[2] + pppDFBlock[z][y][x].sz + z_offset;

				if (pos_sub_l[0]<0 || pos_sub_l[0]>sz_img_ori[0] - 1 ||
					pos_sub_l[1]<0 || pos_sub_l[1]>sz_img_ori[1] - 1 ||
					pos_sub_l[2]<0 || pos_sub_l[2]>sz_img_ori[2] - 1)//sz_img_ori����ԭʼ��ͼ��Ĵ�С
				{
					for (V3DLONG c = 0; c<sz_img_sub[3]; c++)
						p_img_warp_4d[c][pos_warp[2]][pos_warp[1]][pos_warp[0]] = 0;
					continue;
				}
				if (pos_sub[0]<0 || pos_sub[0]>sz_img_sub_read[0] - 1 ||
					pos_sub[1]<0 || pos_sub[1]>sz_img_sub_read[1] - 1 ||
					pos_sub[2]<0 || pos_sub[2]>sz_img_sub_read[2] - 1)
				{
					for (V3DLONG c = 0; c<sz_img_sub[3]; c++)
						p_img_warp_4d[c][pos_warp[2]][pos_warp[1]][pos_warp[0]] = 0;
					continue;
				}

				//nearest neighbor interpolate
				if (i_interpmethod_img == 1)
				{
					V3DLONG pos_sub_nn[3];
					for (int i = 0; i<3; i++)
					{
						pos_sub_nn[i] = pos_sub[i] + 0.5;
						pos_sub_nn[i] = pos_sub_nn[i]<sz_img_sub[i] ? pos_sub_nn[i] : sz_img_sub[i] - 1;
					}
					for (V3DLONG c = 0; c<sz_img_sub[3]; c++)
						p_img_warp_4d[c][pos_warp[2]][pos_warp[1]][pos_warp[0]] = p_img_sub_4d[c][pos_sub_nn[2]][pos_sub_nn[1]][pos_sub_nn[0]];
				}
				//linear interpolate
				else if (i_interpmethod_img == 0)
				{
					//find 8 neighor pixels boundary
					V3DLONG x_s, x_b, y_s, y_b, z_s, z_b;
					x_s = floor(pos_sub[0]);		x_b = ceil(pos_sub[0]);
					y_s = floor(pos_sub[1]);		y_b = ceil(pos_sub[1]);
					z_s = floor(pos_sub[2]);		z_b = ceil(pos_sub[2]);

					//compute weight for left and right, top and bottom -- 4 neighbor pixel's weight in a slice
					double l_w, r_w, t_w, b_w;
					l_w = 1.0 - (pos_sub[0] - x_s);	r_w = 1.0 - l_w;
					t_w = 1.0 - (pos_sub[1] - y_s);	b_w = 1.0 - t_w;
					//compute weight for higer slice and lower slice
					double u_w, d_w;
					u_w = 1.0 - (pos_sub[2] - z_s);	d_w = 1.0 - u_w;

					//linear interpolate each channel
					for (V3DLONG c = 0; c<sz_img_sub[3]; c++)
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
						double intval = (u_w*higher_slice + d_w*lower_slice + 0.5);
						p_img_warp_4d[c][pos_warp[2]][pos_warp[1]][pos_warp[0]] = intval;
					}
				}

			}

	return true;
}





bool STPS_interpolate_sort_CYF(const int GPU_type, Vol3DSimple<DisplaceFieldF3D> *pSubDF, DisplaceFieldF3D ***pppSubDF, Matrix x_bsplinebasis,
    const V3DLONG sz_gridwnd, DisplaceFieldF3D ***&pppDFBlock,  const V3DLONG *sz_img_sub, const V3DLONG szBlock_x,const V3DLONG szBlock_y,
    const V3DLONG szBlock_z, const int i_interpmethod_img, long long x_offset, long long y_offset, long long z_offset, const long long *sz_img_ori,
    float * &sort_x, float * &sort_y, float * &sort_z)
{
    if (GPU_type)
    {
        int gsz2 = pSubDF->sz2();
        int gsz1 = pSubDF->sz1();
        int gsz0 = pSubDF->sz0();
        clock_t stps_interpolation;
        stps_interpolation = clock();



        gpu_interpolation_sort(gsz2, gsz1, gsz0, pppSubDF, x_bsplinebasis, sz_gridwnd, pppDFBlock, sz_img_sub, szBlock_x, szBlock_y, szBlock_z, i_interpmethod_img, x_offset, y_offset, z_offset, sz_img_ori[2],
            sz_img_ori[1], sz_img_ori[0],sort_x,sort_y,sort_z);


        printf("\t>>interpolation time consume %.2f s\n", (float)(clock() - stps_interpolation) / CLOCKS_PER_SEC);
    }
    else
    {
        for (V3DLONG substart_z = 0; substart_z < pSubDF->sz2() - 1 - 2; substart_z++)
            for (V3DLONG substart_y = 0; substart_y < pSubDF->sz1() - 1 - 2; substart_y++)
                for (V3DLONG substart_x = 0; substart_x < pSubDF->sz0() - 1 - 2; substart_x++)
                {
                    ////bspline interpolate the SubDfBlock to DFBlock
                    //q_dfblcokinterp_bspline(pppSubDF, x_bsplinebasis, sz_gridwnd, substart_x, substart_y, substart_z, pppDFBlock);
                    ////warp image block using DFBlock
                    //q_imgblockwarp_CYF(p_img_sub_4d, sz_img_sub, pppDFBlock, szBlock_x, szBlock_y, szBlock_z, i_interpmethod_img, substart_x, substart_y, substart_z, p_img_warp_4d,
                    //	sz_img_sub_ori, x_offset, y_offset, z_offset, sz_img_ori, sz_img_read_block);
                }
    }
    return true;
}
