// Bigwarp.cpp
// Big warp 's utility functions
// by CYF
// 2022-09-25
#include <QtGui>

//#include "q_warp_affine_tps.h"
#include "q_imgwarp_tps_quicksmallmemory.h"
#include "basic_memory.cpp"//note: should not include .h file, since they are template functions
#include <string>
#include <sstream>
#include "stackutil.h"
#include "Bigwarp.h"
#include "dir_make.h"
#include "RawFmtMngr.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include "newmatap.h"
#include "newmatio.h"

#define _FILE_OFFSET_BITS  64  //20140919

#define b_VERBOSE_PRINT 1
#define ZZBIG 320000 //previous I define it as 1500, so that to limit the size of an image is at most 1.5G //change 2010-05-21 // hang 2011-08-25 6000->10000

#ifdef _MSC_VER       //2010-05-21, by PHC
#include <sys/stat.h>
#include <io.h>
#endif

#define DEFINE_NBYTE2G \
  V3DLONG nBytes2G = (V3DLONG(1024)*V3DLONG(1024)*V3DLONG(1024)-1)*V3DLONG(2);

void MkdirWithPath(const string& FolderPath)
{
#ifdef WIN32

	string PathTrans;
	PathTrans.reserve(FolderPath.size() * 2);

	string::const_iterator itSrcPath = FolderPath.begin();
	while (itSrcPath != FolderPath.end())
	{
		if (*itSrcPath == '/')
		{
			PathTrans.push_back('\\');
		}
		else
		{
			PathTrans.push_back(*itSrcPath);
		}
		++itSrcPath;
	}
	string(PathTrans).swap(PathTrans);
	string fileName = "mkdir " + PathTrans;

#else
    string fileName = "mkdir -p " + FolderPath;

#endif // WIN32

    system(fileName.c_str());
}

bool read_2Dtif_BLOCK(string * &path_tif, unsigned char *&p_img_read_block, unsigned char ****&p_img_read_4d_block, long long * &sz_img_read_block, const long long *sz_img_ori,
	long long * &sz_img_read_size, long long per_sz_block, long long * &sz_img_read_end_block)
{

	if (p_img_read_block)
	{
		printf("WARNNING: output image pointer is not null, original memory it point to will lost!\n");
		p_img_read_block = 0;
	}
	if (p_img_read_4d_block)
	{
		printf("WARNNING: output image pointer is not null, original memory it point to will lost!\n");
		p_img_read_4d_block = 0;
	}

	V3DLONG sz0, sz1, sz2, sz3;

	sz0 = sz_img_read_end_block[0] - sz_img_read_block[0];
	sz1 = sz_img_read_end_block[1] - sz_img_read_block[1];
	sz2 = sz_img_read_end_block[2] - sz_img_read_block[2];
	sz3 = 1;

	sz_img_read_size[0] = sz0;
	sz_img_read_size[1] = sz1;
	sz_img_read_size[2] = sz2;
	sz_img_read_size[3] = sz3;
	//sz_img_read_block[4] = 1;
	//allocate memory
	p_img_read_block = new unsigned char[sz0 * sz1 * sz2 * sz3]();
	if (!p_img_read_block)
	{
		printf("ERROR: Fail to allocate memory for p_img_sub2tar.\n");
		return false;
	}

	//unsigned char ****p_img_sub2tar_4d_block = 0;
	if (!new4dpointer(p_img_read_4d_block, sz0, sz1, sz2, sz3, p_img_read_block))
	{
		printf("ERROR: Fail to allocate memory for the 4d pointer of image.\n");
		if (p_img_read_block) 		{ delete[]p_img_read_block;		p_img_read_block = 0; }
		if (p_img_read_4d_block) 	{ delete4dpointer(p_img_read_4d_block, sz0, sz1, sz2, sz3); }
		return false;
	}

	char *err_rawfmt_c;
	int sV0, sV1, sH0, sH1, sD0, sD1;
	long long stridex, stridexy, stridexyz;
    sH0 = sz_img_read_block[0];
    sV0 = sz_img_read_block[1];
    sD0 = 0;

    sH1 = sz_img_read_block[0] + sz0;
    sV1 = sz_img_read_block[1] + sz1;
    sD1 = 1;
    stridex = sH1 - sH0;
    stridexy = (sV1 - sV0)*(sH1 - sH0);
    stridexyz = (sD1 - sD0)*(sV1 - sV0)*(sH1 - sH0);



#pragma omp parallel for
    for (long long i = 0; i < sz2; i++)
	{

        unsigned char *p_img_sub = 0;
        V3DLONG *sz_img_sub = 0;
        int datatype_sub = 0;
        unsigned char ****p_img_sub_4d = 0;

        long long Zpath;
		Zpath = i + sz_img_read_block[2];


		if ((err_rawfmt_c = copyRawFileBlock2Buffer_cyf((char *)path_tif[Zpath].c_str(), sV0, sV1, sH0, sH1, sD0, sD1, p_img_sub, datatype_sub, sz_img_sub, 0, stridex, stridexy, stridexyz)) != 0)
		{
			printf("ERROR: loadImage() return false in loading [%s].\n", path_tif[Zpath].c_str());
//			return false;
		}

		if (!new4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3], p_img_sub))
		{
			printf("ERROR: Fail to allocate memory for the 4d pointer of image.\n");
			if (p_img_sub_4d) 		{ delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
//			return false;
		}

		for (long long x = 0; x < sz0; x++)
		{
			for (long long y = 0; y < sz1; y++)
			{
				p_img_read_4d_block[0][i][y][x] = p_img_sub_4d[0][0][y][x];
//
			}

		}

		//saveImage(save_path.c_str(), p_img_block, sz_img_block, 1);
		//free memory

		if (p_img_sub_4d) 		{ delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
		if (p_img_sub) 				{ delete[]p_img_sub;			p_img_sub = 0; }
		if (sz_img_sub) 			{ delete[]sz_img_sub;			sz_img_sub = 0; }
	}
	//allocate memory

	//saveImage(save_path.c_str(), p_img_read_block, sz_img_read_block, 1);

	//if (p_img_sub2tar_4d_block) 	{ delete4dpointer(p_img_sub2tar_4d_block, sz0, sz1, sz2, sz3); }
	return true;
}


bool conversion_point(const long long *warp_pos, long long * &read_pos, const Matrix &x4x4_affinematrix, const int type_mode)
{

	//printf("warp�����ĵ����꣺%d,%d,%d\n", centre_x, centre_y, centre_z);
	Matrix x_pt_sub2tar_homo(4, 1), x_pt_sub_homo(4, 1);
	//compute the inverse affine projected coordinate in subject image���������ͶӰ����������ͼ��
	x_pt_sub2tar_homo(1, 1) = warp_pos[0];
	x_pt_sub2tar_homo(2, 1) = warp_pos[1];
	x_pt_sub2tar_homo(3, 1) = warp_pos[2];
	x_pt_sub2tar_homo(4, 1) = 1.0;
	x_pt_sub_homo = x4x4_affinematrix*x_pt_sub2tar_homo;

	double cur_centre_read_pos[3];//x,y,z

	cur_centre_read_pos[0] = x_pt_sub_homo(1, 1);
	cur_centre_read_pos[1] = x_pt_sub_homo(2, 1);
	cur_centre_read_pos[2] = x_pt_sub_homo(3, 1);
	long long centre_read_point[3];
	switch (type_mode)//�м�����Զ����0����ʼ����1
	{
	case 0:
		read_pos[0] = ceil(cur_centre_read_pos[0]);
		read_pos[1] = ceil(cur_centre_read_pos[1]);
		read_pos[2] = ceil(cur_centre_read_pos[2]);
		break;
	case 1:
		read_pos[0] = floor(cur_centre_read_pos[0]);
		read_pos[1] = floor(cur_centre_read_pos[1]);
		read_pos[2] = floor(cur_centre_read_pos[2]);
		break;
	default:
		printf("your type_mode is invalid!");

	}

	if (read_pos[0] < 0)
	{
		read_pos[0] = 0;
	}
	if (read_pos[1] < 0)
	{
		read_pos[1] = 0;
	}
	if (read_pos[2] < 0)
	{
		read_pos[2] = 0;
	}
	return true;
}


bool q_imagewarp_stitch_2DRAW_new(const char* imgSrcFile, string imgSrcFile_save, const long long *sz_img_ori, long long per_block_size, const long long *sz_offset_x, const long long *sz_offset_y,
	const long long Z_Slice_control)
{
	string *path_TIF = 0;
	long long TIF_DEPTH_MAX;
	if (!read_2Dtif_dir(imgSrcFile, path_TIF, TIF_DEPTH_MAX))//��ȡTIF�ļ�Ŀ¼�б����path_TIF������
	{
		printf("ERROR: read_2Dtif_dir() return false in read [%s].\n", imgSrcFile);
		return false;
	}
	char *err_rawfmt_c;
	int X_SIZE, Y_SIZE, Z_SIZE;
	unsigned char *p_img_sub = 0;
	unsigned char *p_img_sub_save_2d = 0;//����ƴ�ɵ�ͼ������1d
	unsigned char *p_img_sub_all = 0;//ÿһ�������ͼ������1d
	V3DLONG *sz_img_sub = 0;
	int datatype_sub = 0;
	unsigned char ****p_img_sub_4d = 0;
	unsigned char ****p_img_sub_4d_save_2d = 0;//����ƴ�ɵ�ͼ������4d
	unsigned char ****p_img_sub_all_4d = 0;//ÿһ�������ͼ������4d
	long long Zpath;
	long long Z_Slice;
	long long Z_Slice_max;
	V3DLONG sz_img_sub_save_2d[4];
	std::stringstream base_2D_TIF_path;
	std::stringstream abs_pos_z;
	Z_SIZE = sz_img_ori[2] % per_block_size == 0 ? sz_img_ori[2] / per_block_size : sz_img_ori[2] / per_block_size + 1;
	Y_SIZE = sz_img_ori[1] % per_block_size == 0 ? sz_img_ori[1] / per_block_size : sz_img_ori[1] / per_block_size + 1;
	X_SIZE = sz_img_ori[0] % per_block_size == 0 ? sz_img_ori[0] / per_block_size : sz_img_ori[0] / per_block_size + 1;

	sz_img_sub_save_2d[0] = sz_img_ori[0];
	sz_img_sub_save_2d[1] = sz_img_ori[1];
	sz_img_sub_save_2d[2] = 1;
	sz_img_sub_save_2d[3] = 1;


	//p_img_affine = new unsigned char[sz_img_output[0] * sz_img_output[1] * sz_img_output[2] * sz_img_output[3]]();

	printf("Z:%d,Y:%d,X:%d\n", Z_SIZE, Y_SIZE, X_SIZE);
	int PJ_mode = X_SIZE*Y_SIZE;
	//int PJ_mode = 16;
	if (PJ_mode <= 15)
	{
		for (int i = 0; i < Z_SIZE; i++)//ÿһ��
		{

			Z_Slice_max = (i + 1)*per_block_size > sz_img_ori[2] ? (sz_img_ori[2] - i*per_block_size) : per_block_size;
			if (p_img_sub_all)
			{
				printf("WARNNING: output image pointer is not null, original memory it point to will lost!\n");
				p_img_sub_all = 0;
			}
			p_img_sub_all = new unsigned char[sz_img_sub_save_2d[0] * sz_img_sub_save_2d[1] * Z_Slice_max * sz_img_sub_save_2d[3]]();

			if (!new4dpointer(p_img_sub_all_4d, sz_img_ori[0], sz_img_ori[1], Z_Slice_max, 1, p_img_sub_all))
			{
				printf("ERROR: Fail to allocate memory for the 4d pointer of image.\n");
				if (p_img_sub_all_4d) 	{ delete4dpointer(p_img_sub_all_4d, sz_img_ori[0], sz_img_ori[1], Z_Slice_max, 1); }
				if (p_img_sub_all) 		{ delete[]p_img_sub_all;		p_img_sub_all = 0; }

				return false;
			}

            for (int k = 0; k < Y_SIZE*X_SIZE; k++)//每一层的切片的块
            {
                unsigned char *p_img_sub = 0;
                unsigned char ****p_img_sub_4d = 0;

                if (p_img_sub)
                {
                    printf("WARNNING: output image pointer is not null, original memory it point to will lost!\n");
                    p_img_sub = 0;
                }

                if (!loadImage((char *)path_TIF[k + i*Y_SIZE*X_SIZE].c_str(), p_img_sub, sz_img_sub, datatype_sub))
                {
                    printf("ERROR: loadImage() return false in loading [%s].\n", path_TIF[k + i*Y_SIZE*X_SIZE].c_str());
                    //return false;
                }

                if (!new4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3], p_img_sub)) //负责搬运的数组
                {
                    printf("ERROR: Fail to allocate memory for the 4d pointer of image.\n");
                    if (p_img_sub_4d)               { delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
                    //return false;
                }

                for (long long z = 0; z < sz_img_sub[2]; z++)
                {
                    for (long long y = 0; y < sz_img_sub[1]; y++)
                    {
                        for (long long x = 0; x < sz_img_sub[0]; x++)
                        {
                            p_img_sub_all_4d[0][z][y + sz_offset_y[k + i*Y_SIZE*X_SIZE]][x + sz_offset_x[k + i*Y_SIZE*X_SIZE]] = p_img_sub_4d[0][z][y][x];
                        }

                    }
                }

                if (p_img_sub_4d)               { delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
                if (p_img_sub)                          { delete[]p_img_sub;                    p_img_sub = 0; }
                if (sz_img_sub)                         { delete[]sz_img_sub;                   sz_img_sub = 0; }

            }


            if (!(saveImage_PJ(Z_Slice_max, i, imgSrcFile_save, per_block_size, p_img_sub_all_4d, sz_img_sub_save_2d)))
            {
                printf("ERROR: Fail to save image of stitching");
                return false;
            }

//            for (int j = 0; j < Z_Slice_max; j++)//每一层的切片
//            {
//                unsigned char *p_img_sub_save_2d = 0;//最终拼成的图像数组1d
//                unsigned char ****p_img_sub_4d_save_2d = 0;//最终拼成的图像数组4d
//
//                if (p_img_sub_save_2d)
//                {
//                    printf("WARNNING: output image pointer is not null, original memory it point to will lost!\n");
//                    p_img_sub_save_2d = 0;
//                }
//
//                p_img_sub_save_2d = new unsigned char[sz_img_sub_save_2d[0] * sz_img_sub_save_2d[1] * sz_img_sub_save_2d[2] * sz_img_sub_save_2d[3]]();
//
//                if (!new4dpointer(p_img_sub_4d_save_2d, sz_img_ori[0], sz_img_ori[1], 1, 1, p_img_sub_save_2d))
//                {
//                    printf("ERROR: Fail to allocate memory for the 4d pointer of image.\n");
//                    if (p_img_sub_4d_save_2d) 	{ delete4dpointer(p_img_sub_4d_save_2d, sz_img_ori[0], sz_img_ori[1], 1, 1); }
//                    if (p_img_sub_save_2d) 		{ delete[]p_img_sub_save_2d;		p_img_sub_save_2d = 0; }
//
//                    //return false;
//                }
//
//                for (long long y = 0; y < sz_img_ori[1]; y++)
//                {
//                    for (long long x = 0; x < sz_img_ori[0]; x++)
//                    {
//                        p_img_sub_4d_save_2d[0][0][y][x] = p_img_sub_all_4d[0][j][y][x];
//                    }
//
//                }
//
//                abs_pos_z.width(7);
//                abs_pos_z.fill('0');
//                abs_pos_z << (int)((i*per_block_size + j) * 10);
//                //base_2D_TIF_path << imgSrcFile_save << "/" << (int)(i*per_block_size + j) << ".raw";
//                base_2D_TIF_path << imgSrcFile_save << "/" << abs_pos_z.str() << ".raw";
//                printf("%d is finish!\n", i*per_block_size + j);
//                saveImage_new(base_2D_TIF_path.str().c_str(), p_img_sub_save_2d, sz_img_sub_save_2d, 1);
//
//                //清除地址
//                string numStr;
//                base_2D_TIF_path.clear();
//                base_2D_TIF_path.str("");
//
//                abs_pos_z.clear();
//                abs_pos_z.str("");
//
//                if (p_img_sub_4d_save_2d) 	{ delete4dpointer(p_img_sub_4d_save_2d, sz_img_ori[0], sz_img_ori[1], 1, 1); }
//                if (p_img_sub_save_2d) 		{ delete[]p_img_sub_save_2d;		p_img_sub_save_2d = 0; }
//
//            }

            if (p_img_sub_all_4d) 	{ delete4dpointer(p_img_sub_all_4d, sz_img_ori[0], sz_img_ori[1], Z_Slice_max, 1); }
			if (p_img_sub_all) 		{ delete[]p_img_sub_all;		p_img_sub_all = 0; }
		}
	}
	else
	{
		for (int i = 0; i < Z_SIZE; i++)//ÿһ��
		{

			Z_Slice_max = (i + 1)*per_block_size > sz_img_ori[2] ? (sz_img_ori[2] - i*per_block_size) : per_block_size;
			printf("Z_Slice_max:%d\n", Z_Slice_max);
			if (Z_Slice_max <= Z_Slice_control)
			{
				if (p_img_sub_all)
				{
					printf("WARNNING: output image pointer is not null, original memory it point to will lost!\n");
					p_img_sub_all = 0;
				}
				p_img_sub_all = new unsigned char[sz_img_sub_save_2d[0] * sz_img_sub_save_2d[1] * Z_Slice_max * sz_img_sub_save_2d[3]]();

				if (!new4dpointer(p_img_sub_all_4d, sz_img_ori[0], sz_img_ori[1], Z_Slice_max, 1, p_img_sub_all))
				{
					printf("ERROR: Fail to allocate memory for the 4d pointer of image.\n");
					if (p_img_sub_all_4d) 	{ delete4dpointer(p_img_sub_all_4d, sz_img_ori[0], sz_img_ori[1], Z_Slice_max, 1); }
					if (p_img_sub_all) 		{ delete[]p_img_sub_all;		p_img_sub_all = 0; }

					return false;
				}


                for (int k = 0; k < Y_SIZE*X_SIZE; k++)//每一层的切片的块
                {
                	if (p_img_sub)
                	{
                		printf("WARNNING: output image pointer is not null, original memory it point to will lost!\n");
                		p_img_sub = 0;
                	}


                	if (!loadImage((char *)path_TIF[k + i*Y_SIZE*X_SIZE].c_str(), p_img_sub, sz_img_sub, datatype_sub))
                	{
                		printf("ERROR: loadImage() return false in loading [%s].\n", path_TIF[k + i*Y_SIZE*X_SIZE].c_str());
                		return false;
                	}
                	//printf("check2");
                	if (!new4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3], p_img_sub))	//负责搬运的数组
                	{
                		printf("ERROR: Fail to allocate memory for the 4d pointer of image.\n");
                		if (p_img_sub_4d) 		{ delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
                		return false;
                	}
                	//printf("check3");
                	for (long long z = 0; z < sz_img_sub[2]; z++)
                	{
                		for (long long y = 0; y < sz_img_sub[1]; y++)
                		{
                			for (long long x = 0; x < sz_img_sub[0]; x++)
                			{
                				p_img_sub_all_4d[0][z][y + sz_offset_y[k + i*Y_SIZE*X_SIZE]][x + sz_offset_x[k + i*Y_SIZE*X_SIZE]] = p_img_sub_4d[0][z][y][x];
                			}

                		}
                	}
                	//printf("check4");

                	if (p_img_sub_4d) 		{ delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
                	if (p_img_sub) 				{ delete[]p_img_sub;			p_img_sub = 0; }
                	if (sz_img_sub) 			{ delete[]sz_img_sub;			sz_img_sub = 0; }

                }

//                if (!(loadImage_PJ(path_TIF, i, Y_SIZE, X_SIZE, p_img_sub_all_4d, sz_offset_x, sz_offset_y)))
//                {
//                    printf("ERROR: Fail to load image of stitching");
//                    return false;
//                }

                if (!(saveImage_PJ(Z_Slice_max, i, imgSrcFile_save, per_block_size, p_img_sub_all_4d, sz_img_sub_save_2d)))
                {
                    printf("ERROR: Fail to save image of stitching");
                    return false;
                }

                if (p_img_sub_all_4d) 	{ delete4dpointer(p_img_sub_all_4d, sz_img_ori[0], sz_img_ori[1], Z_Slice_max, 1); }
				if (p_img_sub_all) 		{ delete[]p_img_sub_all;		p_img_sub_all = 0; }
			}
			else//ͨ��Z_Slice_control����ÿ�ζ�ȡ�Ĳ���
			{
				//int Z_s_max = ceil(Z_Slice_max / Z_Slice_control);
				int Z_s_max = Z_Slice_max % Z_Slice_control == 0 ? Z_Slice_max / Z_Slice_control : Z_Slice_max / Z_Slice_control + 1;
				//Z_SIZE = sz_img_ori[2] % per_block_size == 0 ? sz_img_ori[2] / per_block_size : sz_img_ori[2] / per_block_size + 1;
				printf("Z_s_max:%d\n", Z_s_max);
				for (int Z_s = 0; Z_s < Z_s_max; Z_s++)
				{
					long long Z_Slice_cur = Z_Slice_control*(Z_s + 1) > Z_Slice_max ? (Z_Slice_max - Z_Slice_control*Z_s) : Z_Slice_control;
					printf("Z_Slice_cur:%d\n", Z_Slice_cur);

					if (p_img_sub_all)
					{
						printf("WARNNING: output image pointer is not null, original memory it point to will lost!\n");
						p_img_sub_all = 0;
					}
					p_img_sub_all = new unsigned char[sz_img_sub_save_2d[0] * sz_img_sub_save_2d[1] * Z_Slice_cur * sz_img_sub_save_2d[3]]();

					if (!new4dpointer(p_img_sub_all_4d, sz_img_ori[0], sz_img_ori[1], Z_Slice_cur, 1, p_img_sub_all))
					{
						printf("ERROR: Fail to allocate memory for the 4d pointer of image.\n");
						if (p_img_sub_all_4d) 	{ delete4dpointer(p_img_sub_all_4d, sz_img_ori[0], sz_img_ori[1], Z_Slice_cur, 1); }
						if (p_img_sub_all) 		{ delete[]p_img_sub_all;		p_img_sub_all = 0; }

						return false;
					}

                    if (!(loadImage_PJ_control(Z_Slice_cur, Z_s, Z_Slice_control, path_TIF, i, Y_SIZE, X_SIZE, p_img_sub_all_4d, sz_offset_x, sz_offset_y)))
                    {
                        printf("ERROR: Fail to load image of stitching");
                        return false;
                    }

//                    for (int k = 0; k < Y_SIZE*X_SIZE; k++)//每一层的切片的块
//                    {
//                    	if (p_img_sub)
//                    	{
//                    		printf("WARNNING: output image pointer is not null, original memory it point to will lost!\n");
//                    		p_img_sub = 0;
//                    	}
//
//
//                    	long long Z_sD1 = Z_Slice_cur + Z_s*Z_Slice_control;
//                    	if ((err_rawfmt_c = copyRawFileBlock2Buffer_cyf_PINJIE_new((char *)path_TIF[k + i*Y_SIZE*X_SIZE].c_str(), 0, 0, Z_s*Z_Slice_control, Z_sD1, p_img_sub, datatype_sub, sz_img_sub, 0)) != 0)
//                    	{
//                    	   printf("ERROR: loadImage() return false in loading [%s].\n", path_TIF[k + i*Y_SIZE*X_SIZE].c_str());
//                    	   return false;
//                    	}
//                    	//printf("check2");
//                    	if (!new4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3], p_img_sub))	//负责搬运的数组
//                    	{
//                    		printf("ERROR: Fail to allocate memory for the 4d pointer of image.\n");
//                    		if (p_img_sub_4d) 		{ delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
//                    		return false;
//                    	}
//                    	//printf("check3");
//                    	for (long long z = 0; z < sz_img_sub[2]; z++)
//                    	{
//                    		for (long long y = 0; y < sz_img_sub[1]; y++)
//                    		{
//                    			for (long long x = 0; x < sz_img_sub[0]; x++)
//                    			{
//                    				p_img_sub_all_4d[0][z][y + sz_offset_y[k + i*Y_SIZE*X_SIZE]][x + sz_offset_x[k + i*Y_SIZE*X_SIZE]] = p_img_sub_4d[0][z][y][x];
//                    			}
//
//                    		}
//                    	}
//                    	//printf("check4");
//
//                    	if (p_img_sub_4d) 		{ delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
//                    	if (p_img_sub) 				{ delete[]p_img_sub;			p_img_sub = 0; }
//                    	if (sz_img_sub) 			{ delete[]sz_img_sub;			sz_img_sub = 0; }
//
//                    }


                    if (!(saveImage_PJ_control(Z_Slice_cur, Z_s, Z_Slice_control, i, imgSrcFile_save, per_block_size, p_img_sub_all_4d, sz_img_sub_save_2d)))
                    {
                        printf("ERROR: Fail to save image of stitching");
                        return false;
                    }

					if (p_img_sub_all_4d) 	{ delete4dpointer(p_img_sub_all_4d, sz_img_ori[0], sz_img_ori[1], Z_Slice_cur, 1); }
					if (p_img_sub_all) 		{ delete[]p_img_sub_all;		p_img_sub_all = 0; }
				}
				
			}
			
		}
	}


	if (path_TIF){ delete[]path_TIF; }
	return true;
}

void bubbleSort(vector<long long> &q)
{
	for (int i = q.size() - 1; i > 0; i--){
		bool flag = false;
		for (int j = 0; j + 1 <= i; j++){
			if (q[j] > q[j + 1]){
				swap(q[j], q[j + 1]);
				flag = true;
			}
		}
		if (!flag)
			break;
	}
}

bool saveImage_new(const char filename[], const unsigned char * data1d, const V3DLONG * sz, const int datatype)
{
	if (!data1d || !filename || !sz)
	{
		printf("This image data is empty or the file name or the size pointer is invalid. Nothing done.\n");
		return false;
	}

	int dt;
	ImagePixelType curtype;
	switch (datatype)
	{
	case 1:  dt = 1; curtype = V3D_UINT8; break;
	case 2:  dt = 2; curtype = V3D_UINT16; break;
	case 4:  dt = 4; curtype = V3D_FLOAT32; break;
	default:
		printf("The data type is unsupported. Nothing done.\n");
		return false;
		break;
	}

	const char * curFileSuffix = getSuffix((char *)filename);
	//	if (b_VERBOSE_PRINT)
	//		printf("The current input file has the surfix [%s]\n", curFileSuffix);

	//it seems when curFileSuffix is NULL then strcasecmp() has a crashing bug. thus check now. 20120410. by PHC
	if (curFileSuffix && (strcasecmp(curFileSuffix, "tif") == 0 || strcasecmp(curFileSuffix, "tiff") == 0)) //write tiff stacks
	{
		if (saveStack2Tif(filename, data1d, sz, dt))
		{
			printf("Error happens in writing TIF file [%s]. Stop. \n", filename);
			return false;
		}
	}
	else if (curFileSuffix && strcasecmp(curFileSuffix, "raw5") == 0) //write .raw5 data
	{
		if (saveStack2Raw5d(filename, data1d, sz, dt))
		{
			printf("Error happens in writing V3D .raw5 file [%s]. Stop. \n", filename);
			return false;
		}
	}
#ifdef _ALLOW_WORKMODE_MENU_
	else if (curFileSuffix && (strcasecmp(curFileSuffix, "v3dpbd") == 0)) //  v3dpbd - pack-bit-difference encoding for sparse stacks
		// || strcasecmp(curFileSuffix, "mp4")==0) ) //to add mp4 later
	{
		v3d_msg("prepare for pbd file saving", 0);
		ImageLoaderBasic imageLoader;
		if (imageLoader.saveStack2RawPBD(filename, curtype, (unsigned char *)data1d, sz)) {
			printf("Error happens in v3dpbd file saving. Stop. \n");
			return false;
		}
	}
#endif
	else //then assume it is Hanchuan's RAW format
	{
		//if (b_VERBOSE_PRINT)
		//printf("The data is not with a known Vaa3D format, -- now this program assumes it is a Vaa3D RAW format. \n");
		if (saveStack2Raw_cyf(filename, data1d, sz, dt) != 0) //0 is no error //note that as I updated the saveStack2Raw to RAW-4-byte, the actual mask file cannot be read by the old wano program, i.e. the wano must be updated on Windows machine as well. 060921
		{
			printf("Error happens in writing RAW file stack [defined by Hanchuan Peng] [%s].\n", filename);
			return false;
		}
	}

	return true;
}

/* The following is the core function for image stack writing */

int saveStack2Raw_cyf(const char * filename, const unsigned char * img, const V3DLONG * sz, int datatype)
{


	int berror = 0;
	V3DLONG i;

	FILE * fid = fopen(filename, "wb");
	if (!fid)
	{
		printf("Fail to open file for writing.\n");
		berror = 1;
		return berror;
	}



	char formatkey[] = "raw_image_stack_by_hpeng";
	int lenkey = strlen(formatkey);

	V3DLONG nwrite = fwrite(formatkey, 1, lenkey, fid);
	if (nwrite != lenkey)
	{
		printf("File write error.\n");
		berror = 1;
		return berror;
	}

	char endianCodeMachine = checkMachineEndian();
	if (endianCodeMachine != 'B' && endianCodeMachine != 'L')
	{
		printf("This program only supports big- or little- endian but not other format. Cannot save data on this machine.\n");
		berror = 1;
		return berror;
	}

	nwrite = fwrite(&endianCodeMachine, 1, 1, fid);
	if (nwrite != 1)
	{
		printf("Error happened in file writing.\n");
		berror = 1;
		return berror;
	}

	short int dcode = (short int)datatype;
	if (dcode != 1 && dcode != 2 && dcode != 4)
	{
		printf("Unrecognized data type code [%d]. This code is not supported in this version.\n", dcode);
		berror = 1;
		return berror;
	}

	nwrite = fwrite(&dcode, 2, 1, fid);
	if (nwrite != 1)
	{
		printf("Writing file error.\n");
		berror = 1;
		return berror;
	}

	V3DLONG unitSize = datatype;

	//short int mysz[4];
	BIT32_UNIT mysz[4];
	for (i = 0; i<4; i++) mysz[i] = (BIT32_UNIT)sz[i];
	nwrite = fwrite(mysz, 4, 4, fid);
	if (nwrite != 4)
	{
		printf("Writing file error.\n");
		berror = 1;
		return berror;
	}

	V3DLONG totalUnit = 1;
	for (i = 0; i<4; i++)
	{
		totalUnit *= sz[i];
	}

	nwrite = fwrite(img, unitSize, totalUnit, fid);
	if (nwrite != totalUnit)
	{
		printf("Something wrong in file writing. The program wrote [%ld data points] but the file says there should be [%ld data points].\n", nwrite, totalUnit);
		berror = 1;
		return berror;
	}


	fclose(fid);

	return berror;
}

bool loadImage_cyf(char imgSrcFile[], unsigned char *& data1d, V3DLONG * &sz, int & datatype, int start_x, int start_y, int end_x, int end_y)
{
	if (data1d)
	{
		printf("Warning: The pointer for 1d data storage is not empty. This pointer will be freed first and the  reallocated. \n");
		delete[]data1d; data1d = 0;
	}
	if (sz)
	{
		printf("Warning: The pointer for size is not empty. This pointer will be freed first and the  reallocated. \n");
		delete[]sz; sz = 0;
	}

	unsigned char *tmp_data1d = 0;
	V3DLONG * tmp_sz = 0;
	V3DLONG * tmp_sz_ori = 0;
	int tmp_datatype = 0;
	bool b_5d = false;

	const char * curFileSuffix = getSuffix(imgSrcFile);



	if (curFileSuffix && (strcasecmp(curFileSuffix, "tif") == 0 || strcasecmp(curFileSuffix, "tiff") == 0)) //read tiff stacks
	{
		if (!ensure_file_exists_and_size_not_too_big(imgSrcFile, (V3DLONG)1024 * 1024 * ZZBIG)) //tif file at most should be 900M bytes
		{
			printf("The tif file may not exist or may be too big to load [sz threshold=%ld bytes].\n", (V3DLONG)1024 * 1024 * ZZBIG);
			return false;
		}

	}
#ifdef _ALLOW_WORKMODE_MENU_    
	else if (curFileSuffix && ImageLoaderBasic::hasPbdExtension(imgSrcFile)) // read v3dpbd - pack-bit-difference encoding for sparse stacks
	{
		v3d_msg("prepare for pbd file loading", 0);
		Image4DSimple *tmpimg = new Image4DSimple;

		ImageLoaderBasic imageLoader;
		if (!imageLoader.loadImage(tmpimg, imgSrcFile)) {
			printf("Error happens in v3dpbd file reading. Stop. \n");
			return false;
		}
		// The following few lines are to avoid disturbing the existing code below

		tmp_data1d = tmpimg->getRawData();
		tmp_datatype = tmpimg->getDatatype();
		tmp_sz = new V3DLONG[4];
		tmp_sz[0] = tmpimg->getXDim();
		tmp_sz[1] = tmpimg->getYDim();
		tmp_sz[2] = tmpimg->getZDim();
		tmp_sz[3] = tmpimg->getCDim();
	}
#endif
	else if (curFileSuffix && strcasecmp(curFileSuffix, "lsm") == 0) //read lsm stacks
	{
		if (!ensure_file_exists_and_size_not_too_big(imgSrcFile, (V3DLONG)1024 * 1024 * ZZBIG)) //2047 //lsm file at most should be 900M bytes
		{
			printf("The lsm file may not exist or may be too big to load [sz threshold=%ld bytes].\n", (V3DLONG)1024 * 1024 * ZZBIG);
			return false;
		}
		if (loadLsm2Stack(imgSrcFile, tmp_data1d, tmp_sz, tmp_datatype))
		{
			printf("Error happens in LSM file reading. Stop. \n");
			return false;
		}
	}
	else if (curFileSuffix && strcasecmp(curFileSuffix, "mrc") == 0) //read MRC stacks
	{
		if (!ensure_file_exists_and_size_not_too_big(imgSrcFile, (V3DLONG)1024 * 1024 * ZZBIG)) //MRC file at most should be 1.5G bytes
		{
			printf("The MRC file may not exist or may be too big to load.\n");
			return false;
		}

		if (loadMRC2Stack(imgSrcFile, tmp_data1d, tmp_sz, tmp_datatype))
		{
			if (b_VERBOSE_PRINT)
				printf("The data doesn't look like a correct MRC file. \n");
			return false;
		}
	}
	else if (curFileSuffix && strcasecmp(curFileSuffix, "raw5") == 0) //read lsm stacks
	{
		if (!ensure_file_exists_and_size_not_too_big(imgSrcFile, (V3DLONG)1024 * 1024 * ZZBIG)) //
		{
			printf("The lsm file may not exist or may be too big to load [sz threshold=%ld bytes].\n", (V3DLONG)1024 * 1024 * ZZBIG);
			return false;
		}
		if (loadRaw5d2Stack(imgSrcFile, tmp_data1d, tmp_sz, tmp_datatype))
		{
			printf("Error happens in V3D .raw5 (5D) file reading. Stop. \n");
			return false;
		}
		b_5d = true;
	}
	else //then assume it is Hanchuan's RAW format
	{
		if (b_VERBOSE_PRINT)
			printf("The data is not with a TIF/LSM surfix, -- now this program assumes it is RAW format defined by Hanchuan Peng. \n");
		if (!ensure_file_exists_and_size_not_too_big(imgSrcFile, (V3DLONG)1024 * 1024 * ZZBIG)) //RAW file at most should be 1.5G bytes
		{
			printf("The RAW file may not exist or may be too big to load.\n");
			return false;
		}

		if (loadRaw2Stack(imgSrcFile, tmp_data1d, tmp_sz, tmp_datatype))
		{
			if (b_VERBOSE_PRINT)
				printf("The data doesn't look like a correct 4-byte-size RAW file. Try 2-byte-raw. \n");
			if (loadRaw2Stack_2byte(imgSrcFile, tmp_data1d, tmp_sz, tmp_datatype))
			{
				printf("Error happens in reading 2-byte-size RAW file. Stop. \n");
				return false;
			}
		}
	}


	//copy output data

	switch (tmp_datatype)
	{
	case 1:
		datatype = 1;
		break;

	case 2:
		datatype = 2;
		break;

	case 4:
		datatype = 4;
		break;

	default:
		printf("Something wrong with the program, -- should NOT display this message at all. Check your program. \n");
		if (data1d) { delete[]data1d; data1d = 0; }
		if (tmp_sz) { delete[]tmp_sz; tmp_sz = 0; }
		if (tmp_sz_ori) { delete[]tmp_sz_ori; tmp_sz_ori = 0; }
		if (sz) { delete[]sz; sz = 0; }
		return false;
	}

	sz = new V3DLONG[5];
	sz[0] = tmp_sz[0];
	sz[1] = tmp_sz[1];
	sz[2] = tmp_sz[2];
	sz[3] = tmp_sz[3]; //no longer merge the 3rd and 4th dimensions
	sz[4] = (b_5d) ? tmp_sz[4] : 1; //090802

	data1d = tmp_data1d;


	if (tmp_sz) { delete[]tmp_sz; tmp_sz = 0; }
	if (tmp_sz_ori) { delete[]tmp_sz_ori; tmp_sz_ori = 0; }

	return true;
}

bool loadImage_PJ(string *&path_TIF, int block_index, int Y_SIZE, int X_SIZE, unsigned char ****p_img_sub_all_4d, const long long *sz_offset_x, const long long *sz_offset_y)
{
#pragma omp parallel for
    for (int k = 0; k < Y_SIZE*X_SIZE; k++)//每一层的切片的块
    {
        unsigned char *p_img_sub = 0;
        unsigned char ****p_img_sub_4d = 0;
        V3DLONG *sz_img_sub = 0;
        int datatype_sub = 0;

        if (!loadImage((char *)path_TIF[k + block_index*Y_SIZE*X_SIZE].c_str(), p_img_sub, sz_img_sub, datatype_sub))
        {
            printf("ERROR: loadImage() return false in loading [%s].\n", path_TIF[k + block_index*Y_SIZE*X_SIZE].c_str());
            //return false;
        }

        if (!new4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3], p_img_sub)) //负责搬运的数组
        {
            printf("ERROR: Fail to allocate memory for the 4d pointer of image.\n");
            if (p_img_sub_4d)               { delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
            //return false;
        }

        for (long long z = 0; z < sz_img_sub[2]; z++)
        {
            for (long long y = 0; y < sz_img_sub[1]; y++)
            {
                for (long long x = 0; x < sz_img_sub[0]; x++)
                {
                    p_img_sub_all_4d[0][z][y + sz_offset_y[k + block_index*Y_SIZE*X_SIZE]][x + sz_offset_x[k + block_index*Y_SIZE*X_SIZE]] = p_img_sub_4d[0][z][y][x];
                }

            }
        }

        if (p_img_sub_4d)               { delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
        if (p_img_sub)                          { delete[]p_img_sub;                    p_img_sub = 0; }
        if (sz_img_sub)                         { delete[]sz_img_sub;                   sz_img_sub = 0; }

    }

    return true;

}

bool saveImage_PJ(long long Z_Slice_max, int block_index, string imgSrcFile_save, long long per_block_size, unsigned char ****p_img_sub_all_4d, V3DLONG *sz_img_sub_save_2d)
{
#pragma omp parallel for
    for (int j = 0; j < Z_Slice_max; j++)//ÿһ�����Ƭ
    {
        unsigned char *p_img_sub_save_2d = 0;
        unsigned char ****p_img_sub_4d_save_2d = 0;

        p_img_sub_save_2d = new unsigned char[sz_img_sub_save_2d[0] * sz_img_sub_save_2d[1] * sz_img_sub_save_2d[2] * sz_img_sub_save_2d[3]]();

        if (!new4dpointer(p_img_sub_4d_save_2d, sz_img_sub_save_2d[0], sz_img_sub_save_2d[1], sz_img_sub_save_2d[2], sz_img_sub_save_2d[3], p_img_sub_save_2d))
        {
            printf("ERROR: Fail to allocate memory for the 4d pointer of image.\n");
            if (p_img_sub_4d_save_2d)       { delete4dpointer(p_img_sub_4d_save_2d, sz_img_sub_save_2d[0], sz_img_sub_save_2d[1], sz_img_sub_save_2d[2], sz_img_sub_save_2d[3]); }
            if (p_img_sub_save_2d)          { delete[]p_img_sub_save_2d;            p_img_sub_save_2d = 0; }

        }

        for (long long y = 0; y < sz_img_sub_save_2d[1]; y++)
        {
            for (long long x = 0; x < sz_img_sub_save_2d[0]; x++)
            {
                p_img_sub_4d_save_2d[0][0][y][x] = p_img_sub_all_4d[0][j][y][x];
            }
        }

        std::stringstream base_2D_TIF_path;
        std::stringstream abs_pos_z;

        abs_pos_z.width(7);
        abs_pos_z.fill('0');
        abs_pos_z << (int)((block_index*per_block_size + j) * 10);
        //base_2D_TIF_path << imgSrcFile_save << "/" << (int)(i*per_block_size + j) << ".raw";
        base_2D_TIF_path << imgSrcFile_save << "/" << abs_pos_z.str() << ".raw";
//        printf("%d is finish!\n", block_index*per_block_size + j);
        saveImage_new(base_2D_TIF_path.str().c_str(), p_img_sub_save_2d, sz_img_sub_save_2d, 1);

        //�����ַ
        base_2D_TIF_path.clear();
        base_2D_TIF_path.str("");

        abs_pos_z.clear();
        abs_pos_z.str("");

        if (p_img_sub_4d_save_2d)       { delete4dpointer(p_img_sub_4d_save_2d, sz_img_sub_save_2d[0], sz_img_sub_save_2d[1], sz_img_sub_save_2d[2], sz_img_sub_save_2d[3]); }
        if (p_img_sub_save_2d)          { delete[]p_img_sub_save_2d;            p_img_sub_save_2d = 0; }

    }
    return true;

}

bool loadImage_PJ_control(long long Z_Slice_cur, int Z_s, long long Z_Slice_control, string *&path_TIF, int block_index, int Y_SIZE, int X_SIZE, unsigned char ****p_img_sub_all_4d, const long long *sz_offset_x, const long long *sz_offset_y)
{
#pragma omp parallel for
    for (int k = 0; k < Y_SIZE*X_SIZE; k++)//每一层的切片的块
    {
        unsigned char *p_img_sub = 0;
        unsigned char ****p_img_sub_4d = 0;
        V3DLONG *sz_img_sub = 0;
        int datatype_sub = 0;
        char *err_rawfmt_c;

        long long Z_sD1 = Z_Slice_cur + Z_s*Z_Slice_control;
        if ((err_rawfmt_c = copyRawFileBlock2Buffer_cyf_PINJIE_new((char *)path_TIF[k + block_index*Y_SIZE*X_SIZE].c_str(), 0, 0, Z_s*Z_Slice_control, Z_sD1, p_img_sub, datatype_sub, sz_img_sub, 0)) != 0)
        {
            printf("ERROR: loadImage() return false in loading [%s].\n", path_TIF[k + block_index*Y_SIZE*X_SIZE].c_str());
            //return false;
        }

        if (!new4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3], p_img_sub)) //负责搬运的数组
        {
            printf("ERROR: Fail to allocate memory for the 4d pointer of image.\n");
            if (p_img_sub_4d)               { delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
            //return false;
        }

        for (long long z = 0; z < sz_img_sub[2]; z++)
        {
            for (long long y = 0; y < sz_img_sub[1]; y++)
            {
                for (long long x = 0; x < sz_img_sub[0]; x++)
                {
                    p_img_sub_all_4d[0][z][y + sz_offset_y[k + block_index*Y_SIZE*X_SIZE]][x + sz_offset_x[k + block_index*Y_SIZE*X_SIZE]] = p_img_sub_4d[0][z][y][x];
                }

            }
        }

        if (p_img_sub_4d)               { delete4dpointer(p_img_sub_4d, sz_img_sub[0], sz_img_sub[1], sz_img_sub[2], sz_img_sub[3]); }
        if (p_img_sub)                          { delete[]p_img_sub;                    p_img_sub = 0; }
        if (sz_img_sub)                         { delete[]sz_img_sub;                   sz_img_sub = 0; }

    }

    return true;

}

bool saveImage_PJ_control(long long Z_Slice_cur, int Z_s, long long Z_Slice_control, int block_index, string imgSrcFile_save, long long per_block_size, unsigned char ****p_img_sub_all_4d, V3DLONG *sz_img_sub_save_2d)
{
#pragma omp parallel for
    for (int j = 0; j < Z_Slice_cur; j++)//每一层的切片
    {
        unsigned char *p_img_sub_save_2d = 0;
        unsigned char ****p_img_sub_4d_save_2d = 0;



        p_img_sub_save_2d = new unsigned char[sz_img_sub_save_2d[0] * sz_img_sub_save_2d[1] * sz_img_sub_save_2d[2] *
                                              sz_img_sub_save_2d[3]]();

        if (!new4dpointer(p_img_sub_4d_save_2d, sz_img_sub_save_2d[0], sz_img_sub_save_2d[1], sz_img_sub_save_2d[2],
                          sz_img_sub_save_2d[3], p_img_sub_save_2d)) {
            printf("ERROR: Fail to allocate memory for the 4d pointer of image.\n");
            if (p_img_sub_4d_save_2d) {
                delete4dpointer(p_img_sub_4d_save_2d, sz_img_sub_save_2d[0], sz_img_sub_save_2d[1],
                                sz_img_sub_save_2d[2], sz_img_sub_save_2d[3]);
            }
            if (p_img_sub_save_2d) {
                delete[]p_img_sub_save_2d;
                p_img_sub_save_2d = 0;
            }

            //return false;
        }

        for (long long y = 0; y < sz_img_sub_save_2d[1]; y++) {
            for (long long x = 0; x < sz_img_sub_save_2d[0]; x++) {
                p_img_sub_4d_save_2d[0][0][y][x] = p_img_sub_all_4d[0][j][y][x];
            }

        }

        std::stringstream base_2D_TIF_path;
        std::stringstream abs_pos_z;

        abs_pos_z.width(7);
        abs_pos_z.fill('0');
        abs_pos_z << (int) ((block_index * per_block_size + Z_s * Z_Slice_control + j) * 10);
        //base_2D_TIF_path << imgSrcFile_save << "/" << (int)(i*per_block_size + j) << ".raw";
        base_2D_TIF_path << imgSrcFile_save << "/" << abs_pos_z.str() << ".raw";
//        printf("%d is finish!\n", block_index * per_block_size + Z_s * Z_Slice_control + j);
        saveImage_new(base_2D_TIF_path.str().c_str(), p_img_sub_save_2d, sz_img_sub_save_2d, 1);

        //清除地址
        string numStr;
        base_2D_TIF_path.clear();
        base_2D_TIF_path.str("");

        abs_pos_z.clear();
        abs_pos_z.str("");

        if (p_img_sub_4d_save_2d) {
            delete4dpointer(p_img_sub_4d_save_2d, sz_img_sub_save_2d[0], sz_img_sub_save_2d[1], sz_img_sub_save_2d[2],
                            sz_img_sub_save_2d[3]);
        }
        if (p_img_sub_save_2d) {
            delete[]p_img_sub_save_2d;
            p_img_sub_save_2d = 0;
        }

    }
    return true;
}