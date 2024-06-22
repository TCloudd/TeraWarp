#ifdef _WIN32
#include "dirent_win.h"
#else
#include <dirent.h>
#endif
//#include "q_warp_affine_tps.h"
#include <iostream>
#include <string>
#include <list>
#include <sstream>
using namespace std;
bool read_RES_dir(char* &imgSrcFile, string *** &path_Z)
{
	string tmp;
	int DEPTH_Y, DEPTH_X, DEPTH_Z;
	DIR *cur_dir_lev3_Y = opendir(imgSrcFile);
	dirent *entry_lev3;
	dirent *entry_lev3_z;
	list<string> entries_lev3;
	string entry_Y, entry_X, entry_Z;
	string *path_Y;
	string **path_X;
	//string ***path_Z;
	if (path_Z != 0) { return false; } //if the "path_Z" is not empty initially, then do nothing and return un-successful

	//*************************************************��ȡY�ļ�
	while ((entry_lev3 = readdir(cur_dir_lev3_Y)) != 0)//��ȡYĿ¼
	{
		tmp = entry_lev3->d_name;
		if (tmp.compare(".") != 0 && tmp.compare("..") != 0 && tmp != "mdata.bin"&& tmp != ".iim.format")//&&tmp.find(".") != string::npos)
		{
			entries_lev3.push_back(tmp);
		}
	
	}

	entries_lev3.sort();//sort() ��list����
	DEPTH_Y = (int)entries_lev3.size();
	//printf("%d\n", DEPTH_Y);
	closedir(cur_dir_lev3_Y);

	path_Y = new string[DEPTH_Y];
	for (int z = 0; z<DEPTH_Y; z++)
	{
		path_Y[z] = imgSrcFile;
		path_Y[z] += "/";
		entry_Y = entries_lev3.front();//���ص�һ��Ԫ��
		path_Y[z] += entry_Y;
		entries_lev3.pop_front();//ɾ����һ��Ԫ��
		cout << path_Y[z].c_str() << endl;
	}
	entries_lev3.clear();//ɾ������Ԫ��
	//*************************************************

	//*************************************************��ȡX�ļ�
	path_X = new string*[DEPTH_Y];
	path_Z = new string**[DEPTH_Y];

	for (int z = 0; z < DEPTH_Y; z++)
	{
		DIR *cur_dir_lev3_X = opendir(path_Y[z].c_str());
		while ((entry_lev3 = readdir(cur_dir_lev3_X)) != 0)//��ȡXĿ¼
		{
			tmp = entry_lev3->d_name;
			if (tmp.compare(".") != 0 && tmp.compare("..") != 0 && tmp != "mdata.bin"&& tmp != ".iim.format")//&&tmp.find(".") != string::npos)
			{
				entries_lev3.push_back(tmp);
				//cout << tmp.c_str() << endl;
			}

		}
		entries_lev3.sort();//sort() ��list����
		DEPTH_X = (int)entries_lev3.size();
		//printf("%d\n", DEPTH_X);
		closedir(cur_dir_lev3_X);

		path_X[z] = new string[DEPTH_X];
		for (int s = 0; s<DEPTH_X; s++)
		{
			path_X[z][s] = path_Y[z];
			path_X[z][s] += "/";
			entry_X = entries_lev3.front();//���ص�һ��Ԫ��
			path_X[z][s] += entry_X;
			entries_lev3.pop_front();//ɾ����һ��Ԫ��
		}
		entries_lev3.clear();//ɾ������Ԫ��
		//*************************************************

		//*************************************************��ȡZ�ļ�
		path_Z[z] = new string*[DEPTH_X];

		for (int k = 0; k < DEPTH_X; k++)
		{
			DIR *cur_dir_lev3_Z = opendir(path_X[z][k].c_str());
			while ((entry_lev3_z = readdir(cur_dir_lev3_Z)) != 0)//��ȡZĿ¼
			{
				tmp = entry_lev3_z->d_name;
				if (tmp.compare(".") != 0 && tmp.compare("..") != 0 && tmp != "mdata.bin"&& tmp != ".iim.format")//&&tmp.find(".") != string::npos)
				{
					entries_lev3.push_back(tmp);
					//cout << tmp.c_str() << endl;
				}
			}
			entries_lev3.sort();//sort() ��list����
			DEPTH_Z = (int)entries_lev3.size();
			//printf("%d\n", DEPTH_Z);
			closedir(cur_dir_lev3_Z);

			path_Z[z][k] = new string[DEPTH_Z];
		
			for (int p = 0; p<DEPTH_Z; p++)
			{
				path_Z[z][k][p] = path_X[z][k];
				path_Z[z][k][p] += "/";
				entry_Z = entries_lev3.front();//���ص�һ��Ԫ��
				path_Z[z][k][p] += entry_Z;
				entries_lev3.pop_front();//ɾ����һ��Ԫ��
			}
			entries_lev3.clear();//ɾ������Ԫ��
		}
	}
	//*************************************************
	
	

	cout << path_Z[2][1][0].c_str() << endl;
	//cout << name[2].c_str() << endl;
	return true;
	//printf("%s\n", test);
}

void eachnumber(long long *sz)
{
	int x_num[7];
	int y_num[5];
	int z_num[6];
	int x, y, z;
	x = y = z = 0;
	int block_size = 100;
	for (int i = 0; i < 7; i++)
	{
		x_num[i] = x;
//		printf("%d\n", x_num[i]);
		x += block_size;
		if (x > sz[0])
		{
			x = sz[0];
		}
		
	}
	for (int i = 0; i < 5; i++)
	{
		y_num[i] = y;
//		printf("%d\n", y_num[i]);
		y += block_size;
		if (y > sz[1])
		{
			y = sz[1];
		}
		
	}
	for (int i = 0; i < 6; i++)
	{
		z_num[i] = z;
//		printf("%d\n", z_num[i]);
		z += block_size;
		if (z > sz[2])
		{
			z = sz[2];
		}
		
	}
}

bool read_2Dtif_dir(const char* imgSrcFile, string * &path_tif, long long &DEPTH_TIF)
{
	DIR *cur_dir_tif = opendir(imgSrcFile);
	dirent *entry_tif;
	string tmp;
	list<string> entries_tif;
	//int DEPTH_TIF;
	string entry_TIF;
	//dirent *entry_lev3;
	if (path_tif != 0) { return false; } //if the "path_Z" is not empty initially, then do nothing and return un-successful

	//*************************************************��ȡTIF�ļ���
	while ((entry_tif = readdir(cur_dir_tif)) != 0)//��ȡTIF�ļ���Ŀ¼
	{
		tmp = entry_tif->d_name;
		if (tmp.compare(".") != 0 && tmp.compare("..") != 0 && tmp != "mdata.bin"&& tmp != ".iim.format")//&&tmp.find(".") != string::npos)
		{
			entries_tif.push_back(tmp);
		}

	}

	entries_tif.sort();//sort() ��list����
	DEPTH_TIF = (long long)entries_tif.size();
	//printf("%d\n", DEPTH_Y);
	closedir(cur_dir_tif);

	path_tif = new string[DEPTH_TIF];
	for (int z = 0; z<DEPTH_TIF; z++)
	{
		path_tif[z] = imgSrcFile;
		path_tif[z] += "/";
		entry_TIF = entries_tif.front();//���ص�һ��Ԫ��
		path_tif[z] += entry_TIF;
		entries_tif.pop_front();//ɾ����һ��Ԫ��
		//cout << path_tif[z].c_str() << endl;
	}
	entries_tif.clear();//ɾ������Ԫ��
	//*************************************************
	//cout << path_tif[0].c_str() << endl;
	//cout << name[2].c_str() << endl;
	return true;
}