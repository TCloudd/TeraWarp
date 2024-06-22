/* binarization_plugin.cpp
 * This is a test plugin, you can use it as a demo.
 * 2018-10-12 : by Guochanghao
 */
 
#include "v3d_message.h"
#include <iostream>
#include <vector>
#include "newmat.h"
#include "stackutil.h"
#include "q_warp_affine_tps.h"
#include "dir_make.h"
#include "RawFmtMngr.h"
#include "q_paradialog_bigimagewarp.h"
#include "q_imgwarp_tps_quicksmallmemory.h"
#include "plugin_Bigimagewarp.h"

//extern "C" int gpu_A_i_new(int ncpt, const Matrix &A, Matrix &A_i);
using namespace std;
Q_EXPORT_PLUGIN2(Bigimagewarp, BigImageWarPPlugin);

bool bigimagewarp(V3DPluginCallback2 & callback,  QString qs_filename_img_sub,QString qs_filename_tpscpt_sub,QString qs_filename_tpscpt_tar,
                  int i_interpmethod_df,int i_interpmethod_img, int GPU_TYPE, int i_program_mode, int i_downsample_mode, QString qs_filename_img_output);
void printHelp();
 
QStringList BigImageWarPPlugin::menulist() const
{
	return QStringList() 
		<<tr("bigimagewarp")
		<<tr("about");
}

QStringList BigImageWarPPlugin::funclist() const
{
	return QStringList()
		<<tr("bigimagewarp")
		<<tr("help");
}

void BigImageWarPPlugin::domenu(const QString &menu_name, V3DPluginCallback2 &callback, QWidget *parent)
{
    if (menu_name == tr("bigimagewarp"))
    {
        CParaDialog_bigimagewarp DLG_bigimagewarp(callback,parent);
        if(DLG_bigimagewarp.exec()!=QDialog::Accepted)	return;
        QString qs_filename_img_sub = DLG_bigimagewarp.lineEdit_img_sub->text();
        QString qs_filename_img_warp = DLG_bigimagewarp.lineEdit_img_warp->text();
        QString qs_filename_marker_sub = DLG_bigimagewarp.lineEdit_marker_sub->text();
        QString qs_filename_marker_tar = DLG_bigimagewarp.lineEdit_marker_tar->text();
//        int interpmethod_df = DLG_bigimagewarp.radioButton_df_interp_bspline->isChecked();
//        int interpmethod_img = DLG_bigimagewarp.radioButton_img_interp_nn->isChecked();
        int interpmethod_df = 1;
        int interpmethod_img = 0;
        int device_mode = DLG_bigimagewarp.radioButton_device_gpu->isChecked();
        int program_mode = DLG_bigimagewarp.radioButton_affine->isChecked();
//        int downsample_mode = DLG_bigimagewarp.radioButton_affine_2->isChecked();
        int downsample_mode = 0;

        if(!bigimagewarp(callback,  qs_filename_img_sub,qs_filename_marker_sub,qs_filename_marker_tar,
                interpmethod_df, interpmethod_img, device_mode, program_mode, downsample_mode, qs_filename_img_warp))
        {
            v3d_msg(tr("ERROR: bigimagewarp() return false!"));
            return;
        }
    }
	else
	{
        v3d_msg(tr("Warp big image based on given markers. "
            "Developed by , 2023-6-17"));
	}
}

bool BigImageWarPPlugin::dofunc(const QString & func_name, const V3DPluginArgList & input, V3DPluginArgList & output, V3DPluginCallback2 & callback,  QWidget * parent)
{
    if (func_name == tr("bigimagewarp"))
    {
        cout<<"============== Welcome to bigimagewarp function ================="<<endl;
        if(input.size()!=2 || output.size()!=1 || ((vector<char*> *)(input.at(1).p))->size()!=7)
        {
            v3d_msg(tr("ERROR: no enough para!"));
//            printHelp();
            return false;
        }

//        get paras
        QString qs_filename_img_sub=((vector<char*> *)(input.at(0).p))->at(0);
        QString qs_filename_marker_sub=((vector<char*> *)(input.at(0).p))->at(1);
        QString qs_filename_marker_tar=((vector<char*> *)(input.at(0).p))->at(2);
        QString qs_filename_img_warp=((vector<char*> *)(output.at(0).p))->at(0);
        vector<char*> paras = (*(vector<char*> *)(input.at(1).p));
        int program_mode=atoi(paras.at(0));
        int device_mode=atoi(paras.at(1));
        int downsample_mode = atoi(paras.at(2));

        if(!bigimagewarp(callback,  qs_filename_img_sub,qs_filename_marker_sub,qs_filename_marker_tar,
                         1, 0, device_mode, program_mode, downsample_mode, qs_filename_img_warp))
        {
            v3d_msg(tr("ERROR: bigimagewarp() return false!"));
            return false;
        }
        return true;
    }
    else
    {
        printHelp();
    }
}

void printHelp()
{
    printf("\nUsage: v3d -x <bigimagewarp> -f bigimagewarp -i <input_image_sub> <input_marker_sub> <input_marker_tar> -o <output_image_file> -p program_mode device_mode downsample_mode\n");
    printf("\t input_image_sub:         input image file or folder to be warped (subject image)\n");
    printf("\t input_marker_sub:        marker file in subject image (markers in this file define the control points in subject image)\n");
    printf("\t input_marker_tar:        marker file in target image (markers in this file define the control points in target image)\n");
    printf("\t output_image_file:       output warped image file or folder(warped image)\n");
    printf("\t program_mode:         image warp method (0:affine, 1:stps)\n");
    printf("\t device_mode:        device mode (0:cpu, 1:gpu)\n");
    printf("Demo :\t v3d -x bigimagewarp -f bigimagewarp -i /Users/qul/Desktop/testdata/output_global.v3draw /Users/qul/Desktop/testdata/output_subject.marker /Users/qul/Desktop/testdata/output_target.marker -o /Users/qul/Desktop/testdata/output_warp_littlequick.v3draw -p 1 0\n");
    return;
}


bool bigimagewarp(V3DPluginCallback2 & callback,  QString qs_filename_img_sub,QString qs_filename_tpscpt_sub,QString qs_filename_tpscpt_tar,
                  int i_interpmethod_df,int i_interpmethod_img, int GPU_TYPE, int i_program_mode, int i_downsample, QString qs_filename_img_output)
{
    printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
    printf(">>Big image warping:\n");
    printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
    printf(">>input parameters:\n");
    printf(">>  input subject image:           %s\n",qPrintable(qs_filename_img_sub));
    printf(">>  input target  marker:          %s\n",qPrintable(qs_filename_tpscpt_tar));
    printf(">>  input subject marker:          %s\n",qPrintable(qs_filename_tpscpt_sub));
    printf(">>  Program run mode :             %d\n",i_program_mode);
    printf(">>  device mode      :             %d\n",GPU_TYPE);
//    printf(">>  DF  interp method:             %d\n",i_interpmethod_df);
//    printf(">>  img interp method:             %d\n",i_interpmethod_img);
    printf(">>-------------------------\n");
    printf(">>output parameters:\n");


    QList<ImageMarker> ql_marker_tar, ql_marker_sub;
    printf("1.  Read marker files. \n");
    if (qs_filename_tpscpt_tar.endsWith(".marker") && qs_filename_tpscpt_sub.endsWith(".marker"))
    {
        ql_marker_tar = readMarker_file(qs_filename_tpscpt_tar);
        ql_marker_sub = readMarker_file(qs_filename_tpscpt_sub);
        printf("\t>>Target: read %d markers from [%s]\n", ql_marker_tar.size(), qPrintable(qs_filename_tpscpt_tar));
        printf("\t>>Subject:read %d markers from [%s]\n", ql_marker_sub.size(), qPrintable(qs_filename_tpscpt_sub));
    }
    else
    {
        printf("ERROR: at least one marker file is invalid.\n");
        return false;
    }


    printf("2. Read input subject files. \n");
    printf("2-1. Read subject image file. \n");
    // 获取sub路径和文件名
    if (qs_filename_img_sub.endsWith(".v3draw"))
    {
        string s = qs_filename_img_sub.toStdString();
        s = s.substr(s.find_last_of('/') + 1, s.find_last_of('.') - s.find_last_of('/') - 1);
        string TeraWarp_save_base;

        if (i_program_mode)
        {
            TeraWarp_save_base = qs_filename_img_output.toStdString() + '/' + s + "_affine" + ".v3draw";
            printf(">>  output warped image folder:           %s\n", TeraWarp_save_base.c_str());
        }
        else
        {
            TeraWarp_save_base = qs_filename_img_output.toStdString() + '/' + s + "_stps" + ".v3draw";
            printf(">>  output warped image folder:           %s\n", TeraWarp_save_base.c_str());
        }

        unsigned char *p_img_sub = 0;
        V3DLONG *sz_img_sub_ = 0;
        int b_swap;

        loadImage((char *)qPrintable(qs_filename_img_sub), p_img_sub, sz_img_sub_, b_swap);

        long long sz_img_sub[4] = {sz_img_sub_[0], sz_img_sub_[1], sz_img_sub_[2], sz_img_sub_[3]};

        if (i_program_mode)
        {
            if (!q_imagewarp_affine_image(ql_marker_tar, ql_marker_sub, p_img_sub, sz_img_sub, TeraWarp_save_base, GPU_TYPE))
            {
                printf("ERROR: q_imagewarp_affine return false!\n");
                if (p_img_sub) 				{ delete[]p_img_sub;			p_img_sub = 0; }
//                if (sz_img_sub) 			{ delete[]sz_img_sub;			sz_img_sub = 0; }
                return false;
            }
        }
        else
        {
            if (!(imgwarp_smallmemory_image(ql_marker_sub, ql_marker_tar, 4, 4, 4, i_interpmethod_df, i_interpmethod_img, p_img_sub, sz_img_sub_, TeraWarp_save_base, GPU_TYPE)))
            {
                printf("ERROR:imgwarp_smallmemory() return false. \n");
                return false;
            }
        }

        if (p_img_sub) {delete[] p_img_sub; p_img_sub = 0;}





    }
    else
    {
        QByteArray b_qs_filename_img_sub = qs_filename_img_sub.toLatin1();
        char *imgSrcFile = b_qs_filename_img_sub.data();
        // 获取输出路径
        string imgSrcFile_str = qs_filename_img_sub.toStdString();
        string TeraWarp_save_base;
#if defined(_MSC_VER) && (_WIN64)
        int pos = imgSrcFile_str.find_last_of('\\') + 1;
        string filename = imgSrcFile_str.substr(pos, imgSrcFile_str.length() - pos);
        TeraWarp_save_base = qs_filename_img_output.toStdString() + "\\" + filename;
#else
        int pos = imgSrcFile_str.find_last_of('/') + 1;
        string filename = imgSrcFile_str.substr(pos, imgSrcFile_str.length() - pos);
        TeraWarp_save_base = qs_filename_img_output.toStdString() + "/" + filename;
#endif
        if (i_program_mode)
            printf(">>  output warped image folder:           %s%s\n", TeraWarp_save_base.c_str(),"_affine_slice");
        else
            printf(">>  output warped image folder:           %s%s\n", TeraWarp_save_base.c_str(), "_stps_slice");
        printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");



        int b_swap;
        int header_len;
        int tmp_datatype = 0;
        char *err_rawfmt;
        void *fhandle;
        string *path_TIF = 0;
        V3DLONG *sz = 0;
        long long Z_DEPTH_MAX;

        if (!read_2Dtif_dir(imgSrcFile, path_TIF, Z_DEPTH_MAX))//读取TIF文件目录列表放入path_TIF数组
        {
            printf("ERROR: read_2Dtif_dir() return false in read [%s].\n", imgSrcFile);
            return false;
        }
        if ((err_rawfmt = loadMetadata_cyf((char*)path_TIF[0].c_str(), sz, tmp_datatype, b_swap, header_len, fhandle)) != 0) {
            if (sz) delete[] sz;
            printf("ERROR: read_2Dtif_dir() return false in read [%s].\n", imgSrcFile);
            return false;
        }
        long long sz_img_ori[4] = {sz[0], sz[1], Z_DEPTH_MAX, 1};

        if (i_program_mode)
        {
            if (!q_imagewarp_affine_cyf(ql_marker_tar, ql_marker_sub, path_TIF, sz_img_ori, TeraWarp_save_base.c_str(), GPU_TYPE))
            {
                printf("ERROR: q_imagewarp_affine_cyf() return false!\n");
                return false;
            }

            if (i_downsample)
            {
                string Tera = TeraWarp_save_base + "_affine_slice";
                QString s_Tera = QString::fromStdString(Tera);
                QByteArray b_Tera = s_Tera.toLatin1();
                char *imgdownfile = b_Tera.data();
                string *path_TIF_downsample = 0;
                long long Z_DEPTH_MAX_downsample  = 0;
                TeraWarp_save_base += "_affine";

                if (!read_2Dtif_dir(imgdownfile, path_TIF_downsample, Z_DEPTH_MAX_downsample))//读取TIF文件目录列表放入path_TIF数组
                {
                    printf("ERROR: read_2Dtif_dir() return false in read [%s].\n", imgdownfile);
                    return false;
                }

                if (!(q_imagewarp_affine_cyf_downsample(ql_marker_tar, ql_marker_tar, path_TIF_downsample, sz_img_ori, TeraWarp_save_base.c_str(), GPU_TYPE)))
                {
                    printf("ERROR: downsample() return false!\n");
                    if (path_TIF_downsample) {delete[]path_TIF_downsample;  path_TIF_downsample = 0; }
                    return false;

                }

                if (path_TIF_downsample) {delete[]path_TIF_downsample;  path_TIF_downsample = 0; }
            }

        }
        else
        {
            if (!(imgwarp_smallmemory_CYF(ql_marker_sub, ql_marker_tar, 4, 4, 4, i_interpmethod_df, i_interpmethod_img, path_TIF, sz_img_ori, TeraWarp_save_base.c_str(), GPU_TYPE)))
            {
                printf("ERROR:imgwarp_smallmemory_CYF() return false. \n");
                return false;
            }

            if (i_downsample)
            {
                string Tera = TeraWarp_save_base + "_stps_slice";
                QString s_Tera = QString::fromStdString(Tera);
                QByteArray b_Tera = s_Tera.toLatin1();
                char *imgdownfile = b_Tera.data();
                string *path_TIF_downsample = 0;
                long long Z_DEPTH_MAX_downsample  = 0;
                TeraWarp_save_base += "_stps";

                if (!read_2Dtif_dir(imgdownfile, path_TIF_downsample, Z_DEPTH_MAX_downsample))//读取TIF文件目录列表放入path_TIF数组
                {
                    printf("ERROR: read_2Dtif_dir() return false in read [%s].\n", imgdownfile);
                    return false;
                }

                if (!(q_imagewarp_affine_cyf_downsample(ql_marker_tar, ql_marker_tar, path_TIF_downsample, sz_img_ori, TeraWarp_save_base.c_str(), GPU_TYPE)))
                {
                    printf("ERROR: downsample() return false!\n");
                    if (path_TIF_downsample) {delete[]path_TIF_downsample;  path_TIF_downsample = 0; }
                    return false;

                }

                if (path_TIF_downsample) {delete[]path_TIF_downsample;  path_TIF_downsample = 0; }
            }

        }

        if (path_TIF) {delete[]path_TIF;  path_TIF = 0; }

    }


    return true;

}












