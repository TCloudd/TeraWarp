//q_paradialog_littlequickwarp.cpp
//by Lei Qu
//2011-04-11

#include <QtGui>
#include "q_paradialog_bigimagewarp.h"


CParaDialog_bigimagewarp::CParaDialog_bigimagewarp(V3DPluginCallback &callback,QWidget *parent):QDialog(parent)
{
	setupUi(this);
	IniDialog();
}
CParaDialog_bigimagewarp::~CParaDialog_bigimagewarp()
{
	QSettings settings("V3D plugin","Bigimagewarp");

	settings.setValue("img_sub",this->lineEdit_img_sub->text());
	settings.setValue("img_warp",this->lineEdit_img_warp->text());
	settings.setValue("marker_sub",this->lineEdit_marker_sub->text());
	settings.setValue("marker_tar",this->lineEdit_marker_tar->text());


//	settings.setValue("dfinterpmethod",this->radioButton_df_interp_bspline->isChecked());
//	settings.setValue("imginterpmethod",this->radioButton_img_interp_nn->isChecked());

    settings.setValue("devicemode",this->radioButton_device_gpu->isChecked());
    settings.setValue("programerunmode",this->radioButton_affine->isChecked());
}

void CParaDialog_bigimagewarp::IniDialog()
{
	//read settings
	QSettings settings("V3D plugin","Bigimagewarp");
	this->lineEdit_img_sub->setText(settings.value("img_sub").toString());
	this->lineEdit_img_warp->setText(settings.value("img_warp").toString());
	this->lineEdit_marker_sub->setText(settings.value("marker_sub").toString());
	this->lineEdit_marker_tar->setText(settings.value("marker_tar").toString());

//	this->radioButton_df_interp_linear->setChecked(!(settings.value("dfinterpmethod",1).toBool()));
//	this->radioButton_df_interp_bspline->setChecked(settings.value("dfinterpmethod",0).toBool());
//	this->radioButton_img_interp_linear->setChecked(!(settings.value("imginterpmethod",0).toBool()));
//	this->radioButton_img_interp_nn->setChecked(settings.value("imginterpmethod",1).toBool());

    this->radioButton_device_gpu->setChecked(settings.value("devicemode",1).toBool());
    this->radioButton_device_cpu->setChecked(!(settings.value("devicemode",0).toBool()));
    this->radioButton_affine->setChecked(settings.value("programerunmode",0).toBool());
    this->radioButton_stps->setChecked(!(settings.value("programerunmode",1).toBool()));


	connect(pushButton_img_sub,SIGNAL(clicked()),this,SLOT(_slots_openDlg_img_sub()));
	connect(pushButton_img_warp,SIGNAL(clicked()),this,SLOT(_slots_openDlg_img_warp()));
	connect(pushButton_marker_sub,SIGNAL(clicked()),this,SLOT(_slots_openDlg_marker_sub()));
	connect(pushButton_marker_tar,SIGNAL(clicked()),this,SLOT(_slots_openDlg_marker_tar()));
//    connect(comboBox_input, SIGNAL(currentIndexChanged(int)), this, SLOT(on_comboBox_input_currentIndexChanged(int)));

}

void CParaDialog_bigimagewarp::_slots_openDlg_img_sub()
{
//    QString dirpath_sub = QFileDialog::getExistingDirectory(this,"Choose subject image folder","C:",QFileDialog::ShowDirsOnly);
//    this->lineEdit_img_sub->setText(dirpath_sub);
    if (this->comboBox_input->currentIndex() == 0)
    {
        QString dirpath_sub = QFileDialog::getExistingDirectory(this,"Choose subject image folder","C:",QFileDialog::ShowDirsOnly);
        this->lineEdit_img_sub->setText(dirpath_sub);
    }
    else if (this->comboBox_input->currentIndex() == 1)
    {
        QString dirpath_sub = QFileDialog::getOpenFileName(this, "Choose subject image file", "C:", "V3Draw file (*.v3draw)");
        this->lineEdit_img_sub->setText(dirpath_sub);
    }
}
void CParaDialog_bigimagewarp::_slots_openDlg_marker_sub()
{
    QString marker_sub = QFileDialog::getOpenFileName(this, "Choose subject marker file", "C:", "Marker file (*.marker)");
    this->lineEdit_marker_sub->setText(marker_sub);
}
void CParaDialog_bigimagewarp::_slots_openDlg_marker_tar()
{
    QString marker_tar = QFileDialog::getOpenFileName(this, "Choose target marker file", "C:", "Marker file (*.marker)");
    this->lineEdit_marker_tar->setText(marker_tar);
}
void CParaDialog_bigimagewarp::_slots_openDlg_img_warp()
{
    QString dirpath_warp = QFileDialog::getExistingDirectory(this,"Choose warpped image folder","C:",QFileDialog::ShowDirsOnly);
    this->lineEdit_img_warp->setText(dirpath_warp);
}

//void CParaDialog_bigimagewarp::on_comboBox_input_currentIndexChanged(int index)
//{
//    if (index == 0)
//    {
//        this->groupBox_5->setEnabled(true);
//    }
//    else
//    {
//        this->groupBox_5->setEnabled(false);
//    }
//}
