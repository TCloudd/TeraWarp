//q_paradialog_littlequickwarp.h
//by Lei Qu
//2012-07-16

#ifndef __Q_PARADIALOG_BIGIMAGEWARP_H__
#define __Q_PARADIALOG_BIGIMAGEWARP_H__

#include <QDialog>
#include <v3d_interface.h>
#include "ui_q_paradialog_bigimagewarp.h"


class CParaDialog_bigimagewarp : public QDialog, public Ui::Paradialog_bigimagewarp
{
	Q_OBJECT

public:
	CParaDialog_bigimagewarp(V3DPluginCallback &callback,QWidget *parent);
	~CParaDialog_bigimagewarp();

private slots:
	void _slots_openDlg_img_sub();
	void _slots_openDlg_img_warp();
	void _slots_openDlg_marker_sub();
	void _slots_openDlg_marker_tar();
//    void on_comboBox_input_currentIndexChanged(int index);
//    void _slots_openDlg_img_med();

public:
	void IniDialog();

	v3dhandleList h_wndlist;

};


#endif
