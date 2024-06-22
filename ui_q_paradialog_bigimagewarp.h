/********************************************************************************
** Form generated from reading UI file 'q_paradialog_bigimagewarp.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_Q_PARADIALOG_BIGIMAGEWARP_H
#define UI_Q_PARADIALOG_BIGIMAGEWARP_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QComboBox>
#include <QtGui/QDialog>
#include <QtGui/QDialogButtonBox>
#include <QtGui/QGroupBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>

QT_BEGIN_NAMESPACE

class Ui_Paradialog_bigimagewarp
{
public:
    QDialogButtonBox *buttonBox;
    QGroupBox *groupBox_fileio;
    QLabel *label_img_sub_2;
    QLabel *label_img_sub;
    QLineEdit *lineEdit_img_sub;
    QPushButton *pushButton_img_sub;
    QPushButton *pushButton_marker_sub;
    QLineEdit *lineEdit_marker_sub;
    QLineEdit *lineEdit_marker_tar;
    QLabel *label_img_sub2tar;
    QPushButton *pushButton_marker_tar;
    QPushButton *pushButton_img_warp;
    QLabel *label_swc_grid;
    QLineEdit *lineEdit_img_warp;
    QComboBox *comboBox_input;
    QGroupBox *groupBox_paras;
    QGroupBox *groupBox_3;
    QRadioButton *radioButton_device_cpu;
    QRadioButton *radioButton_device_gpu;
    QGroupBox *groupBox_4;
    QRadioButton *radioButton_affine;
    QRadioButton *radioButton_stps;

    void setupUi(QDialog *Paradialog_bigimagewarp)
    {
        if (Paradialog_bigimagewarp->objectName().isEmpty())
            Paradialog_bigimagewarp->setObjectName(QString::fromUtf8("Paradialog_bigimagewarp"));
        Paradialog_bigimagewarp->resize(704, 497);
        buttonBox = new QDialogButtonBox(Paradialog_bigimagewarp);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setGeometry(QRect(530, 460, 171, 32));
        buttonBox->setStyleSheet(QString::fromUtf8(""));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);
        groupBox_fileio = new QGroupBox(Paradialog_bigimagewarp);
        groupBox_fileio->setObjectName(QString::fromUtf8("groupBox_fileio"));
        groupBox_fileio->setGeometry(QRect(10, 0, 691, 271));
        label_img_sub_2 = new QLabel(groupBox_fileio);
        label_img_sub_2->setObjectName(QString::fromUtf8("label_img_sub_2"));
        label_img_sub_2->setGeometry(QRect(10, 30, 331, 20));
        label_img_sub = new QLabel(groupBox_fileio);
        label_img_sub->setObjectName(QString::fromUtf8("label_img_sub"));
        label_img_sub->setGeometry(QRect(10, 90, 161, 20));
        lineEdit_img_sub = new QLineEdit(groupBox_fileio);
        lineEdit_img_sub->setObjectName(QString::fromUtf8("lineEdit_img_sub"));
        lineEdit_img_sub->setGeometry(QRect(200, 50, 341, 22));
        pushButton_img_sub = new QPushButton(groupBox_fileio);
        pushButton_img_sub->setObjectName(QString::fromUtf8("pushButton_img_sub"));
        pushButton_img_sub->setGeometry(QRect(550, 50, 131, 21));
        pushButton_marker_sub = new QPushButton(groupBox_fileio);
        pushButton_marker_sub->setObjectName(QString::fromUtf8("pushButton_marker_sub"));
        pushButton_marker_sub->setGeometry(QRect(360, 110, 131, 21));
        lineEdit_marker_sub = new QLineEdit(groupBox_fileio);
        lineEdit_marker_sub->setObjectName(QString::fromUtf8("lineEdit_marker_sub"));
        lineEdit_marker_sub->setGeometry(QRect(10, 110, 341, 22));
        lineEdit_marker_tar = new QLineEdit(groupBox_fileio);
        lineEdit_marker_tar->setObjectName(QString::fromUtf8("lineEdit_marker_tar"));
        lineEdit_marker_tar->setGeometry(QRect(10, 170, 341, 22));
        label_img_sub2tar = new QLabel(groupBox_fileio);
        label_img_sub2tar->setObjectName(QString::fromUtf8("label_img_sub2tar"));
        label_img_sub2tar->setGeometry(QRect(10, 150, 161, 21));
        pushButton_marker_tar = new QPushButton(groupBox_fileio);
        pushButton_marker_tar->setObjectName(QString::fromUtf8("pushButton_marker_tar"));
        pushButton_marker_tar->setGeometry(QRect(360, 170, 131, 21));
        pushButton_img_warp = new QPushButton(groupBox_fileio);
        pushButton_img_warp->setObjectName(QString::fromUtf8("pushButton_img_warp"));
        pushButton_img_warp->setGeometry(QRect(360, 230, 131, 20));
        label_swc_grid = new QLabel(groupBox_fileio);
        label_swc_grid->setObjectName(QString::fromUtf8("label_swc_grid"));
        label_swc_grid->setGeometry(QRect(10, 210, 241, 20));
        lineEdit_img_warp = new QLineEdit(groupBox_fileio);
        lineEdit_img_warp->setObjectName(QString::fromUtf8("lineEdit_img_warp"));
        lineEdit_img_warp->setGeometry(QRect(10, 230, 341, 22));
        comboBox_input = new QComboBox(groupBox_fileio);
        comboBox_input->setObjectName(QString::fromUtf8("comboBox_input"));
        comboBox_input->setGeometry(QRect(10, 50, 161, 22));
        groupBox_paras = new QGroupBox(Paradialog_bigimagewarp);
        groupBox_paras->setObjectName(QString::fromUtf8("groupBox_paras"));
        groupBox_paras->setGeometry(QRect(10, 300, 691, 141));
        groupBox_3 = new QGroupBox(groupBox_paras);
        groupBox_3->setObjectName(QString::fromUtf8("groupBox_3"));
        groupBox_3->setGeometry(QRect(360, 30, 321, 81));
        radioButton_device_cpu = new QRadioButton(groupBox_3);
        radioButton_device_cpu->setObjectName(QString::fromUtf8("radioButton_device_cpu"));
        radioButton_device_cpu->setGeometry(QRect(10, 20, 91, 18));
        radioButton_device_cpu->setChecked(true);
        radioButton_device_gpu = new QRadioButton(groupBox_3);
        radioButton_device_gpu->setObjectName(QString::fromUtf8("radioButton_device_gpu"));
        radioButton_device_gpu->setGeometry(QRect(10, 50, 91, 18));
        groupBox_4 = new QGroupBox(groupBox_paras);
        groupBox_4->setObjectName(QString::fromUtf8("groupBox_4"));
        groupBox_4->setGeometry(QRect(10, 30, 321, 80));
        radioButton_affine = new QRadioButton(groupBox_4);
        radioButton_affine->setObjectName(QString::fromUtf8("radioButton_affine"));
        radioButton_affine->setGeometry(QRect(10, 20, 91, 20));
        radioButton_affine->setChecked(true);
        radioButton_stps = new QRadioButton(groupBox_4);
        radioButton_stps->setObjectName(QString::fromUtf8("radioButton_stps"));
        radioButton_stps->setGeometry(QRect(10, 50, 91, 18));
        groupBox_paras->raise();
        buttonBox->raise();
        groupBox_fileio->raise();

        retranslateUi(Paradialog_bigimagewarp);
        QObject::connect(buttonBox, SIGNAL(accepted()), Paradialog_bigimagewarp, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), Paradialog_bigimagewarp, SLOT(reject()));

        QMetaObject::connectSlotsByName(Paradialog_bigimagewarp);
    } // setupUi

    void retranslateUi(QDialog *Paradialog_bigimagewarp)
    {
        Paradialog_bigimagewarp->setWindowTitle(QApplication::translate("Paradialog_bigimagewarp", "BigImageWarp v1.0", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        buttonBox->setToolTip(QApplication::translate("Paradialog_bigimagewarp", "<html><head/><body><p><br/></p></body></html>", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        groupBox_fileio->setTitle(QApplication::translate("Paradialog_bigimagewarp", "File IO: ", 0, QApplication::UnicodeUTF8));
        label_img_sub_2->setText(QApplication::translate("Paradialog_bigimagewarp", "Input subject image format and folder(to be warpped): ", 0, QApplication::UnicodeUTF8));
        label_img_sub->setText(QApplication::translate("Paradialog_bigimagewarp", "Input subject marker file:", 0, QApplication::UnicodeUTF8));
        pushButton_img_sub->setText(QApplication::translate("Paradialog_bigimagewarp", "Browse for dir...", 0, QApplication::UnicodeUTF8));
        pushButton_marker_sub->setText(QApplication::translate("Paradialog_bigimagewarp", "Browse for file...", 0, QApplication::UnicodeUTF8));
        label_img_sub2tar->setText(QApplication::translate("Paradialog_bigimagewarp", "Input target marker file:", 0, QApplication::UnicodeUTF8));
        pushButton_marker_tar->setText(QApplication::translate("Paradialog_bigimagewarp", "Browse for file...", 0, QApplication::UnicodeUTF8));
        pushButton_img_warp->setText(QApplication::translate("Paradialog_bigimagewarp", "Browse for dir...", 0, QApplication::UnicodeUTF8));
        label_swc_grid->setText(QApplication::translate("Paradialog_bigimagewarp", "Output warped image folder:", 0, QApplication::UnicodeUTF8));
        lineEdit_img_warp->setText(QString());
        comboBox_input->clear();
        comboBox_input->insertItems(0, QStringList()
         << QApplication::translate("Paradialog_bigimagewarp", "Vaa3D raw\357\274\210series\357\274\2142D\357\274\211", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("Paradialog_bigimagewarp", "Vaa3D raw", 0, QApplication::UnicodeUTF8)
        );
        groupBox_paras->setTitle(QApplication::translate("Paradialog_bigimagewarp", "Paras:", 0, QApplication::UnicodeUTF8));
        groupBox_3->setTitle(QApplication::translate("Paradialog_bigimagewarp", "Device mode:", 0, QApplication::UnicodeUTF8));
        radioButton_device_cpu->setText(QApplication::translate("Paradialog_bigimagewarp", "CPU", 0, QApplication::UnicodeUTF8));
        radioButton_device_gpu->setText(QApplication::translate("Paradialog_bigimagewarp", "GPU", 0, QApplication::UnicodeUTF8));
        groupBox_4->setTitle(QApplication::translate("Paradialog_bigimagewarp", " Program run mode:", 0, QApplication::UnicodeUTF8));
        radioButton_affine->setText(QApplication::translate("Paradialog_bigimagewarp", "Affine", 0, QApplication::UnicodeUTF8));
        radioButton_stps->setText(QApplication::translate("Paradialog_bigimagewarp", "Non-rigid", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class Paradialog_bigimagewarp: public Ui_Paradialog_bigimagewarp {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_Q_PARADIALOG_BIGIMAGEWARP_H
