/*
* Author: Frederic GARCIA BECERRO
* Email: frederic.garcia.becerro@gmail.com
* Website: http://www.frederic-garcia-becerro.com
*/

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

// Std
#include <fstream>
#include <string>
// QT
#include <QMainWindow>
#include <QFileDialog>
#include <QMessageBox>
// OpenCV
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/contrib/contrib.hpp"
// Include Data Fusion class
#include "c_datafusion.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_ui_sl_sigmaSpatial_valueChanged(int value);
    void on_ui_sl_sigmaSpatial_sliderReleased();
    void on_ui_sl_sigmaRange_valueChanged(int value);
    void on_ui_sl_sigmaRange_sliderReleased();
    void on_ui_sl_sigmaCredMap_valueChanged(int value);
    void on_ui_sl_sigmaCredMap_sliderReleased();
    void on_ui_sl_maxPalette_valueChanged(int value);
    void on_ui_sl_maxPalette_sliderReleased();
    void on_ui_sl_minPalette_valueChanged(int value);
    void on_ui_sl_minPalette_sliderReleased();
    void on_ui_pb_LoadDepthMap_clicked();
    void on_ui_pb_LoadGuidImage_clicked();
    void on_ui_pb_ApplyFilter_clicked();
    void on_ui_cb_dataFusion_Display_currentIndexChanged(int index);
    void on_ui_pb_SaveDisplayedImage_clicked();
    void on_ui_ActionAbout_clicked();

private:
    Ui::MainWindow *ui;

    c_DataFusion *m_DataFusion; // Instance to Data Fusion class

    cv::Mat *m_cvDepthImage; // Depth map
    cv::Mat *m_cvRGBImage; // Guidance image
    bool    m_bDepthImageLoaded; // Flag to indicate loaded depth map
    bool    m_bRGBImageLoaded; // Flag to indicate loaded RGB image

    short   m_sXRes; // X Resolution
    short   m_sYRes; // Y Resolution
    short   m_sMaxPaletteValue; // Max palette value
    short   m_sMinPaletteValue; // Min palette value
};

#endif // MAINWINDOW_H
