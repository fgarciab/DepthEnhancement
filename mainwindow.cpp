/*
* Author: Frederic GARCIA BECERRO
* Email: frederic.garcia.becerro@gmail.com
* Website: http://www.frederic-garcia-becerro.com
*/

#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->statusBar->showMessage(tr("Frederic Garcia - http://www.frederic-garcia-becerro.com - frederic.garcia.becerro@gmail.com"));
    connect(ui->actionExit, SIGNAL(triggered()), this, SLOT(close()));
    connect(ui->actionLoad_Guid_Image_I, SIGNAL(triggered()), this, SLOT(on_ui_pb_LoadGuidImage_clicked()));
    connect(ui->actionLoad_Depth_Map_D, SIGNAL(triggered()), this, SLOT(on_ui_pb_LoadDepthMap_clicked()));
    connect(ui->actionSave_Displayed_Image, SIGNAL(triggered()), this, SLOT(on_ui_pb_SaveDisplayedImage_clicked()));    

    m_sMaxPaletteValue = ui->ui_sl_maxPalette->value()*10;
    m_sMinPaletteValue = ui->ui_sl_minPalette->value()*10;
    m_sXRes = 640; // Default res
    m_sYRes = 480;
    m_sScaleFactor = ui->ui_cb_dataFusion_Downsampling->currentIndex();

    // Instance to the Data Fusion class, downsampling factor
    m_DataFusion = new c_DataFusion(D_Fusion_2D_3D_2DW, D_Fusion_2D_3D_2DH, m_sScaleFactor); // D_Fusion_2D_3D_2DW = 640; D_Fusion_2D_3D_2DH = 480;
    // Initialize depth map
    m_cvDepthImage = new cv::Mat(m_sYRes, m_sXRes, CV_16UC1, cv::Scalar(0));
    m_cvRGBImage = new cv::Mat(m_sYRes, m_sXRes, CV_8UC1, cv::Scalar(0));
    m_bDepthImageLoaded = false;
    m_bRGBImageLoaded = false;
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_ui_sl_sigmaSpatial_valueChanged(int value)
{
    QString l_qStr = "Sigma Spatial: ";

    if (!ui->ui_sl_sigmaSpatial->isSliderDown())
    {
        l_qStr.append(QString::number(value));
        ui->ui_lbl_sigmaSpatial->setText(l_qStr);
        m_DataFusion->SetSigmaS(value);
    }
}

void MainWindow::on_ui_sl_sigmaSpatial_sliderReleased()
{
    QString l_qStr = "Sigma Spatial: ";

    l_qStr.append(QString::number(ui->ui_sl_sigmaSpatial->value()));
    ui->ui_lbl_sigmaSpatial->setText(l_qStr);
    m_DataFusion->SetSigmaS(ui->ui_sl_sigmaSpatial->value());
}

void MainWindow::on_ui_sl_sigmaRange_valueChanged(int value)
{
    QString l_qStr = "Sigma Range: ";

    if (!ui->ui_sl_sigmaRange->isSliderDown())
    {
        l_qStr.append(QString::number(value));
        ui->ui_lbl_sigmaRange->setText(l_qStr);
        m_DataFusion->SetSigmaR(value);
    }
}

void MainWindow::on_ui_sl_sigmaRange_sliderReleased()
{
    QString l_qStr = "Sigma Range: ";

    l_qStr.append(QString::number(ui->ui_sl_sigmaRange->value()));
    ui->ui_lbl_sigmaRange->setText(l_qStr);
    m_DataFusion->SetSigmaR(ui->ui_sl_sigmaRange->value());
}

void MainWindow::on_ui_sl_sigmaCredMap_valueChanged(int value)
{
    QString l_qStr = "Sigma Credibility Map: ";

    if (!ui->ui_sl_sigmaCredMap->isSliderDown())
    {
        l_qStr.append(QString::number(value));
        ui->ui_lbl_sigmaCredMap->setText(l_qStr);
        m_DataFusion->SetSigmaQ(value);
    }
}

void MainWindow::on_ui_sl_sigmaCredMap_sliderReleased()
{
    QString l_qStr = "Sigma Credibility Map: ";

    l_qStr.append(QString::number(ui->ui_sl_sigmaCredMap->value()));
    ui->ui_lbl_sigmaCredMap->setText(l_qStr);
    m_DataFusion->SetSigmaQ(ui->ui_sl_sigmaCredMap->value());
}

void MainWindow::on_ui_sl_maxPalette_valueChanged(int value)
{
    QString l_qStr = "Max palette: ";

    if (!ui->ui_sl_maxPalette->isSliderDown())
    {
        l_qStr.append(QString::number(value*10));
        l_qStr.append(" mm");
        ui->ui_lbl_maxPalette->setText(l_qStr);
        m_sMaxPaletteValue = value*10;
        // Update GUI
        on_ui_cb_dataFusion_Display_currentIndexChanged(ui->ui_cb_dataFusion_Display->currentIndex());
        on_ui_cb_dataFusion_Display_currentIndexChanged(ui->ui_cb_dataFusion_Display->currentIndex());
    }
}

void MainWindow::on_ui_sl_maxPalette_sliderReleased()
{
    QString l_qStr = "Max palette: ";

    l_qStr.append(QString::number(ui->ui_sl_maxPalette->value()*10));
    l_qStr.append(" mm");
    ui->ui_lbl_maxPalette->setText(l_qStr);
    m_sMaxPaletteValue = ui->ui_sl_maxPalette->value()*10;
    // Update GUI
    on_ui_cb_dataFusion_Display_currentIndexChanged(ui->ui_cb_dataFusion_Display->currentIndex());
}

void MainWindow::on_ui_sl_minPalette_valueChanged(int value)
{
    QString l_qStr = "Min palette: ";

    if (!ui->ui_sl_minPalette->isSliderDown())
    {
        l_qStr.append(QString::number(value*10));
        l_qStr.append(" mm");
        ui->ui_lbl_minPalette->setText(l_qStr);
        m_sMinPaletteValue = value*10;
        // Update GUI
        on_ui_cb_dataFusion_Display_currentIndexChanged(ui->ui_cb_dataFusion_Display->currentIndex());
    }
}

void MainWindow::on_ui_sl_minPalette_sliderReleased()
{
    QString l_qStr = "Min palette: ";

    l_qStr.append(QString::number(ui->ui_sl_minPalette->value()*10));
    l_qStr.append(" mm");
    ui->ui_lbl_minPalette->setText(l_qStr);
    m_sMinPaletteValue = ui->ui_sl_minPalette->value()*10;
    // Update GUI
    on_ui_cb_dataFusion_Display_currentIndexChanged(ui->ui_cb_dataFusion_Display->currentIndex());
}

void MainWindow::on_ui_pb_LoadDepthMap_clicked()
{
    QFileInfo l_qFileInfo(QFileDialog::getOpenFileName(this, tr("Load depth map"), "./", tr("Depth map Files (*.bmp *.png *.pgm *.txt)")));
    if (!l_qFileInfo.fileName().isEmpty()) // Go ahead
    {
        QString l_QStrFilename = l_qFileInfo.absoluteFilePath();
        // Depth data to be filtered
        unsigned short* l_pBuffer = (unsigned short*)m_cvDepthImage->data;        
        if (l_QStrFilename.contains(QString(".txt")))
        {
            std::ifstream l_fp_in;  // Declaration of input stream
            l_fp_in.open(l_QStrFilename.toUtf8().constData(), std::ios::in);  // Open the input file stream
            if (l_fp_in.is_open())
            {
                unsigned short l_usAux;

                for (short y=0; y<m_sYRes; y++)
                {
                    for (short x=0; x<m_sXRes; x++)
                    {
                        l_fp_in >> l_usAux; // Depth
                        *l_pBuffer++ = l_usAux;
                    }
                }
                l_fp_in.close();   // close the streams
            }            
        }
        else
        {
            *m_cvDepthImage = cv::imread(l_QStrFilename.toUtf8().constData(), -1);
        }
        // Check image dimensions
        if ((m_cvDepthImage->cols != m_sXRes) || (m_cvDepthImage->rows != m_sYRes))
        {
            QMessageBox::information(0, QString("Warning"), QString(QString("Depth map resolution must be (")+
                                                                    QString::number(m_sXRes)+QString("x")+
                                                                    QString::number(m_sYRes)+QString(")\nDepth map has been resized!")));
            cv::resize(*m_cvDepthImage, *m_cvDepthImage, cv::Size(m_sXRes, m_sYRes), 0, 0, cv::INTER_NEAREST);
        }
        // Display depth map
        QImage *l_pQImage;
        QPixmap l_QPixmap;
        cv::Mat l_cvDepthImage = m_cvDepthImage->clone();

        // Note: the followind 16-bit to 8-bit conversion not only chop off the high order byte instead but it scales the values to an 8-bit range
        l_cvDepthImage.convertTo(l_cvDepthImage, CV_8U, 255.0/(m_sMaxPaletteValue-m_sMinPaletteValue), -m_sMinPaletteValue*255.0/(m_sMaxPaletteValue-m_sMinPaletteValue));
        if (l_cvDepthImage.channels() == 1)
            cv::cvtColor(l_cvDepthImage, l_cvDepthImage, CV_GRAY2RGB);
        cv::applyColorMap(l_cvDepthImage, l_cvDepthImage, cv::COLORMAP_JET);
        l_pQImage = new QImage((uchar*)l_cvDepthImage.data, l_cvDepthImage.cols, l_cvDepthImage.rows, QImage::Format_RGB888);
        l_QPixmap.convertFromImage(*l_pQImage);
        l_QPixmap = l_QPixmap.scaledToWidth(ui->ui_lbl_dataFusion_smallWin_Depth->width());
        ui->ui_lbl_dataFusion_smallWin_Depth->setPixmap(l_QPixmap);
        // Update flag
        m_bDepthImageLoaded = true;
        if (m_bDepthImageLoaded && m_bRGBImageLoaded)
            ui->ui_pb_ApplyFilter->setEnabled(true);
    }
}

void MainWindow::on_ui_pb_LoadGuidImage_clicked()
{
    QFileInfo l_qFileInfo(QFileDialog::getOpenFileName(this, tr("Load guidance image"), "./", tr("Guidance image Files (*.bmp *.png *.pgm)")));
    if (!l_qFileInfo.fileName().isEmpty()) // Go ahead
    {
        QString l_QStrFilename = l_qFileInfo.absoluteFilePath();
        *m_cvRGBImage = cv::imread(l_QStrFilename.toUtf8().constData());
        // Check image dimensions
        if ((m_cvRGBImage->cols != m_sXRes) || (m_cvRGBImage->rows != m_sYRes))
        {
            QMessageBox::information(0, QString("Warning"), QString(QString("Image resolution must be (")+
                                                                    QString::number(m_sXRes)+QString("x")+
                                                                    QString::number(m_sYRes)+QString(")\nImage has been resized!")));
            cv::resize(*m_cvRGBImage, *m_cvRGBImage, cv::Size(m_sXRes, m_sYRes), 0, 0, cv::INTER_NEAREST);
        }
        // Display guidance image
        QImage *l_pQImage;
        QPixmap l_QPixmap;

        if (m_cvRGBImage->channels() == 1)
            cv::cvtColor(*m_cvRGBImage, *m_cvRGBImage, CV_GRAY2RGB);
        else
            cv::cvtColor(*m_cvRGBImage, *m_cvRGBImage, CV_BGR2RGB);
        l_pQImage = new QImage((uchar*)m_cvRGBImage->data, m_cvRGBImage->cols, m_cvRGBImage->rows, QImage::Format_RGB888);
        l_QPixmap.convertFromImage(*l_pQImage);
        l_QPixmap = l_QPixmap.scaledToWidth(ui->ui_lbl_dataFusion_smallWin_RGB->width());
        ui->ui_lbl_dataFusion_smallWin_RGB->setPixmap(l_QPixmap);
        // Update flag
        m_bRGBImageLoaded = true;
        ui->ui_pb_LoadDepthMap->setEnabled(true);
        ui->ui_lbl_dataFusion_smallWin_Depth->setEnabled(true);
        ui->actionLoad_Depth_Map_D->setEnabled(true);
    }
}

void MainWindow::on_ui_pb_ApplyFilter_clicked()
{
    switch (ui->ui_cb_dataFusion_Filter->currentIndex())
    {
        case 0: // RGB-D Filter
            m_DataFusion->DataProcessing((unsigned short*)m_cvDepthImage->data, (unsigned char*)m_cvRGBImage->data, RGBD_FILTER);
            m_DataFusion->UML_RGB_Filter();
        break;
        case 1: // UML Filter (RGB)
            m_DataFusion->DataProcessing((unsigned short*)m_cvDepthImage->data, (unsigned char*)m_cvRGBImage->data, UML_RGB_FILTER);
            m_DataFusion->UML_RGB_Filter();
            break;
        case 2: // UML Filter (Grayscale)
            m_DataFusion->DataProcessing((unsigned short*)m_cvDepthImage->data, (unsigned char*)m_cvRGBImage->data, UML_FILTER);
            m_DataFusion->UML_Filter();
            break;
        case 3: // PWAS Filter (RGB)
            m_DataFusion->DataProcessing((unsigned short*)m_cvDepthImage->data, (unsigned char*)m_cvRGBImage->data, PWAS_RGB_FILTER);
            m_DataFusion->PWAS_RGB_Filter();
            break;
        case 4: // PWAS Filter (Grayscale)
            m_DataFusion->DataProcessing((unsigned short*)m_cvDepthImage->data, (unsigned char*)m_cvRGBImage->data, PWAS_FILTER);
            m_DataFusion->PWAS_Filter();
            break;
        case 5: // JBU Filter
            m_DataFusion->DataProcessing((unsigned short*)m_cvDepthImage->data, (unsigned char*)m_cvRGBImage->data, JBU_FILTER);
            m_DataFusion->JBU_Filter();
            break;
        case 6: // Bilateral Filter
            m_DataFusion->DataProcessing((unsigned short*)m_cvDepthImage->data, NULL, BF_FILTER);
            m_DataFusion->BF_Filter();
            break;
        case 8: // JBU Filter (Kopf et al.)
            m_DataFusion->DataProcessing((unsigned short*)m_cvDepthImage->data, (unsigned char*)m_cvRGBImage->data, LITERATURE_FILTERS);
            m_DataFusion->JBU_Filter_Kopf();
            break;
        case 9: // JBU Filter (Kopf et al.)
            m_DataFusion->DataProcessing((unsigned short*)m_cvDepthImage->data, (unsigned char*)m_cvRGBImage->data, LITERATURE_FILTERS);
            m_DataFusion->NJBU_Filter_Kim();
            break;
        case 10: // NAFDU Filter (Chan et al.)
            m_DataFusion->DataProcessing((unsigned short*)m_cvDepthImage->data, (unsigned char*)m_cvRGBImage->data, LITERATURE_FILTERS);
            m_DataFusion->NAFDU_Filter_Chan();
            break;
    }
    // Refresh output
    ui->ui_cb_dataFusion_Display->setEnabled(true);
    on_ui_cb_dataFusion_Display_currentIndexChanged(ui->ui_cb_dataFusion_Display->currentIndex());
    ui->ui_lbl_maxPalette->setEnabled(true);
    ui->ui_sl_maxPalette->setEnabled(true);
    ui->ui_lbl_minPalette->setEnabled(true);
    ui->ui_sl_minPalette->setEnabled(true);
    ui->actionSave_Displayed_Image->setEnabled(true);
    ui->ui_pb_SaveDisplayedImage->setEnabled(true);
}

void MainWindow::on_ui_cb_dataFusion_Display_currentIndexChanged(int index)
{
    cv::Mat l_cvImage;
    unsigned char* l_ucData;
    unsigned short* l_usData;
    unsigned char* l_ucData_cv;
    unsigned short* l_usData_cv;
    long l_lSize;

    switch (index)
    {
    case 0: // Guidance Image I_RGB
        l_cvImage.create(cv::Size(m_sXRes, m_sYRes), CV_8UC3);
        l_ucData = m_DataFusion->GetGuidanceImageRGB();
        l_ucData_cv = (unsigned char*)l_cvImage.data;
        l_lSize = m_sXRes*m_sYRes*3;
        while (l_lSize--)
            *l_ucData_cv++ = *l_ucData++;
        break;
    case 1: // Guidance Image I
        l_cvImage.create(cv::Size(m_sXRes, m_sYRes), CV_8UC1);
        l_ucData = m_DataFusion->GetGuidanceImage();
        l_ucData_cv = (unsigned char*)l_cvImage.data;
        l_lSize = m_sXRes*m_sYRes;
        while (l_lSize--)
            *l_ucData_cv++ = *l_ucData++;
        cv::cvtColor(l_cvImage, l_cvImage, CV_GRAY2RGB);
        break;
    case 2: // Downsampled Guidance Image I_ds
        l_cvImage.create(cv::Size(m_sXRes>>m_sScaleFactor, m_sYRes>>m_sScaleFactor), CV_8UC1);
        l_ucData = m_DataFusion->GetGuidanceImage_ds();
        l_ucData_cv = (unsigned char*)l_cvImage.data;
        l_lSize = (m_sXRes>>m_sScaleFactor)*(m_sYRes>>m_sScaleFactor);
        while (l_lSize--)
            *l_ucData_cv++ = *l_ucData++;
        cv::cvtColor(l_cvImage, l_cvImage, CV_GRAY2RGB);
        break;
    case 3: // Depth Map D
        l_cvImage.create(cv::Size(m_sXRes, m_sYRes), CV_16UC1);
        l_usData = m_DataFusion->GetDepthMap();
        l_usData_cv = (unsigned short*)l_cvImage.data;
        l_lSize = m_sXRes*m_sYRes;
        while (l_lSize--)
        {
            if (*l_usData == MLF_APP_MAX_DISTANCE)
                *l_usData_cv++ = 0;
            else
                *l_usData_cv++ = *l_usData;
            l_usData++;
        }
        // Note: the followind 16-bit to 8-bit conversion not only chop off the high order byte instead but it scales the values to an 8-bit range
        l_cvImage.convertTo(l_cvImage, CV_8U, 255.0/(m_sMaxPaletteValue-m_sMinPaletteValue), -m_sMinPaletteValue*255.0/(m_sMaxPaletteValue-m_sMinPaletteValue));
        //l_cvImage.convertTo(l_cvImage, CV_8U, 255.0/3500, 0); // Dist. clipping at 3.5m
        cv::cvtColor(l_cvImage, l_cvImage, CV_GRAY2RGB);
        cv::applyColorMap(l_cvImage, l_cvImage, cv::COLORMAP_JET);
        break;
    case 4: // Downsampled Depth Map D_ds
        l_cvImage.create(cv::Size(m_sXRes>>m_sScaleFactor, m_sYRes>>m_sScaleFactor), CV_16UC1);
        l_usData = m_DataFusion->GetDepthMap_ds();
        l_usData_cv = (unsigned short*)l_cvImage.data;
        l_lSize = (m_sXRes>>m_sScaleFactor)*(m_sYRes>>m_sScaleFactor);
        while (l_lSize--)
        {
            if (*l_usData == MLF_APP_MAX_DISTANCE)
                *l_usData_cv++ = 0;
            else
                *l_usData_cv++ = *l_usData;
            l_usData++;
        }
        // Note: the followind 16-bit to 8-bit conversion not only chop off the high order byte instead but it scales the values to an 8-bit range
        l_cvImage.convertTo(l_cvImage, CV_8U, 255.0/(m_sMaxPaletteValue-m_sMinPaletteValue), -m_sMinPaletteValue*255.0/(m_sMaxPaletteValue-m_sMinPaletteValue));
        //l_cvImage.convertTo(l_cvImage, CV_8U, 255.0/3500, 0); // Dist. clipping at 3.5m
        cv::cvtColor(l_cvImage, l_cvImage, CV_GRAY2RGB);
        cv::applyColorMap(l_cvImage, l_cvImage, cv::COLORMAP_JET);
        break;
    case 5: // Credibility Map Q
        l_cvImage.create(cv::Size(m_sXRes, m_sYRes), CV_8UC1);
        l_ucData = m_DataFusion->GetCredibilityMap();
        l_ucData_cv = (unsigned char*)l_cvImage.data;
        l_lSize = m_sXRes*m_sYRes;
        while (l_lSize--)
            *l_ucData_cv++ = 255-*l_ucData++; // Invert the image for display reasons
        cv::cvtColor(l_cvImage, l_cvImage, CV_GRAY2RGB);
        cv::applyColorMap(l_cvImage, l_cvImage, cv::COLORMAP_JET);
        break;
    case 6: // Credibility Map Q_ds
        l_cvImage.create(cv::Size(m_sXRes>>m_sScaleFactor, m_sYRes>>m_sScaleFactor), CV_8UC1);
        l_ucData = m_DataFusion->GetCredibilityMap_ds();
        l_ucData_cv = (unsigned char*)l_cvImage.data;
        l_lSize = (m_sXRes>>m_sScaleFactor)*(m_sYRes>>m_sScaleFactor);
        while (l_lSize--)
            *l_ucData_cv++ = 255-*l_ucData++; // Invert the image for display reasons
        cv::cvtColor(l_cvImage, l_cvImage, CV_GRAY2RGB);
        cv::applyColorMap(l_cvImage, l_cvImage, cv::COLORMAP_JET);
        break;
    case 7: // Blending Mask B
        l_cvImage.create(cv::Size(m_sXRes, m_sYRes), CV_8UC1);
        l_ucData = m_DataFusion->GetBlendingMask();
        l_ucData_cv = (unsigned char*)l_cvImage.data;
        l_lSize = m_sXRes*m_sYRes;
        while (l_lSize--)
            *l_ucData_cv++ = *l_ucData++;
        cv::cvtColor(l_cvImage, l_cvImage, CV_GRAY2RGB);
        cv::applyColorMap(l_cvImage, l_cvImage, cv::COLORMAP_JET);
        break;
    case 8: // Enhanced Depth Map J
        l_cvImage.create(cv::Size(m_sXRes, m_sYRes), CV_16UC1);
        l_usData = m_DataFusion->GetEnhancedDepthData();
        l_usData_cv = (unsigned short*)l_cvImage.data;
        l_lSize = m_sXRes*m_sYRes;
        while (l_lSize--)
        {
            if (*l_usData == MLF_APP_MAX_DISTANCE)
                *l_usData_cv++ = 0;
            else
                *l_usData_cv++ = *l_usData;
            l_usData++;
        }
        // Note: the followind 16-bit to 8-bit conversion not only chop off the high order byte instead but it scales the values to an 8-bit range
        l_cvImage.convertTo(l_cvImage, CV_8U, 255.0/(m_sMaxPaletteValue-m_sMinPaletteValue), -m_sMinPaletteValue*255.0/(m_sMaxPaletteValue-m_sMinPaletteValue));
        cv::cvtColor(l_cvImage, l_cvImage, CV_GRAY2RGB);
        cv::applyColorMap(l_cvImage, l_cvImage, cv::COLORMAP_JET);
        break;
    }
    QImage *l_pqImage;
    l_pqImage = new QImage((uchar*)l_cvImage.data,l_cvImage.cols,l_cvImage.rows,QImage::Format_RGB888);
    QPixmap l_qPixmap;
    l_qPixmap.convertFromImage(*l_pqImage);
    l_qPixmap = l_qPixmap.scaledToWidth(ui->ui_lbl_dataFusion_mainWin->width());
    ui->ui_lbl_dataFusion_mainWin->setPixmap(l_qPixmap);

    // Update depth map thumb
    l_cvImage.create(cv::Size(m_sXRes, m_sYRes), CV_16UC1);
    l_usData = m_DataFusion->GetDepthMap();
    l_usData_cv = (unsigned short*)l_cvImage.data;
    l_lSize = m_sXRes*m_sYRes;
    while (l_lSize--)
    {
        if (*l_usData == MLF_APP_MAX_DISTANCE)
            *l_usData_cv++ = 0;
        else
            *l_usData_cv++ = *l_usData;
        l_usData++;
    }
    // Note: the followind 16-bit to 8-bit conversion not only chop off the high order byte instead but it scales the values to an 8-bit range
    l_cvImage.convertTo(l_cvImage, CV_8U, 255.0/(m_sMaxPaletteValue-m_sMinPaletteValue), -m_sMinPaletteValue*255.0/(m_sMaxPaletteValue-m_sMinPaletteValue));
    //l_cvImage.convertTo(l_cvImage, CV_8U, 255.0/3500, 0); // Dist. clipping at 3.5m
    cv::cvtColor(l_cvImage, l_cvImage, CV_GRAY2RGB);
    cv::applyColorMap(l_cvImage, l_cvImage, cv::COLORMAP_JET);
    l_pqImage = new QImage((uchar*)l_cvImage.data,l_cvImage.cols,l_cvImage.rows,QImage::Format_RGB888);
    l_qPixmap.convertFromImage(*l_pqImage);
    l_qPixmap = l_qPixmap.scaledToWidth(ui->ui_lbl_dataFusion_smallWin_Depth->width());
    ui->ui_lbl_dataFusion_smallWin_Depth->setPixmap(l_qPixmap);
}

void MainWindow::on_ui_pb_SaveDisplayedImage_clicked()
{
    cv::Mat l_cvImage;
    unsigned char* l_ucData;
    unsigned short* l_usData;
    unsigned char* l_ucData_cv;
    unsigned short* l_usData_cv;
    long l_lSize;

    QFileInfo l_qFileInfo(QFileDialog::getSaveFileName(this, tr("Save displayed image"),"./",tr("RAW Image File (*.pgm)")));

    if (!l_qFileInfo.fileName().isEmpty()) // Go ahead
    {
        switch (ui->ui_cb_dataFusion_Display->currentIndex())
        {
        case 0: // Guidance Image I_RGB
            l_cvImage.create(cv::Size(m_sXRes, m_sYRes), CV_8UC3);
            l_ucData = m_DataFusion->GetGuidanceImageRGB();
            l_ucData_cv = (unsigned char*)l_cvImage.data;
            l_lSize = m_sXRes*m_sYRes*3;
            while (l_lSize--)
                *l_ucData_cv++ = *l_ucData++;
            break;
        case 1: // Guidance Image I
            l_cvImage.create(cv::Size(m_sXRes, m_sYRes), CV_8UC1);
            l_ucData = m_DataFusion->GetGuidanceImage();
            l_ucData_cv = (unsigned char*)l_cvImage.data;
            l_lSize = m_sXRes*m_sYRes;
            while (l_lSize--)
                *l_ucData_cv++ = *l_ucData++;
            break;
        case 2: // Downsampled Guidance Image I_ds
            l_cvImage.create(cv::Size(m_sXRes>>m_sScaleFactor, m_sYRes>>m_sScaleFactor), CV_8UC1);
            l_ucData = m_DataFusion->GetGuidanceImage_ds();
            l_ucData_cv = (unsigned char*)l_cvImage.data;
            l_lSize = (m_sXRes>>m_sScaleFactor)*(m_sYRes>>m_sScaleFactor);
            while (l_lSize--)
                *l_ucData_cv++ = *l_ucData++;
            break;
        case 3: // Depth Map D
            l_cvImage.create(cv::Size(m_sXRes, m_sYRes), CV_16UC1);
            l_usData = m_DataFusion->GetDepthMap();
            l_usData_cv = (unsigned short*)l_cvImage.data;
            l_lSize = m_sXRes*m_sYRes;
            while (l_lSize--)
            {
                if (*l_usData == MLF_APP_MAX_DISTANCE)
                    *l_usData_cv++ = 0;
                else
                    *l_usData_cv++ = *l_usData;
                l_usData++;
            }
            break;
        case 4: // Downsampled Depth Map D_ds
            l_cvImage.create(cv::Size(m_sXRes>>m_sScaleFactor, m_sYRes>>m_sScaleFactor), CV_16UC1);
            l_usData = m_DataFusion->GetDepthMap_ds();
            l_usData_cv = (unsigned short*)l_cvImage.data;
            l_lSize = (m_sXRes>>m_sScaleFactor)*(m_sYRes>>m_sScaleFactor);
            while (l_lSize--)
            {
                if (*l_usData == MLF_APP_MAX_DISTANCE)
                    *l_usData_cv++ = 0;
                else
                    *l_usData_cv++ = *l_usData;
                l_usData++;
            }
            break;
        case 5: // Credibility Map Q
            l_cvImage.create(cv::Size(m_sXRes, m_sYRes), CV_8UC1);
            l_ucData = m_DataFusion->GetCredibilityMap();
            l_ucData_cv = (unsigned char*)l_cvImage.data;
            l_lSize = m_sXRes*m_sYRes;
            while (l_lSize--)
                *l_ucData_cv++ = 255-*l_ucData++; // Invert the image for display reasons
            break;
        case 6: // Credibility Map Q_ds
            l_cvImage.create(cv::Size(m_sXRes>>m_sScaleFactor, m_sYRes>>m_sScaleFactor), CV_8UC1);
            l_ucData = m_DataFusion->GetCredibilityMap_ds();
            l_ucData_cv = (unsigned char*)l_cvImage.data;
            l_lSize = (m_sXRes>>m_sScaleFactor)*(m_sYRes>>m_sScaleFactor);
            while (l_lSize--)
                *l_ucData_cv++ = 255-*l_ucData++; // Invert the image for display reasons
            break;
        case 7: // Blending Mask B
            l_cvImage.create(cv::Size(m_sXRes, m_sYRes), CV_8UC1);
            l_ucData = m_DataFusion->GetBlendingMask();
            l_ucData_cv = (unsigned char*)l_cvImage.data;
            l_lSize = m_sXRes*m_sYRes;
            while (l_lSize--)
                *l_ucData_cv++ = *l_ucData++;
            break;
        case 8: // Enhanced Depth Map J
            l_cvImage.create(cv::Size(m_sXRes, m_sYRes), CV_16UC1);
            l_usData = m_DataFusion->GetEnhancedDepthData();
            l_usData_cv = (unsigned short*)l_cvImage.data;
            l_lSize = m_sXRes*m_sYRes;
            while (l_lSize--)
            {
                if (*l_usData == MLF_APP_MAX_DISTANCE)
                    *l_usData_cv++ = 0;
                else
                    *l_usData_cv++ = *l_usData;
                l_usData++;
            }
            break;
        }
        // Save image
        cv::imwrite(QString(l_qFileInfo.absolutePath()+QString("/")+l_qFileInfo.baseName()+".pgm").toUtf8().constData(), l_cvImage);
    }
}

void MainWindow::on_actionAbout_triggered()
{
    QMessageBox::information(0, QString("About"), QString(
                                QString("This tool implements the RGB-D, UML and PWAS filters presented at:\n\n")+
                                QString("Real-Time Depth Enhancement by Fusion for RGB-D Cameras\n")+
                                QString("Garcia, F. and Aouada, D. and Solignac, T. and Mirbach, B. and Ottersten, B.\n")+
                                QString("IET Computer Vision, v.7:5, pp.1-11, 2013\n\n")+
                                QString("A New Multi-lateral Filter for Real-Time Depth Enhancement\n")+
                                QString("Garcia, F. and Aouada, D. and Mirbach, B. and Solignac, T. and Ottersten, B.\n")+
                                QString("IEEE International Conference on Advanced Video and Signal-Based Surveillance (AVSS), pp.42-47, 2011\n\n")+
                                QString("Pixel Weighted Average Strategy for Depth Sensor Data Fusion\n")+
                                QString("Garcia, F. and Mirbach, B. and Ottersten, B. and Grandidier, F. and Cuesta, A.\n")+
                                QString("IEEE International Conference on Image Processing (ICIP), pp.2805-2808, 2010\n\n")+
                                QString("Download the papers at http://www.frederic-garcia-becerro.com")));
}

void MainWindow::on_ui_cb_dataFusion_Downsampling_currentIndexChanged(int index)
{
    m_sScaleFactor = index;
    delete m_DataFusion;
    m_DataFusion = NULL;
    m_DataFusion = new c_DataFusion(D_Fusion_2D_3D_2DW, D_Fusion_2D_3D_2DH, m_sScaleFactor);
}
