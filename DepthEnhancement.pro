#-------------------------------------------------
#
# Project created by QtCreator 2015-03-17T15:48:19
# Author: Frederic GARCIA BECERRO
# Email: frederic.garcia.becerro@gmail.com
# Website: http://www.frederic-garcia-becerro.com
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = DepthEnhancement
TEMPLATE = app

CONFIG += console

SOURCES += main.cpp\
        mainwindow.cpp \
        c_datafusion.cpp


HEADERS  += mainwindow.h \
         app_types.h \
         c_datafusion.h

FORMS    += mainwindow.ui

win32:CONFIG(release, debug|release): DEFINES += QT_NO_DEBUG_OUTPUT
win32:CONFIG(debug, debug|release): DEFINES += QT_DEBUG_MODE

#add OpenCV Libraries
INCLUDEPATH += C:/OpenCv/OpenCv2_4_4/include
INCLUDEPATH += C:/OpenCv/OpenCv2_4_4/include/opencv
win32:CONFIG(release, debug|release): LIBS += -LC:/OpenCv/OpenCv2_4_4/Lib \
        -lopencv_core244 \
        -lopencv_calib3d244 \
        -lopencv_contrib244 \
        -lopencv_features2d244 \
        -lopencv_flann244 \
        -lopencv_gpu244 \
        -lopencv_highgui244 \
        -lopencv_imgproc244 \
        -lopencv_legacy244 \
        -lopencv_ml244 \
        -lopencv_nonfree244 \
        -lopencv_objdetect244 \
        -lopencv_photo244 \
        -lopencv_stitching244 \
        -lopencv_video244 \
        -lopencv_videostab244
else: LIBS += -LC:/OpenCv/OpenCv2_4_4/Lib \
        -lopencv_core244d \
        -lopencv_calib3d244d \
        -lopencv_contrib244d \
        -lopencv_features2d244d \
        -lopencv_flann244d \
        -lopencv_gpu244d \
        -lopencv_highgui244d \
        -lopencv_imgproc244d \
        -lopencv_legacy244d \
        -lopencv_ml244d \
        -lopencv_nonfree244d \
        -lopencv_objdetect244d \
        -lopencv_photo244d \
        -lopencv_stitching244d \
        -lopencv_video244d \
        -lopencv_videostab244d

#add destination directory
win32:CONFIG(release, debug|release): DESTDIR = D:\\Documents\\Research\\src_code\\DepthEnhancement\\bin\\r
else: DESTDIR = D:\\Documents\\Research\\src_code\\DepthEnhancement\\bin\\d

win32:CONFIG(release, debug|release): QMAKE_POST_LINK  += copy C:\\Qt\\Qt5.3.2\\5.3\\msvc2010_opengl\\bin\\Qt5Core.dll $${DESTDIR}\\Qt5Core.dll \
                       & copy C:\\Qt\\Qt5.3.2\\5.3\\msvc2010_opengl\\bin\\Qt5Gui.dll $${DESTDIR}\\Qt5Gui.dll \
                       & copy C:\\Qt\\Qt5.3.2\\5.3\\msvc2010_opengl\\bin\\Qt5OpenGL.dll $${DESTDIR}\\Qt5OpenGL.dll \
                       & copy C:\\Qt\\Qt5.3.2\\5.3\\msvc2010_opengl\\bin\\icuin52.dll $${DESTDIR}\\icuin52.dll \
                       & copy C:\\Qt\\Qt5.3.2\\5.3\\msvc2010_opengl\\bin\\icudt52.dll $${DESTDIR}\\icudt52.dll \
                       & copy C:\\Qt\\Qt5.3.2\\5.3\\msvc2010_opengl\\bin\\icuuc52.dll $${DESTDIR}\\icuuc52.dll \
                       & copy C:\\Qt\\Qt5.3.2\\5.3\\msvc2010_opengl\\bin\\Qt5Widgets.dll $${DESTDIR}\\Qt5Widgets.dll \
                       & mkdir $${DESTDIR}\\platforms \
                       & copy C:\\Qt\\Qt5.3.2\\5.3\\msvc2010_opengl\\plugins\\platforms\\qwindows.dll $${DESTDIR}\\platforms\\qwindows.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_core244.dll $${DESTDIR}\\opencv_core244.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_calib3d244.dll $${DESTDIR}\\opencv_calib3d244.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_contrib244.dll $${DESTDIR}\\opencv_contrib244.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_features2d244.dll $${DESTDIR}\\opencv_features2d244.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_flann244.dll $${DESTDIR}\\opencv_flann244.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_gpu244.dll $${DESTDIR}\\opencv_gpu244.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_highgui244.dll $${DESTDIR}\\opencv_highgui244.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_imgproc244.dll $${DESTDIR}\\opencv_imgproc244.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_legacy244.dll $${DESTDIR}\\opencv_legacy244.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_ml244.dll $${DESTDIR}\\opencv_ml244.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_nonfree244.dll $${DESTDIR}\\opencv_nonfree244.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_objdetect244.dll $${DESTDIR}\\opencv_objdetect244.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_photo244.dll $${DESTDIR}\\opencv_photo244.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_stitching244.dll $${DESTDIR}\\opencv_stitching244.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_video244.dll $${DESTDIR}\\opencv_video244.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_videostab244.dll $${DESTDIR}\\opencv_videostab244.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_ffmpeg244.dll $${DESTDIR}\\opencv_ffmpeg244.dll
else: QMAKE_POST_LINK  += copy C:\\Qt\\Qt5.3.2\\5.3\\msvc2010_opengl\\bin\\Qt5Cored.dll $${DESTDIR}\\Qt5Cored.dll \
                       & copy C:\\Qt\\Qt5.3.2\\5.3\\msvc2010_opengl\\bin\\Qt5Guid.dll $${DESTDIR}\\Qt5Guid.dll \
                       & copy C:\\Qt\\Qt5.3.2\\5.3\\msvc2010_opengl\\bin\\Qt5OpenGLd.dll $${DESTDIR}\\Qt5OpenGLd.dll \
                       & copy C:\\Qt\\Qt5.3.2\\5.3\\msvc2010_opengl\\bin\\icuin52.dll $${DESTDIR}\\icuin52.dll \
                       & copy C:\\Qt\\Qt5.3.2\\5.3\\msvc2010_opengl\\bin\\icudt52.dll $${DESTDIR}\\icudt52.dll \
                       & copy C:\\Qt\\Qt5.3.2\\5.3\\msvc2010_opengl\\bin\\icuuc52.dll $${DESTDIR}\\icuuc52.dll \
                       & copy C:\\Qt\\Qt5.3.2\\5.3\\msvc2010_opengl\\bin\\Qt5Widgetsd.dll $${DESTDIR}\\Qt5Widgetsd.dll \
                       & mkdir $${DESTDIR}\\platforms \
                       & copy C:\\Qt\\Qt5.3.2\\5.3\\msvc2010_opengl\\plugins\\platforms\\qwindowsd.dll $${DESTDIR}\\platforms\\qwindowsd.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_core244d.dll $${DESTDIR}\\opencv_core244d.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_calib3d244d.dll $${DESTDIR}\\opencv_calib3d244d.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_contrib244d.dll $${DESTDIR}\\opencv_contrib244d.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_features2d244d.dll $${DESTDIR}\\opencv_features2d244d.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_flann244d.dll $${DESTDIR}\\opencv_flann244d.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_gpu244d.dll $${DESTDIR}\\opencv_gpu244d.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_highgui244d.dll $${DESTDIR}\\opencv_highgui244d.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_imgproc244d.dll $${DESTDIR}\\opencv_imgproc244d.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_legacy244d.dll $${DESTDIR}\\opencv_legacy244d.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_ml244d.dll $${DESTDIR}\\opencv_ml244d.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_nonfree244d.dll $${DESTDIR}\\opencv_nonfree244d.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_objdetect244d.dll $${DESTDIR}\\opencv_objdetect244d.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_photo244d.dll $${DESTDIR}\\opencv_photo244d.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_stitching244d.dll $${DESTDIR}\\opencv_stitching244d.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_video244d.dll $${DESTDIR}\\opencv_video244d.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_videostab244d.dll $${DESTDIR}\\opencv_videostab244d.dll \
                       & copy C:\\OpenCv\\OpenCv2_4_4\\bin\\opencv_ffmpeg244.dll $${DESTDIR}\\opencv_ffmpeg244.dll
