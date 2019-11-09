QT += core
QT -= gui

CONFIG += c++11

TARGET = difequation
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    matrix_3diagonal.cpp

HEADERS += \
    matrix_3diagonal.h \
    functions.h
