# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'mainwindow.ui'
#
# Created: Fri May 11 14:10:03 2012
#      by: PyQt4 UI code generator 4.8.5
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(776, 635)
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "ColorComposer", None, QtGui.QApplication.UnicodeUTF8))
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.gridLayout = QtGui.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setSizeConstraint(QtGui.QLayout.SetMaximumSize)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.label = QtGui.QLabel(self.centralwidget)
        self.label.setText(QtGui.QApplication.translate("MainWindow", "Preview", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setObjectName(_fromUtf8("label"))
        self.verticalLayout_2.addWidget(self.label)
        self.graphicsView = QtGui.QGraphicsView(self.centralwidget)
        self.graphicsView.setObjectName(_fromUtf8("graphicsView"))
        self.verticalLayout_2.addWidget(self.graphicsView)
        self.gridLayout.addLayout(self.verticalLayout_2, 0, 0, 1, 1)
        self.verticalLayout_3 = QtGui.QVBoxLayout()
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.fileChooseLabel = QtGui.QLabel(self.centralwidget)
        self.fileChooseLabel.setText(QtGui.QApplication.translate("MainWindow", "File List", None, QtGui.QApplication.UnicodeUTF8))
        self.fileChooseLabel.setObjectName(_fromUtf8("fileChooseLabel"))
        self.verticalLayout_3.addWidget(self.fileChooseLabel)
        self.fileWidget = QtGui.QListWidget(self.centralwidget)
        self.fileWidget.setObjectName(_fromUtf8("fileWidget"))
        self.verticalLayout_3.addWidget(self.fileWidget)
        self.addFileButton = QtGui.QPushButton(self.centralwidget)
        self.addFileButton.setText(QtGui.QApplication.translate("MainWindow", "Add File...", None, QtGui.QApplication.UnicodeUTF8))
        self.addFileButton.setObjectName(_fromUtf8("addFileButton"))
        self.verticalLayout_3.addWidget(self.addFileButton)
        self.doTransformationButton = QtGui.QPushButton(self.centralwidget)
        self.doTransformationButton.setText(QtGui.QApplication.translate("MainWindow", "Transformation", None, QtGui.QApplication.UnicodeUTF8))
        self.doTransformationButton.setObjectName(_fromUtf8("doTransformationButton"))
        self.verticalLayout_3.addWidget(self.doTransformationButton)
        spacerItem = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_3.addItem(spacerItem)
        self.deleteBeadButton = QtGui.QPushButton(self.centralwidget)
        self.deleteBeadButton.setText(QtGui.QApplication.translate("MainWindow", "Delete Bead", None, QtGui.QApplication.UnicodeUTF8))
        self.deleteBeadButton.setObjectName(_fromUtf8("deleteBeadButton"))
        self.verticalLayout_3.addWidget(self.deleteBeadButton)
        self.horizontalLayout_4 = QtGui.QHBoxLayout()
        self.horizontalLayout_4.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout_4.setObjectName(_fromUtf8("horizontalLayout_4"))
        self.addRedBeadButton = QtGui.QPushButton(self.centralwidget)
        self.addRedBeadButton.setText(QtGui.QApplication.translate("MainWindow", "Add Red Bead", None, QtGui.QApplication.UnicodeUTF8))
        self.addRedBeadButton.setObjectName(_fromUtf8("addRedBeadButton"))
        self.horizontalLayout_4.addWidget(self.addRedBeadButton)
        spacerItem1 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_4.addItem(spacerItem1)
        self.addGreenBeadButton = QtGui.QPushButton(self.centralwidget)
        self.addGreenBeadButton.setText(QtGui.QApplication.translate("MainWindow", "Add Green Bead", None, QtGui.QApplication.UnicodeUTF8))
        self.addGreenBeadButton.setObjectName(_fromUtf8("addGreenBeadButton"))
        self.horizontalLayout_4.addWidget(self.addGreenBeadButton)
        self.verticalLayout_3.addLayout(self.horizontalLayout_4)
        spacerItem2 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_3.addItem(spacerItem2)
        self.label_3 = QtGui.QLabel(self.centralwidget)
        self.label_3.setText(QtGui.QApplication.translate("MainWindow", "Percent of occurrency", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.verticalLayout_3.addWidget(self.label_3)
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.PercentOfOccurrencySlider = QtGui.QSlider(self.centralwidget)
        self.PercentOfOccurrencySlider.setMaximum(1000)
        self.PercentOfOccurrencySlider.setPageStep(5)
        self.PercentOfOccurrencySlider.setProperty("value", 950)
        self.PercentOfOccurrencySlider.setSliderPosition(950)
        self.PercentOfOccurrencySlider.setOrientation(QtCore.Qt.Horizontal)
        self.PercentOfOccurrencySlider.setObjectName(_fromUtf8("PercentOfOccurrencySlider"))
        self.horizontalLayout_3.addWidget(self.PercentOfOccurrencySlider)
        self.PercentOfOccurrencySpinBox = QtGui.QDoubleSpinBox(self.centralwidget)
        self.PercentOfOccurrencySpinBox.setDecimals(1)
        self.PercentOfOccurrencySpinBox.setMaximum(100.0)
        self.PercentOfOccurrencySpinBox.setSingleStep(1.0)
        self.PercentOfOccurrencySpinBox.setProperty("value", 95.0)
        self.PercentOfOccurrencySpinBox.setObjectName(_fromUtf8("PercentOfOccurrencySpinBox"))
        self.horizontalLayout_3.addWidget(self.PercentOfOccurrencySpinBox)
        self.verticalLayout_3.addLayout(self.horizontalLayout_3)
        self.label_2 = QtGui.QLabel(self.centralwidget)
        self.label_2.setText(QtGui.QApplication.translate("MainWindow", "Cursor Radius:", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.verticalLayout_3.addWidget(self.label_2)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.cursorRadiusSlider = QtGui.QSlider(self.centralwidget)
        self.cursorRadiusSlider.setMaximum(1000)
        self.cursorRadiusSlider.setProperty("value", 100)
        self.cursorRadiusSlider.setOrientation(QtCore.Qt.Horizontal)
        self.cursorRadiusSlider.setObjectName(_fromUtf8("cursorRadiusSlider"))
        self.horizontalLayout.addWidget(self.cursorRadiusSlider)
        self.cursorRadiusSpinBox = QtGui.QDoubleSpinBox(self.centralwidget)
        self.cursorRadiusSpinBox.setDecimals(1)
        self.cursorRadiusSpinBox.setMaximum(1000.0)
        self.cursorRadiusSpinBox.setProperty("value", 10.0)
        self.cursorRadiusSpinBox.setObjectName(_fromUtf8("cursorRadiusSpinBox"))
        self.horizontalLayout.addWidget(self.cursorRadiusSpinBox)
        self.verticalLayout_3.addLayout(self.horizontalLayout)
        self.label_5 = QtGui.QLabel(self.centralwidget)
        self.label_5.setText(QtGui.QApplication.translate("MainWindow", "AutoDetection:", None, QtGui.QApplication.UnicodeUTF8))
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.verticalLayout_3.addWidget(self.label_5)
        self.label_4 = QtGui.QLabel(self.centralwidget)
        self.label_4.setText(QtGui.QApplication.translate("MainWindow", "Maximal Variance                 Maximal Distance", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.verticalLayout_3.addWidget(self.label_4)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.maximalVarianceSlider = QtGui.QSlider(self.centralwidget)
        self.maximalVarianceSlider.setMaximum(50)
        self.maximalVarianceSlider.setProperty("value", 10)
        self.maximalVarianceSlider.setOrientation(QtCore.Qt.Horizontal)
        self.maximalVarianceSlider.setObjectName(_fromUtf8("maximalVarianceSlider"))
        self.horizontalLayout_2.addWidget(self.maximalVarianceSlider)
        self.maximalVarianceSpinBox = QtGui.QDoubleSpinBox(self.centralwidget)
        self.maximalVarianceSpinBox.setDecimals(1)
        self.maximalVarianceSpinBox.setMinimum(0.0)
        self.maximalVarianceSpinBox.setMaximum(5.0)
        self.maximalVarianceSpinBox.setSingleStep(0.5)
        self.maximalVarianceSpinBox.setProperty("value", 1.0)
        self.maximalVarianceSpinBox.setObjectName(_fromUtf8("maximalVarianceSpinBox"))
        self.horizontalLayout_2.addWidget(self.maximalVarianceSpinBox)
        self.maximalDistanceSlider = QtGui.QSlider(self.centralwidget)
        self.maximalDistanceSlider.setMaximum(100)
        self.maximalDistanceSlider.setProperty("value", 30)
        self.maximalDistanceSlider.setOrientation(QtCore.Qt.Horizontal)
        self.maximalDistanceSlider.setObjectName(_fromUtf8("maximalDistanceSlider"))
        self.horizontalLayout_2.addWidget(self.maximalDistanceSlider)
        self.maximalDistanceSpinBox = QtGui.QSpinBox(self.centralwidget)
        self.maximalDistanceSpinBox.setMaximum(10)
        self.maximalDistanceSpinBox.setProperty("value", 3)
        self.maximalDistanceSpinBox.setObjectName(_fromUtf8("maximalDistanceSpinBox"))
        self.horizontalLayout_2.addWidget(self.maximalDistanceSpinBox)
        self.verticalLayout_3.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_5 = QtGui.QHBoxLayout()
        self.horizontalLayout_5.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout_5.setObjectName(_fromUtf8("horizontalLayout_5"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setContentsMargins(-1, 0, -1, -1)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.statsNumber = QtGui.QLabel(self.centralwidget)
        self.statsNumber.setText(QtGui.QApplication.translate("MainWindow", "Statistics: no beads", None, QtGui.QApplication.UnicodeUTF8))
        self.statsNumber.setObjectName(_fromUtf8("statsNumber"))
        self.verticalLayout.addWidget(self.statsNumber)
        self.statsIntensity = QtGui.QLabel(self.centralwidget)
        self.statsIntensity.setText(_fromUtf8(""))
        self.statsIntensity.setObjectName(_fromUtf8("statsIntensity"))
        self.verticalLayout.addWidget(self.statsIntensity)
        self.horizontalLayout_5.addLayout(self.verticalLayout)
        self.verticalLayout_4 = QtGui.QVBoxLayout()
        self.verticalLayout_4.setSpacing(6)
        self.verticalLayout_4.setSizeConstraint(QtGui.QLayout.SetDefaultConstraint)
        self.verticalLayout_4.setObjectName(_fromUtf8("verticalLayout_4"))
        self.label_6 = QtGui.QLabel(self.centralwidget)
        self.label_6.setText(QtGui.QApplication.translate("MainWindow", "Frames:", None, QtGui.QApplication.UnicodeUTF8))
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.verticalLayout_4.addWidget(self.label_6)
        self.statsNumberOfFrames = QtGui.QLabel(self.centralwidget)
        self.statsNumberOfFrames.setText(_fromUtf8(""))
        self.statsNumberOfFrames.setObjectName(_fromUtf8("statsNumberOfFrames"))
        self.verticalLayout_4.addWidget(self.statsNumberOfFrames)
        self.horizontalLayout_5.addLayout(self.verticalLayout_4)
        self.verticalLayout_3.addLayout(self.horizontalLayout_5)
        spacerItem3 = QtGui.QSpacerItem(253, 13, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_3.addItem(spacerItem3)
        self.zoomOutButton = QtGui.QPushButton(self.centralwidget)
        self.zoomOutButton.setText(QtGui.QApplication.translate("MainWindow", "zoom -", None, QtGui.QApplication.UnicodeUTF8))
        self.zoomOutButton.setObjectName(_fromUtf8("zoomOutButton"))
        self.verticalLayout_3.addWidget(self.zoomOutButton)
        self.zoomInButton = QtGui.QPushButton(self.centralwidget)
        self.zoomInButton.setText(QtGui.QApplication.translate("MainWindow", "zoom +", None, QtGui.QApplication.UnicodeUTF8))
        self.zoomInButton.setObjectName(_fromUtf8("zoomInButton"))
        self.verticalLayout_3.addWidget(self.zoomInButton)
        self.gridLayout.addLayout(self.verticalLayout_3, 0, 1, 1, 1)
        self.gridLayout.setColumnStretch(0, 20)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 776, 25))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setTitle(QtGui.QApplication.translate("MainWindow", "File", None, QtGui.QApplication.UnicodeUTF8))
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        self.menuView = QtGui.QMenu(self.menubar)
        self.menuView.setTitle(QtGui.QApplication.translate("MainWindow", "View", None, QtGui.QApplication.UnicodeUTF8))
        self.menuView.setObjectName(_fromUtf8("menuView"))
        self.menuExtras = QtGui.QMenu(self.menubar)
        self.menuExtras.setTitle(QtGui.QApplication.translate("MainWindow", "Extras", None, QtGui.QApplication.UnicodeUTF8))
        self.menuExtras.setObjectName(_fromUtf8("menuExtras"))
        self.menuHelp = QtGui.QMenu(self.menubar)
        self.menuHelp.setTitle(QtGui.QApplication.translate("MainWindow", "Help", None, QtGui.QApplication.UnicodeUTF8))
        self.menuHelp.setObjectName(_fromUtf8("menuHelp"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)
        self.actionTransformation_Matrix = QtGui.QAction(MainWindow)
        self.actionTransformation_Matrix.setText(QtGui.QApplication.translate("MainWindow", "Transformation Matrix...", None, QtGui.QApplication.UnicodeUTF8))
        self.actionTransformation_Matrix.setObjectName(_fromUtf8("actionTransformation_Matrix"))
        self.actionExport_Composed_image = QtGui.QAction(MainWindow)
        self.actionExport_Composed_image.setText(QtGui.QApplication.translate("MainWindow", "Export Composed image...", None, QtGui.QApplication.UnicodeUTF8))
        self.actionExport_Composed_image.setObjectName(_fromUtf8("actionExport_Composed_image"))
        self.actionExport_as_stack_of_images = QtGui.QAction(MainWindow)
        self.actionExport_as_stack_of_images.setText(QtGui.QApplication.translate("MainWindow", "Export as stack of images...", None, QtGui.QApplication.UnicodeUTF8))
        self.actionExport_as_stack_of_images.setObjectName(_fromUtf8("actionExport_as_stack_of_images"))
        self.actionAuto_detect_beads = QtGui.QAction(MainWindow)
        self.actionAuto_detect_beads.setText(QtGui.QApplication.translate("MainWindow", "Auto-detect beads", None, QtGui.QApplication.UnicodeUTF8))
        self.actionAuto_detect_beads.setObjectName(_fromUtf8("actionAuto_detect_beads"))
        self.actionAbout = QtGui.QAction(MainWindow)
        self.actionAbout.setText(QtGui.QApplication.translate("MainWindow", "About", None, QtGui.QApplication.UnicodeUTF8))
        self.actionAbout.setObjectName(_fromUtf8("actionAbout"))
        self.actionClear_all = QtGui.QAction(MainWindow)
        self.actionClear_all.setText(QtGui.QApplication.translate("MainWindow", "Discard loaded images (Reset)", None, QtGui.QApplication.UnicodeUTF8))
        self.actionClear_all.setObjectName(_fromUtf8("actionClear_all"))
        self.actionExport_transformed_coordinates = QtGui.QAction(MainWindow)
        self.actionExport_transformed_coordinates.setText(QtGui.QApplication.translate("MainWindow", "Export transformed Coordinates as List...", None, QtGui.QApplication.UnicodeUTF8))
        self.actionExport_transformed_coordinates.setObjectName(_fromUtf8("actionExport_transformed_coordinates"))
        self.actionBig_bead_circles = QtGui.QAction(MainWindow)
        self.actionBig_bead_circles.setText(QtGui.QApplication.translate("MainWindow", "Big bead circles", None, QtGui.QApplication.UnicodeUTF8))
        self.actionBig_bead_circles.setObjectName(_fromUtf8("actionBig_bead_circles"))
        self.menuFile.addAction(self.actionExport_Composed_image)
        self.menuFile.addAction(self.actionExport_as_stack_of_images)
        self.menuFile.addAction(self.actionExport_transformed_coordinates)
        self.menuView.addAction(self.actionTransformation_Matrix)
        self.menuExtras.addAction(self.actionAuto_detect_beads)
        self.menuExtras.addAction(self.actionClear_all)
        self.menuHelp.addAction(self.actionAbout)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuView.menuAction())
        self.menubar.addAction(self.menuExtras.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        pass

