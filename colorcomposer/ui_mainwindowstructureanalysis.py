# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'mainwindowstructureanalysis.ui'
#
# Created: Tue Jul 17 14:42:55 2012
#      by: PyQt4 UI code generator 4.8.5
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_MainWindowStructureAnalysis(object):
    def setupUi(self, MainWindowStructureAnalysis):
        MainWindowStructureAnalysis.setObjectName(_fromUtf8("MainWindowStructureAnalysis"))
        MainWindowStructureAnalysis.resize(900, 600)
        MainWindowStructureAnalysis.setWindowTitle(QtGui.QApplication.translate("MainWindowStructureAnalysis", "MainWindow", None, QtGui.QApplication.UnicodeUTF8))
        self.centralwidget = QtGui.QWidget(MainWindowStructureAnalysis)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.horizontalLayout = QtGui.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.graphicsView = QtGui.QGraphicsView(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(200)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.graphicsView.sizePolicy().hasHeightForWidth())
        self.graphicsView.setSizePolicy(sizePolicy)
        self.graphicsView.setObjectName(_fromUtf8("graphicsView"))
        self.horizontalLayout.addWidget(self.graphicsView)
        self.verticalLayout_4 = QtGui.QVBoxLayout()
        self.verticalLayout_4.setContentsMargins(-1, -1, 0, -1)
        self.verticalLayout_4.setObjectName(_fromUtf8("verticalLayout_4"))
        self.fileWidget = QtGui.QListWidget(self.centralwidget)
        self.fileWidget.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(100)
        sizePolicy.setHeightForWidth(self.fileWidget.sizePolicy().hasHeightForWidth())
        self.fileWidget.setSizePolicy(sizePolicy)
        self.fileWidget.setObjectName(_fromUtf8("fileWidget"))
        self.verticalLayout_4.addWidget(self.fileWidget)
        self.addFileButton = QtGui.QPushButton(self.centralwidget)
        self.addFileButton.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "add File", None, QtGui.QApplication.UnicodeUTF8))
        self.addFileButton.setObjectName(_fromUtf8("addFileButton"))
        self.verticalLayout_4.addWidget(self.addFileButton)
        spacerItem = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_4.addItem(spacerItem)
        self.selectChannelcomboBox = QtGui.QComboBox(self.centralwidget)
        self.selectChannelcomboBox.setObjectName(_fromUtf8("selectChannelcomboBox"))
        self.verticalLayout_4.addWidget(self.selectChannelcomboBox)
        self.verticalLayout_5 = QtGui.QVBoxLayout()
        self.verticalLayout_5.setContentsMargins(-1, 0, -1, -1)
        self.verticalLayout_5.setObjectName(_fromUtf8("verticalLayout_5"))
        self.label = QtGui.QLabel(self.centralwidget)
        self.label.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "set sigma for gaussian smoothing", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setObjectName(_fromUtf8("label"))
        self.verticalLayout_5.addWidget(self.label)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.setSigmaSlider = QtGui.QSlider(self.centralwidget)
        self.setSigmaSlider.setMaximum(100)
        self.setSigmaSlider.setPageStep(5)
        self.setSigmaSlider.setProperty("value", 10)
        self.setSigmaSlider.setOrientation(QtCore.Qt.Horizontal)
        self.setSigmaSlider.setObjectName(_fromUtf8("setSigmaSlider"))
        self.horizontalLayout_2.addWidget(self.setSigmaSlider)
        self.setSigmaSpinBox = QtGui.QDoubleSpinBox(self.centralwidget)
        self.setSigmaSpinBox.setDecimals(1)
        self.setSigmaSpinBox.setMaximum(100.0)
        self.setSigmaSpinBox.setSingleStep(0.1)
        self.setSigmaSpinBox.setProperty("value", 1.0)
        self.setSigmaSpinBox.setObjectName(_fromUtf8("setSigmaSpinBox"))
        self.horizontalLayout_2.addWidget(self.setSigmaSpinBox)
        self.doSmoothingButton = QtGui.QPushButton(self.centralwidget)
        self.doSmoothingButton.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "Do Smoothing", None, QtGui.QApplication.UnicodeUTF8))
        self.doSmoothingButton.setObjectName(_fromUtf8("doSmoothingButton"))
        self.horizontalLayout_2.addWidget(self.doSmoothingButton)
        self.verticalLayout_5.addLayout(self.horizontalLayout_2)
        self.verticalLayout_4.addLayout(self.verticalLayout_5)
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setContentsMargins(-1, 0, -1, -1)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.label_2 = QtGui.QLabel(self.centralwidget)
        self.label_2.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "thresholding parameter", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.verticalLayout.addWidget(self.label_2)
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setContentsMargins(-1, 0, 0, -1)
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.thresholdSlider = QtGui.QSlider(self.centralwidget)
        self.thresholdSlider.setMaximum(1000)
        self.thresholdSlider.setPageStep(5)
        self.thresholdSlider.setProperty("value", 500)
        self.thresholdSlider.setOrientation(QtCore.Qt.Horizontal)
        self.thresholdSlider.setObjectName(_fromUtf8("thresholdSlider"))
        self.horizontalLayout_3.addWidget(self.thresholdSlider)
        self.thresholdSpinBox = QtGui.QDoubleSpinBox(self.centralwidget)
        self.thresholdSpinBox.setDecimals(1)
        self.thresholdSpinBox.setSingleStep(1.0)
        self.thresholdSpinBox.setProperty("value", 50.0)
        self.thresholdSpinBox.setObjectName(_fromUtf8("thresholdSpinBox"))
        self.horizontalLayout_3.addWidget(self.thresholdSpinBox)
        self.thresholdButton = QtGui.QPushButton(self.centralwidget)
        self.thresholdButton.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "Threshold", None, QtGui.QApplication.UnicodeUTF8))
        self.thresholdButton.setObjectName(_fromUtf8("thresholdButton"))
        self.horizontalLayout_3.addWidget(self.thresholdButton)
        self.verticalLayout.addLayout(self.horizontalLayout_3)
        self.verticalLayout_4.addLayout(self.verticalLayout)
        self.findConnectedComponentsButton = QtGui.QPushButton(self.centralwidget)
        self.findConnectedComponentsButton.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "find Connected Components", None, QtGui.QApplication.UnicodeUTF8))
        self.findConnectedComponentsButton.setObjectName(_fromUtf8("findConnectedComponentsButton"))
        self.verticalLayout_4.addWidget(self.findConnectedComponentsButton)
        self.horizontalLayout_5 = QtGui.QHBoxLayout()
        self.horizontalLayout_5.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout_5.setObjectName(_fromUtf8("horizontalLayout_5"))
        self.cutOffButton = QtGui.QPushButton(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(40)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.cutOffButton.sizePolicy().hasHeightForWidth())
        self.cutOffButton.setSizePolicy(sizePolicy)
        self.cutOffButton.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "Cut off CC smaller than", None, QtGui.QApplication.UnicodeUTF8))
        self.cutOffButton.setObjectName(_fromUtf8("cutOffButton"))
        self.horizontalLayout_5.addWidget(self.cutOffButton)
        self.numberCutOffspinBox = QtGui.QSpinBox(self.centralwidget)
        self.numberCutOffspinBox.setMaximum(999)
        self.numberCutOffspinBox.setProperty("value", 10)
        self.numberCutOffspinBox.setObjectName(_fromUtf8("numberCutOffspinBox"))
        self.horizontalLayout_5.addWidget(self.numberCutOffspinBox)
        self.verticalLayout_4.addLayout(self.horizontalLayout_5)
        self.verticalLayout_8 = QtGui.QVBoxLayout()
        self.verticalLayout_8.setSpacing(0)
        self.verticalLayout_8.setContentsMargins(-1, 0, -1, -1)
        self.verticalLayout_8.setObjectName(_fromUtf8("verticalLayout_8"))
        self.horizontalLayout_4 = QtGui.QHBoxLayout()
        self.horizontalLayout_4.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout_4.setObjectName(_fromUtf8("horizontalLayout_4"))
        self.verticalLayout_6 = QtGui.QVBoxLayout()
        self.verticalLayout_6.setObjectName(_fromUtf8("verticalLayout_6"))
        self.label_4 = QtGui.QLabel(self.centralwidget)
        self.label_4.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "size", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.verticalLayout_6.addWidget(self.label_4)
        self.sizelabel = QtGui.QLabel(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sizelabel.sizePolicy().hasHeightForWidth())
        self.sizelabel.setSizePolicy(sizePolicy)
        self.sizelabel.setText(_fromUtf8(""))
        self.sizelabel.setObjectName(_fromUtf8("sizelabel"))
        self.verticalLayout_6.addWidget(self.sizelabel)
        self.horizontalLayout_4.addLayout(self.verticalLayout_6)
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setContentsMargins(-1, 0, -1, -1)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.cov1label = QtGui.QLabel(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(10)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.cov1label.sizePolicy().hasHeightForWidth())
        self.cov1label.setSizePolicy(sizePolicy)
        self.cov1label.setText(_fromUtf8(""))
        self.cov1label.setObjectName(_fromUtf8("cov1label"))
        self.verticalLayout_2.addWidget(self.cov1label)
        self.cov2label = QtGui.QLabel(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.cov2label.sizePolicy().hasHeightForWidth())
        self.cov2label.setSizePolicy(sizePolicy)
        self.cov2label.setText(_fromUtf8(""))
        self.cov2label.setObjectName(_fromUtf8("cov2label"))
        self.verticalLayout_2.addWidget(self.cov2label)
        self.horizontalLayout_4.addLayout(self.verticalLayout_2)
        self.verticalLayout_8.addLayout(self.horizontalLayout_4)
        self.verticalLayout_7 = QtGui.QVBoxLayout()
        self.verticalLayout_7.setContentsMargins(-1, 0, -1, -1)
        self.verticalLayout_7.setObjectName(_fromUtf8("verticalLayout_7"))
        self.label_7 = QtGui.QLabel(self.centralwidget)
        self.label_7.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "mean", None, QtGui.QApplication.UnicodeUTF8))
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.verticalLayout_7.addWidget(self.label_7)
        self.meanlabel = QtGui.QLabel(self.centralwidget)
        self.meanlabel.setText(_fromUtf8(""))
        self.meanlabel.setObjectName(_fromUtf8("meanlabel"))
        self.verticalLayout_7.addWidget(self.meanlabel)
        self.verticalLayout_8.addLayout(self.verticalLayout_7)
        self.verticalLayout_4.addLayout(self.verticalLayout_8)
        self.horizontalLayout_6 = QtGui.QHBoxLayout()
        self.horizontalLayout_6.setContentsMargins(-1, 0, -1, -1)
        self.horizontalLayout_6.setObjectName(_fromUtf8("horizontalLayout_6"))
        self.label_3 = QtGui.QLabel(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_3.sizePolicy().hasHeightForWidth())
        self.label_3.setSizePolicy(sizePolicy)
        self.label_3.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "Blob Index(es):", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.horizontalLayout_6.addWidget(self.label_3)
        self.blobIndexlabel = QtGui.QLabel(self.centralwidget)
        self.blobIndexlabel.setText(_fromUtf8(""))
        self.blobIndexlabel.setObjectName(_fromUtf8("blobIndexlabel"))
        self.horizontalLayout_6.addWidget(self.blobIndexlabel)
        self.verticalLayout_4.addLayout(self.horizontalLayout_6)
        spacerItem1 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_4.addItem(spacerItem1)
        self.zoomOutButton = QtGui.QPushButton(self.centralwidget)
        self.zoomOutButton.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "zoom -", None, QtGui.QApplication.UnicodeUTF8))
        self.zoomOutButton.setObjectName(_fromUtf8("zoomOutButton"))
        self.verticalLayout_4.addWidget(self.zoomOutButton)
        self.zoomInButton = QtGui.QPushButton(self.centralwidget)
        self.zoomInButton.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "zoom +", None, QtGui.QApplication.UnicodeUTF8))
        self.zoomInButton.setObjectName(_fromUtf8("zoomInButton"))
        self.verticalLayout_4.addWidget(self.zoomInButton)
        self.horizontalLayout.addLayout(self.verticalLayout_4)
        MainWindowStructureAnalysis.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindowStructureAnalysis)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 900, 25))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setTitle(QtGui.QApplication.translate("MainWindowStructureAnalysis", "File", None, QtGui.QApplication.UnicodeUTF8))
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        self.menuReset = QtGui.QMenu(self.menubar)
        self.menuReset.setTitle(QtGui.QApplication.translate("MainWindowStructureAnalysis", "Reset", None, QtGui.QApplication.UnicodeUTF8))
        self.menuReset.setObjectName(_fromUtf8("menuReset"))
        self.menuOptions = QtGui.QMenu(self.menubar)
        self.menuOptions.setTitle(QtGui.QApplication.translate("MainWindowStructureAnalysis", "Options", None, QtGui.QApplication.UnicodeUTF8))
        self.menuOptions.setObjectName(_fromUtf8("menuOptions"))
        MainWindowStructureAnalysis.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindowStructureAnalysis)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindowStructureAnalysis.setStatusBar(self.statusbar)
        self.actionOpen_File = QtGui.QAction(MainWindowStructureAnalysis)
        self.actionOpen_File.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "Open File", None, QtGui.QApplication.UnicodeUTF8))
        self.actionOpen_File.setObjectName(_fromUtf8("actionOpen_File"))
        self.actionExit = QtGui.QAction(MainWindowStructureAnalysis)
        self.actionExit.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "Exit", None, QtGui.QApplication.UnicodeUTF8))
        self.actionExit.setObjectName(_fromUtf8("actionExit"))
        self.actionReset = QtGui.QAction(MainWindowStructureAnalysis)
        self.actionReset.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "Reset", None, QtGui.QApplication.UnicodeUTF8))
        self.actionReset.setObjectName(_fromUtf8("actionReset"))
        self.actionShow_covariance_matrix = QtGui.QAction(MainWindowStructureAnalysis)
        self.actionShow_covariance_matrix.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "Show covariance matrix", None, QtGui.QApplication.UnicodeUTF8))
        self.actionShow_covariance_matrix.setObjectName(_fromUtf8("actionShow_covariance_matrix"))
        self.actionCalculate_distance_matrix = QtGui.QAction(MainWindowStructureAnalysis)
        self.actionCalculate_distance_matrix.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "Calculate distance matrix", None, QtGui.QApplication.UnicodeUTF8))
        self.actionCalculate_distance_matrix.setObjectName(_fromUtf8("actionCalculate_distance_matrix"))
        self.actionSave_File = QtGui.QAction(MainWindowStructureAnalysis)
        self.actionSave_File.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "Save File", None, QtGui.QApplication.UnicodeUTF8))
        self.actionSave_File.setObjectName(_fromUtf8("actionSave_File"))
        self.actionSave_Features = QtGui.QAction(MainWindowStructureAnalysis)
        self.actionSave_Features.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "Save Features", None, QtGui.QApplication.UnicodeUTF8))
        self.actionSave_Features.setObjectName(_fromUtf8("actionSave_Features"))
        self.actionShow_Heatmatrix = QtGui.QAction(MainWindowStructureAnalysis)
        self.actionShow_Heatmatrix.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "show Heatmatrix", None, QtGui.QApplication.UnicodeUTF8))
        self.actionShow_Heatmatrix.setObjectName(_fromUtf8("actionShow_Heatmatrix"))
        self.actionToggle_colocalization_heatmap = QtGui.QAction(MainWindowStructureAnalysis)
        self.actionToggle_colocalization_heatmap.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "Toggle colocalization heatmap", None, QtGui.QApplication.UnicodeUTF8))
        self.actionToggle_colocalization_heatmap.setObjectName(_fromUtf8("actionToggle_colocalization_heatmap"))
        self.actionColocalization = QtGui.QAction(MainWindowStructureAnalysis)
        self.actionColocalization.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "Colocalization", None, QtGui.QApplication.UnicodeUTF8))
        self.actionColocalization.setObjectName(_fromUtf8("actionColocalization"))
        self.actionToggle_Images = QtGui.QAction(MainWindowStructureAnalysis)
        self.actionToggle_Images.setText(QtGui.QApplication.translate("MainWindowStructureAnalysis", "Toggle Images", None, QtGui.QApplication.UnicodeUTF8))
        self.actionToggle_Images.setObjectName(_fromUtf8("actionToggle_Images"))
        self.menuFile.addAction(self.actionOpen_File)
        self.menuFile.addAction(self.actionSave_File)
        self.menuFile.addAction(self.actionSave_Features)
        self.menuReset.addAction(self.actionReset)
        self.menuOptions.addAction(self.actionShow_covariance_matrix)
        self.menuOptions.addAction(self.actionCalculate_distance_matrix)
        self.menuOptions.addAction(self.actionShow_Heatmatrix)
        self.menuOptions.addAction(self.actionColocalization)
        self.menuOptions.addAction(self.actionToggle_colocalization_heatmap)
        self.menuOptions.addAction(self.actionToggle_Images)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuReset.menuAction())
        self.menubar.addAction(self.menuOptions.menuAction())

        self.retranslateUi(MainWindowStructureAnalysis)
        QtCore.QMetaObject.connectSlotsByName(MainWindowStructureAnalysis)

    def retranslateUi(self, MainWindowStructureAnalysis):
        pass

