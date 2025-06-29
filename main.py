# This file is part of TidyTool
# Copyright (c) 2025 Aleksandra Liszka, Artur Stołowski, Aleksandra Marcisz
# Licensed under the MIT License

import sys
sys.path.append('./libraries') 
sys.path.append('./resources')
from PyQt6.QtWidgets import QApplication
from gui import MainWindow

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.showMaximized()
    sys.exit(app.exec())
