# This file is part of BLASTnBRUSH
# Copyright (c) 2025 Aleksandra Liszka, Aleksandra Marcisz, Artur Stołowski 
# Licensed under the GPL v3.0 License

from PyQt6.QtCore import Qt

def move_checked_items(source_list, target_list):
    for i in range(source_list.count()):
        item = source_list.item(i)
        if item.checkState() == Qt.CheckState.Checked:
            text = item.text()
            if not any(target_list.item(j).text() == text for j in range(target_list.count())):
                target_list.addItem(text)
