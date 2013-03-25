/*
 * Copyright (C) 2011 Joachim Schleicher <J.Schleicher@stud.uni-heidelberg.de>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */

#include "filenamelineedit.h"
#include "stormparams.h"

#include <QDragEnterEvent>
#include <QDropEvent>
#include <QUrl>
#include <QList>
#include <QString>
#include <QFileInfo>
#include <QValidator>

FilenameLineEdit::FilenameLineEdit(QWidget * parent)
: QLineEdit(parent)
{
    setAcceptDrops(true);
}

FilenameLineEdit::~FilenameLineEdit()
{
}

void FilenameLineEdit::dragEnterEvent(QDragEnterEvent* event)
{
    // accept just text/uri-list mime format
    if (event->mimeData()->hasFormat("text/uri-list") && event->mimeData()->hasUrls())
    {
        QString fName = matchingFile(event->mimeData()->urls());
        if (!fName.isEmpty()) {
            event->acceptProposedAction();
        }
    }
}

void FilenameLineEdit::dropEvent(QDropEvent* event)
{
    QString fName = matchingFile(event->mimeData()->urls());
    if (!fName.isEmpty()) {
        setText(fName);
        event->acceptProposedAction();
        emit textEdited(fName);
    }
}

QString FilenameLineEdit::matchingFile(const QList<QUrl> &urlList)
{
    const QValidator *valid = validator();
    for (const QUrl &url : urlList) {
        QString fName = urlList[0].toLocalFile();
        QFileInfo fInfo(fName);
        int length = fName.size();
        if (fInfo.isFile() && fInfo.exists() && (!valid || valid->validate(fName, length) == QValidator::Acceptable)) {
            return fName;
        }
    }
    return QString();
}
