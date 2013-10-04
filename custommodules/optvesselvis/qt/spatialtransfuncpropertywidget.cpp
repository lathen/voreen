/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2013 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

#include "spatialtransfuncpropertywidget.h"

#include "voreen/qt/widgets/voreentoolwindow.h"
#include "spatialtransfuncplugin.h"
#include "voreen/core/processors/processor.h"
#include "../properties/spatialtransfuncproperty.h"

#include <QPushButton>

namespace voreen {

SpatialTransFuncPropertyWidget::SpatialTransFuncPropertyWidget(SpatialTransFuncProperty* prop, QWidget* parent)
    : QPropertyWidgetWithEditorWindow(prop, parent, false)
    , plugin_(0)
    , property_(prop)
    , editBt_(new QPushButton(tr("edit")))
{

    tgtAssert(prop, "No property passed");

    if (!prop->getLazyEditorInstantiation() || editorVisibleOnStartup())
        createEditorWindow(Qt::RightDockWidgetArea);

    addWidget(editBt_);

    connect(editBt_, SIGNAL(clicked()), this, SLOT(setProperty()));
    connect(editBt_, SIGNAL(clicked()), this, SIGNAL(widgetChanged()));

    addVisibilityControls();
    QFontInfo fontInfo(font());
    editBt_->setFont(QFont(fontInfo.family(), QPropertyWidget::fontSize_));
}

void SpatialTransFuncPropertyWidget::updateFromPropertySlot() {
    if (plugin_)
        plugin_->updateFromProperty();
}

void SpatialTransFuncPropertyWidget::setProperty() {
    if (!disconnected_) {
        // lazy instantiation of transfunc editor window
        if (!editorWindow_) {
            createEditorWindow(Qt::RightDockWidgetArea);
            tgtAssert(editorWindow_, "Transfunc editor not instantiated");
        }

        if (editorWindow_->isVisible()) {
            //close widget
            editorWindow_->close();
        }
        else {
            //open Widget
            editorWindow_->showNormal();
            plugin_->updateFromProperty();
        }
    }
}

void SpatialTransFuncPropertyWidget::disconnect() {
    disconnected_ = true;
    if (plugin_)
        plugin_->disconnect();
}

QWidget* SpatialTransFuncPropertyWidget::createEditorWindowWidget() {
    plugin_ = new SpatialTransFuncPlugin(property_, parentWidget(), Qt::Horizontal);
    plugin_->createWidgets();
    plugin_->createConnections();
    connect(plugin_, SIGNAL(transferFunctionChanged()), this, SIGNAL(modified()));

    return plugin_;
}

void SpatialTransFuncPropertyWidget::customizeEditorWindow() {
    editorWindow_->setAllowedAreas(Qt::RightDockWidgetArea);
    editorWindow_->setFloating(true);
}

Property* SpatialTransFuncPropertyWidget::getProperty() {
    return property_;
}

} // namespace voreen
