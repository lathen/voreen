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

#include "linkevaluatorspatialtransfuncid.h"
#include "../spatialtransfuncproperty.h"
#include "voreen/core/datastructures/transfunc/transfunc.h"

namespace voreen {

void LinkEvaluatorSpatialTransFuncId::eval(Property* src, Property* dst) throw (VoreenException) {
    SpatialTransFuncProperty* dstCast = static_cast<SpatialTransFuncProperty*>(dst);
    SpatialTransFuncProperty* srcCast = static_cast<SpatialTransFuncProperty*>(src);

    if(!srcCast->get()) {
        LERRORC("voreen.LinkEvaluatorSpatialTransFuncId", "src is has no TF");
        return;
    }

    TransFunc* tf = srcCast->get()->clone();

    dstCast->set(tf);
    dstCast->setKeysLocked(srcCast->getKeysLocked());
    dstCast->setParamMagnitude(srcCast->getParamMagnitude());
    dstCast->setWindowMode(srcCast->getWindowMode());
}

bool LinkEvaluatorSpatialTransFuncId::arePropertiesLinkable(const Property* p1, const Property* p2) const {
    tgtAssert(p1, "null pointer");
    tgtAssert(p2, "null pointer");

    return (dynamic_cast<const SpatialTransFuncProperty*>(p1) && dynamic_cast<const SpatialTransFuncProperty*>(p2));
}

} // namespace
