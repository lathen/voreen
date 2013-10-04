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

#ifndef VRN_LINKEVALUATOSPATIALTRANSFUNCRID_H
#define VRN_LINKEVALUATOSPATIALTRANSFUNCRID_H

#include "voreen/core/properties/link/linkevaluatoridgeneric.h"

namespace voreen {

class VRN_CORE_API LinkEvaluatorSpatialTransFuncId : public LinkEvaluatorBase {
public:
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorSpatialTransFuncId(); }
    virtual std::string getClassName() const  { return "LinkEvaluatorSpatialTransFuncId"; }
    virtual std::string getGuiName() const    { return "Identity"; }

    virtual void eval(Property* src, Property* dst) throw (VoreenException);

    bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

} // namespace

#endif // VRN_LINKEVALUATOSPATIALTRANSFUNCRID_H
