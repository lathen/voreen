/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2006-2008 Visualization and Computer Graphics Group, *
 * Department of Computer Science, University of Muenster, Germany.   *
 * <http://viscg.uni-muenster.de>                                     *
 *                                                                    *
 * This file is part of the tgt library. This library is free         *
 * software; you can redistribute it and/or modify it under the terms *
 * of the GNU Lesser General Public License version 2.1 as published  *
 * by the Free Software Foundation.                                   *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       *
 * GNU Lesser General Public License for more details.                *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License in the file "LICENSE.txt" along with this library.         *
 * If not, see <http://www.gnu.org/licenses/>.                        *
 *                                                                    *
 **********************************************************************/

#ifndef TGT_INIT_H
#define TGT_INIT_H

#include "tgt/config.h"

namespace tgt {

class InitFeature {
public:
    enum Features {
        NONE            =      0,
        LOG_MANAGER     = 1 << 0,
        FILE_SYSTEM     = 1 << 1,
        GPU_PROPERTIES  = 1 << 2,
        SCRIPT_MANAGER  = 1 << 3,
        SHADER_MANAGER  = 1 << 4,
        MODEL_MANAGER   = 1 << 5,
        TEXTURE_MANAGER = 1 << 6,
        TESSELATOR      = 1 << 7,
        LOG_TO_CONSOLE  = 1 << 30,

        NO_SHADER_CACHING  = 1 << 31,

        ALL             = ( 
                            LOG_MANAGER | FILE_SYSTEM | GPU_PROPERTIES | 
                            SCRIPT_MANAGER | SHADER_MANAGER |
                            MODEL_MANAGER | TEXTURE_MANAGER | TESSELATOR |
                            LOG_TO_CONSOLE | NO_SHADER_CACHING                            
                            )        
    };
};

/// init general purpose singletons of tgt
void init(InitFeature::Features featureset = InitFeature::ALL);
/// init GLEW and OpenGL-dependent singletons of tgt.
/// to be called when OpenGL context already exists.
void initGL(InitFeature::Features featureset = InitFeature::ALL);

/// deinit the singletons of tgt
void deinit();
/// deinit the singletons of tgt
void deinitGL();

};

#endif //TGT_INIT_H
