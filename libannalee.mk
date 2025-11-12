################################################################################
#    This file is part of the NeHeP library.                                   #
#                                                                              #
#    Copyright (C) 1998-2002 Marko Grönroos <magi@iki.fi>                      #
#                                                                              #
################################################################################
#                                                                              #
#   This library is free software; you can redistribute it and/or              #
#   modify it under the terms of the GNU Library General Public                #
#   License as published by the Free Software Foundation; either               #
#   version 2 of the License, or (at your option) any later version.           #
#                                                                              #
#   This library is distributed in the hope that it will be useful,            #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          #
#   Library General Public License for more details.                           #
#                                                                              #
#   You should have received a copy of the GNU Library General Public          #
#   License along with this library; see the file COPYING.LIB.  If             #
#   not, write to the Free Software Foundation, Inc., 59 Temple Place          #
#   - Suite 330, Boston, MA 02111-1307, USA.                                   #
#                                                                              #
################################################################################

################################################################################
# Define root directory of the source tree
################################################################################
export SRCDIR ?= ..

################################################################################
# Define module name and compilation type
################################################################################
modname = libannalee
compile_library = 1
makedox         = 1

################################################################################
# Include build framework
################################################################################
include $(SRCDIR)/build/magicdef.mk

################################################################################
# Source files
################################################################################

sources =	anngenes.cc cangelosi.cc \
		kitano.cc layered.cc \
		learningenv.cc miller.cc nolfi.cc nolfinet.cc puredirect.cc \
		neat.cc

headers =	anngenes.h cangelosi.h cangelosinet.h \
		kitano.h layered.h learningenv.h miller.h nolfi.h nolfinet.h \
		neat.h

headersubdir =	annalee

################################################################################
# Recursively compile some subprojects
################################################################################
makemodules = projects

################################################################################
# Compile
################################################################################
include $(SRCDIR)/build/magiccmp.mk

