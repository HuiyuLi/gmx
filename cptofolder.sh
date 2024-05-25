#!/bin/bash

# Copyright (c) 2023, Huiyu Li and Ao Ma, University of Illinois at Chicago.
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

dir=gromacs-2019.2-code
cp ${dir}/readir.cpp ./src/gromacs/gmxpreprocess/
cp ${dir}/tpxio.cpp ./src/gromacs/fileio/
cp ${dir}/inputrec.h ./src/gromacs/mdtypes/
cp ${dir}/grompp.cpp ./src/gromacs/gmxpreprocess
cp ${dir}/gen_maxwell_velocities.h ./src/gromacs/gmxpreprocess
cp ${dir}/gen_maxwell_velocities.cpp ./src/gromacs/gmxpreprocess
cp ${dir}/md.cpp ./src/gromacs/mdrun/
cp ${dir}/mdoutf.cpp ./src/gromacs/mdlib/
cp ${dir}/mdoutf.h ./src/gromacs/mdlib/

