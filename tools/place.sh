#!/bin/bash
#
# place the GNU header in a file
#

x=$1

HEAD="!**************************************************************************************
!    Copyright 2013-2017 Holger Kruse                                                 *
!                                                                                     *
!    This file is part of Xopt.                                                       *
!                                                                                     *
!    Xopt is free software: you can redistribute it and/or modify                     *
!    it under the terms of the GNU Lesser General Public License as published by      *
!    the Free Software Foundation, either version 3 of the License, or                *
!    (at your option) any later version.                                              *
!                                                                                     *
!    Xopt is distributed in the hope that it will be useful,                          *
!    but WITHOUT ANY WARRANTY; without even the implied warranty of                   *
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                    *
!    GNU Lesser General Public License for more details.                              *
!                                                                                     *
!    You should have received a copy of the GNU Lesser General Public License         *
!    along with xopt.  If not, see <http://www.gnu.org/licenses/>.                    *
!                                                                                     *
!                                                                                     *
!  Feel free to contact me: Holger Kruse (mail2holger@gmail.com)                      *
!                                                                                     *
!**************************************************************************************
~                                                                                          "


for x in $1; do
cp $x .x.tmp
echo $HEAD > $x
cat .x.tmp >> $x
rm .x.tmp
done
