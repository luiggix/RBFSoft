
/*------------------------------------------------------------------------
 *  Copyright (C) 2008  Luis M. de la Cruz
 *
 *  This file is part of TUNA::RBF
 *
 *  TUNA::RBF is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  TUNA::RBF is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ------------------------------------------------------------------------*/

#ifndef _TRAITS_H_
#define _TRAITS_H_

/*!
 ***************************************************************************
 *  \file Traits.hpp
 *  Definition of the variable types used in the whole TUNA::RBF.
 *  In principle, this header should be included in all files of TUNA::RBF.
 ***************************************************************************
 */
#include <ctime>
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
/*****
#include <flens/flens.h>
*****/

typedef double prec_t; ///< Definition of the precision for all calculations.

enum axis_t {X = 0, Y = 1, Z = 2};  // Axis
enum random_t {U = 0, R = 1, RU = 2}; // Random type
enum cartDir_t { UP = 0, DOWN = 1, LEFT = 2, RIGHT = 3 }; // Neighbors in Cart

const prec_t PI = 4.0 * atan(1.0);

/*!
 *  Tool for time measurements (borrowed from FLENS).
 */
struct timer
{
/*!
 *  Start the clock.
 */
    void tic() { _time = std::clock(); }
/*!
 *  Stop the clock.
 */ 
    prec_t toc() const { return prec_t(std::clock() - _time)/CLOCKS_PER_SEC; }
    std::clock_t _time;
};


/*******
// FLENS Matrix and vector declaration.
using namespace flens;
typedef DenseVector<Array<prec_t> > Vec;  ///< Dense vectors from FLENS.
typedef GeMatrix<FullStorage<prec_t, ColMajor> > Mat; ///< Dense matrix from FLENS.
typedef DenseVector<Array<int> > iVec;  ///< Integer dense vectors from FLENS.
typedef GeMatrix<FullStorage<int, ColMajor> > iMat; ///< Integer dense matrix from FLENS.
*******/
//typedef GeMatrix<FullStorage<prec_t, RowMajor> >  Mat;

using namespace Eigen;
typedef MatrixXd Mat; ///< Dense matrix from EIGEN.
typedef VectorXd Vec; ///< Dense vectors from EIGEN.
typedef MatrixXi iMat; ///< Integer dense matrix from EIGEN.
typedef VectorXi iVec; ///< Integer dense vectors from EIGEN.

#endif
