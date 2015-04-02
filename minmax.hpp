#ifndef  minmax_INC
#define  minmax_INC

/*
 * =====================================================================================
 *
 *       Filename:  minmax.hpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  03/17/2015 10:32:34 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */

template < class T >
constexpr const T& Min(const T& a, const T& b) {
  return a < b ? a : b;
}

template < class T >
constexpr const T& Max(const T& a, const T& b) {
  return a > b ? a : b;
}

#endif   /* ----- #ifndef minmax_INC  ----- */
