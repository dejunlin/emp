#ifndef  emp_INC
#define  emp_INC
/*
 * =====================================================================================
 *
 *       Filename:  emp.hpp
 *
 *    Description:  Electric MultiPole tensor class
 *
 *        Version:  1.0
 *        Created:  02/16/2015 10:21:50 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */
#include <array>
#include <type_traits>
#include "../generic/Arithmetics/multinom.hpp"
using namespace std;

template <class UINT, UINT R, class FL, bool Detraced>
class EMP {
  static_assert(is_integral<UINT>::value && is_unsigned<UINT>::value, "Rank of EMP must be unsigned integral");
  static_assert(is_floating_point<FL>::value, "EMP requires floating point");
  //! Number of unique elements of the EMP moments
  static constexpr UINT N = binomial(R+UINT(2), UINT(2));
  //! Maximal rank of trace applicable to this EMP moments
  static constexpr UINT Ntr = R % UINT(2) ? (R-UINT(1)/UINT(2) : R/UINT(2);
  //! The EMP moments array
  array<FL,N> _EMP;
  //! The coefficient matrix that holds the coefficients in EMP interaction
  /**
  * C[i][j] is the coefficient  
  */
  array<FL,Ntr+1> 
};

#endif   /* ----- #ifndef emp_INC  ----- */
