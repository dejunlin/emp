#ifndef  multifac_INC
#define  multifac_INC
/*
 * =====================================================================================
 *
 *       Filename:  multifac.hpp
 *
 *    Description:  compile-time multifactorial function 
 *
 *        Version:  1.0
 *        Created:  02/17/2015 10:49:18 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */
#include <type_traits>

/* ----------------------------------------------------------------------
  N*(N-Step)*(N-2*Step)*...*1 
------------------------------------------------------------------------- */
template <class TN, class TS> 
constexpr 
typename std::enable_if
<
  std::is_integral<typename std::common_type<TN, TS>::type>::value && 
  std::is_unsigned<typename std::common_type<TN, TS>::type>::value, 
  typename std::common_type<TN, TS>::type
>::type
multifac (const TN N, const TS Step) 
{
  return N < Step ? 1 : N*multifac(N-Step, Step);
}		/* -----  end of template function multifact  ----- */

/* ----------------------------------------------------------------------
  N*(N-Step)*(N-2*Step)*...*(N-(K-1)*Step) 
------------------------------------------------------------------------- */
template <class TN, class TS, class TK> 
constexpr 
typename std::enable_if
<
  std::is_integral<typename std::common_type<TN, TS, TK>::type>::value && 
  std::is_unsigned<typename std::common_type<TN, TS, TK>::type>::value, 
  typename std::common_type<TN, TS, TK>::type
>::type
multifac (const TN N, const TS Step, const TK K) 
{
  return K == 0 
            ? 1 
	    : (K == 1 
	          ? N 
		  : (N <= Step ? 0 : N * multifac(N-Step, Step, K-1)));
}		/* -----  end of template function multifact  ----- */

#endif   /* ----- #ifndef multifac_INC  ----- */


