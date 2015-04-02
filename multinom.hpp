#ifndef  multinom_INC
#define  multinom_INC

/*
 * =====================================================================================
 *
 *       Filename:  multinom.hpp
 *
 *    Description:  Multinomial coefficients
 *
 *        Version:  1.0
 *        Created:  02/17/2015 04:10:58 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */

#include <type_traits>
#include "multifac.hpp"

/**
* @brief N choose K
*
* @tparam TN type of N
* @tparam TK type of K
* @param N N
* @param K K
*
* @return N!/( (N-K)!*K! )
*/
template <class TN, class TK>
constexpr 
typename std::enable_if
<
  std::is_integral<typename std::common_type<TN, TK>::type>::value && 
  std::is_unsigned<typename std::common_type<TN, TK>::type>::value, 
  typename std::common_type<TN, TK>::type
>::type
binomial ( const TN& N, const TK& K )
{
  return multifac(N, 1, K)/multifac(K, 1) ;
}		/* -----  end of template function binomial  ----- */

/**
* @brief Trinomial coefficient 
*
* @tparam TN Type of N
* @tparam TK1 Type of K1
* @tparam TK2 Type of K2
* @param N N
* @param K1 K1 
* @param K2 K2
*
* @return N/((N-K1-K2)!*K1!*K2!)
*/
template <class TN, class TK1, class TK2>
constexpr auto trinomial ( const TN& N, const TK1& K1, const TK2& K2 ) 
-> typename std::enable_if
<
  std::is_integral<typename std::common_type<TN, TK1, TK2>::type>::value && 
  std::is_unsigned<typename std::common_type<TN, TK1, TK2>::type>::value,
  typename std::common_type<TN, TK1, TK2>::type
>::type
{
  return multifac(N, 1, K1+K2)/(multifac(K1, 1)*multifac(K2,1));
}		/* -----  end of template function binomial  ----- */

/**
* @brief The number of ways to choose K pairs from N ( N - 2*K >= 0 )
*
* @tparam TN Type of N
* @tparam TK Type of K
* @param N N
* @param K K
*
*/
template <class TN, class TK>
constexpr 
typename std::enable_if
<
  std::is_integral<typename std::common_type<TN, TK>::type>::value && 
  std::is_unsigned<typename std::common_type<TN, TK>::type>::value, 
  typename std::common_type<TN, TK>::type
>::type
NpairK ( const TN& N, const TK& K )
{
  using T = typename std::common_type<TN, TK>::type;
  return multifac(N, 1, 2*K) / ( (T{1} << K) * multifac(K,1) );
}		

#endif   /* ----- #ifndef multinom_INC  ----- */
