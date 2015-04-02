#ifndef  quaternion_INC
#define  quaternion_INC
/*
 * =====================================================================================
 *
 *       Filename:  quaternion.hpp
 *
 *    Description:  quaternion arithmetics
 *
 *        Version:  1.0
 *        Created:  03/30/2015 08:28:02 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */
#include <type_traits>

/**
* @brief convert a quaternion to a 3x3 rotation matrix
*
* @tparam Arr array type 
* @tparam T element type
* @param quat the quaternion
*
* @return the matrix
*/
template < template<class,size_t> class Arr, class T >
typename std::enable_if<
  std::is_arithmetic<T>::value,
  Arr<Arr<T,3>,3> 
>::type
quat_to_mat(const Arr<T, 4>& quat)
{
  double w2 = quat[0]*quat[0];
  double i2 = quat[1]*quat[1];
  double j2 = quat[2]*quat[2];
  double k2 = quat[3]*quat[3];
  double twoij = 2.0*quat[1]*quat[2];
  double twoik = 2.0*quat[1]*quat[3];
  double twojk = 2.0*quat[2]*quat[3];
  double twoiw = 2.0*quat[1]*quat[0];
  double twojw = 2.0*quat[2]*quat[0];
  double twokw = 2.0*quat[3]*quat[0];
  
  Arr<Arr<T,3>,3> ans;

  ans[0][0] = w2+i2-j2-k2;
  ans[0][1] = twoij-twokw;
  ans[0][2] = twojw+twoik;

  ans[1][0] = twoij+twokw;
  ans[1][1] = w2-i2+j2-k2;
  ans[1][2] = twojk-twoiw;

  ans[2][0] = twoik-twojw;
  ans[2][1] = twojk+twoiw;
  ans[2][2] = w2-i2-j2+k2;

  return ans;
}

template < template<class,size_t> class Arr, class T >
typename std::enable_if<
  std::is_arithmetic<T>::value,
  Arr<T,4> 
>::type
quatquat(const Arr<T,4>& a, const Arr<T,4>& b)
{
  Arr<T,4> ans;
  ans[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
  ans[1] = a[0]*b[1] + b[0]*a[1] + a[2]*b[3] - a[3]*b[2];
  ans[2] = a[0]*b[2] + b[0]*a[2] + a[3]*b[1] - a[1]*b[3];
  ans[3] = a[0]*b[3] + b[0]*a[3] + a[1]*b[2] - a[2]*b[1];
  return ans;
}

#endif   /* ----- #ifndef quaternion_INC  ----- */
