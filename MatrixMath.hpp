#ifndef  MatrixMath_INC
#define  MatrixMath_INC

/*
 * =====================================================================================
 *
 *       Filename:  MatrixMath.hpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  03/15/2015 07:11:03 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */
#include <cstdlib>
#include <array>
#include <utility>
#include "Typetraits.hpp"

template<class Elem, size_t N> using TVec = std::array<Elem,N>;
template<class Elem, size_t Nrow, size_t Ncol> using TMatrix = std::array<std::array<Elem,Ncol>,Nrow>;

/**
* @brief multiply a vector by a scalar
*
* @tparam sElem element type of s
* @tparam aElem element type of a
* @tparam N size of the vector
* @param s scalar
* @param a vector
*
* @return vector scaled by s
*/
template<class sElem, class aElem, size_t N>
typename std::enable_if
<
  std::is_arithmetic<sElem>::value && std::is_arithmetic<aElem>::value, 
  TVec<decltype(std::declval<sElem>()*std::declval<aElem>()),N>
>::type 
operator*(const sElem& s, const TVec<aElem,N>& a) {
  using T = decltype(std::declval<sElem>()*std::declval<aElem>());
  TVec<T,N> ans;
  for(size_t i = 0; i < N; ++i) {
    ans[i] = s*a[i];
  }
  return ans;
}

/**
* @brief multiply a scalar by a vector 
*
* @tparam aElem element type of a
* @tparam sElem element type of s
* @tparam N size of the vector
* @param s scalar
* @param a vector
*
* @return vector scaled by s
*/
template<class aElem, class sElem, size_t N>
typename std::enable_if
<
  std::is_arithmetic<sElem>::value && std::is_arithmetic<aElem>::value, 
  TVec<decltype(std::declval<sElem>()*std::declval<aElem>()),N>
>::type 
operator*(const TVec<aElem,N>& a, const sElem& s) 
{
  return operator*(s,a);
}

/**
* @brief translate a vector by another vector
* @details only if bElem can be promoted to aElem
* @tparam aElem element type of a
* @tparam bElem element type of b
* @tparam N size of the vectors
* @param a the vector to be translated
* @param b the amount to translate 
*
* @return the translated vector 
*/
template<class aElem, class bElem, size_t N>
typename std::enable_if<
  std::is_arithmetic<aElem>::value && std::is_arithmetic<bElem>::value &&  
  are_same<aElem, decltype(std::declval<aElem>()+std::declval<bElem>())>::value,
  TVec<aElem,N>&>::type 
operator+=(TVec<aElem,N>& a, const TVec<bElem,N>& b) {
  for(size_t i = 0; i < N; ++i) {
    a[i] += b[i];
  }
  return a;
}

/**
* @brief sum of two vectors
*
* @tparam aElem element type of a
* @tparam bElem element type of b
* @tparam N size of the vectors
* @param a left-hand operand
* @param b right-hand operand
*
* @return sum of the two input vectors
*/
template<class aElem, class bElem, size_t N>
typename std::enable_if
<
  std::is_arithmetic<aElem>::value && std::is_arithmetic<bElem>::value, 
  TVec<decltype(std::declval<aElem>()+std::declval<bElem>()),N> 
>::type 
operator+(const TVec<aElem,N>& a, const TVec<bElem,N>& b) {
  using T = decltype(std::declval<aElem>()+std::declval<bElem>());
  TVec<T,N> ans;
  for(size_t i = 0; i < N; ++i) {
    ans[i] = a[i] + b[i];
  }
  return ans;
}

/**
* @brief difference between two vectors
*
* @tparam aElem element type of a
* @tparam bElem element type of b
* @tparam N size of the vectors
* @param a left-hand operand
* @param b right-hand operand
*
* @return difference between the two vectors 
*/
template<class aElem, class bElem, size_t N>
typename std::enable_if
<
  std::is_arithmetic<aElem>::value && std::is_arithmetic<bElem>::value, 
  TVec<decltype(std::declval<aElem>()-std::declval<bElem>()),N> 
>::type 
operator-(const TVec<aElem,N>& a, const TVec<bElem,N>& b) {
  using T = decltype(std::declval<aElem>()-std::declval<bElem>());
  TVec<T,N> ans;
  for(size_t i = 0; i < N; ++i) {
    ans[i] = a[i] - b[i];
  }
  return ans;
}

/**
* @brief dot product of two vectors
*
* @tparam aElem element type of a
* @tparam bElem element type of b
* @tparam N size of the vectors
* @param a left-hand operand
* @param b right-hand operand
*
* @return the dot product
*/
template<class aElem, class bElem, size_t N>
typename std::enable_if
<
  std::is_arithmetic<aElem>::value && std::is_arithmetic<bElem>::value, 
  decltype(std::declval<aElem>()*std::declval<bElem>()) 
>::type 
operator*(const TVec<aElem,N>& a, const TVec<bElem,N>& b) {
  using T = decltype(std::declval<aElem>()*std::declval<bElem>());
  T ans{0};
  for(size_t i = 0; i < N; ++i) {
    ans += a[i]*b[i];
  }
  return ans;
}

/**
* @brief cross product of two vectors of 3
*
* @tparam aElem element type of a
* @tparam bElem element type of b
* @param a left-hand operand
* @param b right-hand operand
*
* @return the cross product
*/
template<class aElem, class bElem>
typename std::enable_if
<
  std::is_arithmetic<aElem>::value && std::is_arithmetic<bElem>::value, 
  TVec<decltype(std::declval<aElem>()*std::declval<bElem>()),3> 
>::type 
operator^(const TVec<aElem,3>& a, const TVec<bElem,3>& b) {
  using T = TVec<decltype(std::declval<aElem>()*std::declval<bElem>()),3>;
  const T ans = {
    a[1]*b[2] - a[2]*b[1],
    a[2]*b[0] - a[0]*b[2],
    a[0]*b[1] - a[1]*b[0],
  };
  return ans;
}

/**
* @brief element-wise product of two vectors
*
* @tparam aElem element type of a
* @tparam bElem element type of b
* @tparam N size of the vectors
* @param a left-hand operand
* @param b right-hand operand
*
* @return the dot product
*/
template<class aElem, class bElem, size_t N>
typename std::enable_if
<
  std::is_arithmetic<aElem>::value && std::is_arithmetic<bElem>::value, 
  TVec<decltype(std::declval<aElem>()*std::declval<bElem>()),N> 
>::type 
operator&(const TVec<aElem,N>& a, const TVec<bElem,N>& b) {
  using T = decltype(std::declval<aElem>()*std::declval<bElem>());
  TVec<T,N> ans;
  for(size_t i = 0; i < N; ++i) {
    ans[i] = a[i]*b[i];
  }
  return ans;
}

/**
* @brief element-wise product assignment 
*
* @tparam aElem element type of a
* @tparam bElem element type of b
* @tparam N size of the vectors
* @param a left-hand operand
* @param b right-hand operand
*
* @return the left operand
*/
template<class aElem, class bElem, size_t N>
typename std::enable_if
<
  std::is_arithmetic<aElem>::value && std::is_arithmetic<bElem>::value &&
  are_same<aElem, decltype(std::declval<aElem>()*std::declval<bElem>())>::value, 
  TVec<aElem,N>& 
>::type 
operator&=(TVec<aElem,N>& a, const TVec<bElem,N>& b) {
  for(size_t i = 0; i < N; ++i) {
    a[i] *= b[i];
  }
  return a;
}


/**
* @brief multiply a row-major matrix by a row-vector
* @details This always requires the vector is multiplied to the LEFT of the matrix
* @tparam vElem element type of v
* @tparam mElem element type of m
* @tparam Nrow number of rows
* @tparam Ncol number of columns
* @param v vector
* @param m matrix
*
* @return a vector with each element being the dot product between the input vector and one colume of the matrix
*/
template<class vElem, class mElem, size_t Nrow, size_t Ncol>
typename std::enable_if
<
  std::is_arithmetic<vElem>::value && std::is_arithmetic<mElem>::value,
  TVec<decltype(std::declval<vElem>()*std::declval<mElem>()),Ncol>
>::type 
operator*(const TVec<vElem,Nrow>& v, const TMatrix<mElem,Nrow,Ncol>& m)
{
  using T = decltype(std::declval<vElem>()*std::declval<mElem>());
  TVec<T,Ncol> ans{};
  for(size_t i = 0; i < Nrow; ++i) {
    ans += v[i] * m[i];
  }
  return ans;
}

/**
* @brief multiply a row-major matrix by a column-vector
* @details This always requires the vector is multiplied to the RIGHT of the matrix
* @tparam vElem element type of v
* @tparam mElem element type of m
* @tparam Nrow number of rows
* @tparam Ncol number of columns
* @param v vector
* @param m matrix
*
* @return a vector with each element being the dot product between the input vector and one row of the matrix
*/
template<class mElem, class vElem, size_t Nrow, size_t Ncol>
typename std::enable_if
<
  std::is_arithmetic<mElem>::value && std::is_arithmetic<vElem>::value,
  TVec<decltype(std::declval<mElem>()*std::declval<vElem>()),Nrow>
>::type 
operator*(const TMatrix<mElem,Nrow,Ncol>& m, const TVec<vElem,Ncol>& v)
{
  using T = decltype(std::declval<mElem>()*std::declval<vElem>());
  TVec<T,Nrow> ans;
  for(size_t i = 0; i < Nrow; ++i) {
    ans[i] = m[i] * v;
  }
  return ans;
}

/**
* @brief 2-fold contraction of 2 matrices
*
* @tparam Arr array type
* @tparam Elem element type
*/
template <template<class,size_t> class Arr, class ... Elem> struct MatrixDoubleDotProduct;

template <class Elem>
struct MatrixDoubleDotProduct<std::array, Elem> {
  template<size_t Nrow, size_t Ncol>
  typename std::enable_if<std::is_arithmetic<Elem>::value, Elem>::type
  operator()(const TMatrix<Elem,Nrow,Ncol>& m1, const TMatrix<Elem,Nrow,Ncol>& m2) {
    Elem ans{0};
    for(size_t i = 0; i < Nrow; ++i) {
      ans += m1[i] * m2[i];
    }
    return ans;
  }
};

template <class Elem1, class Elem2>
struct MatrixDoubleDotProduct<std::array, Elem1, Elem2> {
  template<size_t Nrow, size_t Ncol>
  typename std::enable_if<
    std::is_arithmetic<Elem1>::value && std::is_arithmetic<Elem2>::value,
    decltype(std::declval<Elem1>()*std::declval<Elem2>()) 
  >::type
  operator()(const TMatrix<Elem1,Nrow,Ncol>& m1, const TMatrix<Elem2,Nrow,Ncol>& m2) {
    using T = decltype(std::declval<Elem1>()*std::declval<Elem2>());
    T ans{0};
    for(size_t i = 0; i < Nrow; ++i) {
      ans += m1[i] * m2[i];
    }
    return ans;
  }
};

/**
* @brief matrix-matrix multiplication
*
* @tparam Elem elemnt type 
* @tparam Nrow1 the number of rows of the left-hand operand
* @tparam Ncol1 the number of columns of the left-hand operand or the number of rows of the right-hand operand
* @tparam Ncol2 the number of columns of the right-hand operand
* @param m1 left-hand operand
* @param m2 right-hand operand
*
* @return product of the two matrices 
*/
template<class Elem1, class Elem2, size_t Nrow1, size_t Ncol1, size_t Ncol2>
typename std::enable_if<
  std::is_arithmetic<Elem1>::value && std::is_arithmetic<Elem2>::value,
  TMatrix<decltype(std::declval<Elem1>()*std::declval<Elem2>()),Nrow1,Ncol2> 
>::type
operator*(const TMatrix<Elem1,Nrow1,Ncol1>& m1, const TMatrix<Elem2,Ncol1,Ncol2>& m2) {
  using T = decltype(std::declval<Elem1>()*std::declval<Elem2>());
  TMatrix<T,Nrow1,Ncol2> ans;
  for(size_t i = 0; i < Nrow1; ++i) {
    ans[i] = m1[i] * m2;
  }
  return ans;
}

template<class Elem, size_t N>
typename std::enable_if<
  std::is_arithmetic<Elem>::value,
  Elem
>::type
sum(const TVec<Elem,N>& m) {
  Elem ans{0};
  for(size_t i = 0; i < N; ++i) {
    ans += m[i];
  }
  return ans;
}


/**
* @brief sum all the elements of a matrix
*
* @tparam Elem element type
* @tparam Nrow number of row
* @tparam Ncol number of column
* @param m input matrix
*
* @return the sum
*/
template<class Elem, size_t Nrow, size_t Ncol>
typename std::enable_if<
  std::is_arithmetic<Elem>::value,
  Elem
>::type
sum(const TMatrix<Elem,Nrow,Ncol>& m) {
  Elem ans{0};
  for(size_t i = 0; i < Nrow; ++i) {
    ans += sum(m[i]);
  }
  return ans;
}

/**
* @brief sum the rows of a matrix
*
* @tparam Elem element type
* @tparam Nrow number of row
* @tparam Ncol number of column
* @param m input matrix
*
* @return the sum
*/
template<class Elem, size_t Nrow, size_t Ncol>
typename std::enable_if<
  std::is_arithmetic<Elem>::value,
  TVec<Elem,Ncol> 
>::type
sumrows(const TMatrix<Elem,Nrow,Ncol>& m) {
  TVec<Elem,Ncol> ans{};
  for(size_t i = 0; i < Nrow; ++i) {
    ans += m[i];
  }
  return ans;
}
#endif   /* ----- #ifndef MatrixMath_INC  ----- */
