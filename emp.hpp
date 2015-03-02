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
#include <vector>
#include <tuple>
#include <type_traits>
#include "../generic/Arithmetics/multinom.hpp"
#include "../generic/Arithmetics/NDArray.hpp"
#include "../generic/Arithmetics/Tuple.hpp"
#include "../generic/Arithmetics/indexseq_generator.hpp"
using namespace std;

/**
* @brief EMP moment class
*
* @tparam UINT unsigned integral type of the rank of the EMP moment
* @tparam R Rank of the EMP moment
* @tparam FL floating point type for the EMP moment
* @tparam Detraced if the moment is detraced
*/
template <class UINT, UINT R, class FL, bool Detraced>
class EMP {
  public:

  static_assert(is_integral<UINT>::value && is_unsigned<UINT>::value, "Rank of EMP must be unsigned integral");
  static_assert(is_floating_point<FL>::value, "EMP requires floating point");
  //! Number of unique elements of the EMP moments
  static constexpr UINT N = binomial(R+2u, 2u);
  //! Maximal rank of trace applicable to this EMP moments
  static constexpr UINT Maxlm = Detraced ? 0 : ( R % 2 ? (R-1)/2 : R/2 );
  //! Min lc
  static constexpr UINT Minlc = 0;
  //! Max lc
  static constexpr UINT Maxlc = R;
  //! Initialize the EMP moment array EMP::_emp to zero
  template<size_t lm> struct init_emp {
    //! The rank after taking lm-fold trace
    static constexpr UINT rank = R-2*lm;
    static constexpr UINT n = binomial(rank+2u, 2u);
    using Twrapper = NDArray<array, init_to_zeros, FL, IndexSeq<n>>;
    using type = typename Twrapper::type;
    static constexpr type value = Twrapper::value; 
  };
  
  //! The EMP moment array and its traces
  using tmaker_for_emp = Tuple_Maker<tuple, init_emp, container_index<Maxlm+1>>;
  using tuple_for_emp = typename tmaker_for_emp::type;
  tuple_for_emp _emp = tmaker_for_emp::value;
  
  //! Initialize a tuple of 2d arrays to zero in order to store the member array EMP::A
  template<size_t lc> struct init_A {
    //maxlm is the maximum of lm given lc
    static constexpr UINT maxlm = (R - lc) % 2 ? (R-lc-1) / 2 : (R-lc)/2;
    static constexpr UINT Nx = maxlm + 1;
    static constexpr UINT Ny = binomial(lc+2u, 2u);
    using Twrapper = NDArray<array, init_to_zeros, FL, IndexSeq<Nx, Ny>>;
    using type = typename Twrapper::type; 
    static constexpr type value = Twrapper::value;
  };
  using tmaker_for_A = Tuple_Maker<tuple, init_A, container_index<Maxlc+1>>;
  using tuple_for_A = typename tmaker_for_A::type; 
  //! Array holding A^(m:lm) .Om. r^(Om), where Om = m-lm-lc, for each lc (first dimension) and lm (second dimension)
  tuple_for_A A = tmaker_for_A::value; 

  /**
  * @brief constructor 
  *
  * @param emp EMP moment array in global frame
  */
  EMP(const array<FL, N>& emp) {
    get<0>(_emp) = emp;
  };

  
};

#endif   /* ----- #ifndef emp_INC  ----- */
