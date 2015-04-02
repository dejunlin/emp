#ifndef  SymmetricCartesianTensor_INC
#define  SymmetricCartesianTensor_INC

/*
 * =====================================================================================
 *
 *       Filename:  SymmetricCartesianTensor.hpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  03/02/2015 10:18:29 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */
#include "multinom.hpp"
#include "multifac.hpp"
#include "NDArray.hpp"
#include "indexseq_generator.hpp"
#include "MatrixMath.hpp"
#include "Tuple.hpp"
#include "quaternion.hpp"
#include <type_traits>
#include <initializer_list>
#include <cmath>
#include <functional>
#include <tuple>
#include <algorithm>
#include <utility>

/**
* @brief Arithmetics of totally symmetric Cartesian tensor
*
* @tparam Arr Arr template class
* @tparam Elem Type of underlying element
* @tparam R Rank of the tensor
* @tparam Detraced If a detracer operator (eq 5.1 in Applequist, J., J. of Phys. A: Math. Gen., 1989, 22, 4303) has been applied to the tensor
*/
template < template<class, size_t> class Arr, class Elem, size_t R, bool Detraced > 
struct SymmetricCartesianTensor {
  static_assert( !Detraced || R >= 2, "Tensor whose rank < 2 can't be detraced"); 
  template <size_t Rank> using ThisType = SymmetricCartesianTensor<Arr, Elem, Rank, Detraced >;
  template <size_t Rank, bool ifDetraced> using ThisType_tr = SymmetricCartesianTensor<Arr, Elem, Rank, ifDetraced >;
  //! just short hand for matrices
  template<size_t Nrow, size_t Ncol> using TeMatrix = Arr<Arr<Elem,Ncol>,Nrow>;
  //! use this to access nested members/routines/classes that doesn't depend on the trace
  template <size_t Rank> using ThisType_igtr = ThisType_tr<Rank,false>;
  //! Number of unique elements
  static constexpr size_t N = binomial(R+2u, 2u);
  //! short hand for R
  static constexpr size_t _rank = R;
  //! if detraced
  static constexpr bool _detraced = Detraced;
  //! check if a size N (e.g., size of an array) can be treated as a tensor
  template<size_t S>
  struct isTensorSize {
    private:
    //! can't use _Rank directly unless value is true
    static constexpr size_t _Rank = ceil((sqrt(8*S+1)-3)/2);
    public:
    static constexpr bool value = (binomial(_Rank+2u,2u) == S);
    static constexpr typename std::enable_if<value,size_t>::type Rank = _Rank;
  };

  using Tinit = const Arr<Elem,N>&;
  using Trefelem = std::reference_wrapper<const Elem>;
  //! Array of all the unique elements
  Arr<Elem, N> _data;
  
  //! default constructor
  SymmetricCartesianTensor() : _data( NDArray<Arr, init_to_zeros, Elem, IndexSeq<N>>::value ) {};
  //! initializer list constructor (should take precedence over the forwarding constructor)
  SymmetricCartesianTensor(const std::initializer_list<Elem>& list) {
    std::copy_n(list.begin(), std::min(N, list.size()), _data.begin());
  };
  //! forwarding constructor (To avoid it taking precedence over the copy constructor or the Tinit constructor, we use 
  //SFINAE to enable this IFF. there's at least one argument is not of type ThisType<R> or Tinit)
  template 
  < 
    class...T, 
    typename std::enable_if<!same_underlying_type<ThisType<R>, T...>::value && 
                            !same_underlying_type<Tinit, T...>::value
			   >::type
  > 
  explicit SymmetricCartesianTensor(T&&... elems) : _data{std::forward<T>(elems)...} {};
  //! initialize from an Arr object
  explicit SymmetricCartesianTensor(Tinit data) : _data(data) {};
  //! move from an Arr object
  explicit SymmetricCartesianTensor(Arr<Elem,N>&& data) : _data(std::move(data)) {};
  //! copy constructor
  SymmetricCartesianTensor(const ThisType<R>& src) : _data(src._data) {};
  //! move constructor
  SymmetricCartesianTensor(ThisType<R>&& src) : _data(std::move(src._data)) {};
  //! access the elements
  using iterator = typename Arr<Elem,N>::iterator;
  using const_iterator = typename Arr<Elem,N>::const_iterator;
  iterator begin() { return _data.begin(); }
  const_iterator begin() const { return _data.begin(); }
  const_iterator cbegin() const { return _data.cbegin(); }
  iterator end() { return _data.end(); }
  const_iterator end() const { return _data.end(); }
  const_iterator cend() const { return _data.cend(); }
  const Elem& operator[](size_t i) const { return _data[i]; }
  Elem& operator[](size_t i) { return _data[i]; }

  //! This helper class maps the index of each element to k -- k = R - nx
  template <class, size_t> struct kMaker;
  template <size_t i> struct kMaker<size_t, i> {
    static constexpr size_t value = ceil( (sqrt(9+8*i)-3)/2 );
  };
  //! Array of k for all the elememnts -- k = R - nx
  static constexpr Arr<size_t, N> _k = NDArray<Arr, kMaker, size_t, IndexSeq<N>>::value;

  //! This helper class maps the indices of _data to the corresponding tricon array
  using Tricon = Arr<size_t, 3>;
  template<class,size_t> struct TriconMaker;
  template<size_t i> struct TriconMaker<Tricon, i> {
    static constexpr size_t k = _k[i];
    static constexpr size_t subN = binomial(k+1u,2u);
    static constexpr size_t nx = R-k;
    static constexpr size_t ny = k % 2 ? k-i+subN : i-subN;
    static constexpr size_t nz = R-nx-ny;
    static constexpr Tricon value = { nx, ny, nz }; 
  };
  //! Array of 'tricon' which are the indices of the compressed symmetric Cartesian tensor
  using Tricons = Arr<Tricon, N>;
  static constexpr Tricons _tricon = NDArray<Arr, TriconMaker, Tricon, IndexSeq<N>>::value; 
  //! This helper function maps the indices of _data to the corresponding tricon array
  static constexpr Tricon index2tricon(const size_t i) { return _tricon[i]; }
  //! This helper function maps the tricon to the indices of _data
  static constexpr size_t tricon2index(const Tricon& tricon) {
    //nx = tricon[0], ny = tricon[1], nz = tricon[2], k = R - nx
    return binomial(R-tricon[0]+1u, 2u) + ((R-tricon[0]) % 2 ? tricon[2] : tricon[1]); 
  };
  
  //! Helper class for initializing _g
  template <class, size_t> struct gMaker;
  template <size_t i> struct gMaker<size_t, i> {
    static constexpr auto tricon = _tricon[i];
    static constexpr auto nx = tricon[0];
    static constexpr auto ny = tricon[1];
    static constexpr size_t value = trinomial(R, nx, ny); 
  };
  //! The degenerancy of each element in _data with respect to the full tensor
  using Tdegen = Arr<size_t, N>;
  static constexpr Tdegen _g = NDArray<Arr, gMaker, size_t, IndexSeq<N>>::value;

  //! Helper class for initializing MatricesForContraction
  template <size_t i, class Type> struct MakeMatrixForContraction {
    static constexpr size_t Rrhs = i;
    static constexpr size_t Rans = R - Rrhs;
    static constexpr size_t Nrhs = binomial(Rrhs+2u,2u);
    static constexpr size_t Nans = binomial(Rans+2u,2u);
    //! Helper class for matrix initialization
    template<class T, size_t irhs, size_t ians> struct Fn {
      Tinit __arr;
      static constexpr auto triconans = ThisType_igtr<Rans>::index2tricon(ians);
      static constexpr auto triconrhs = ThisType_igtr<Rrhs>::index2tricon(irhs);
      static constexpr Tricon tricon = {triconans[0]+triconrhs[0],triconans[1]+triconrhs[1],triconans[2]+triconrhs[2]};
      static constexpr size_t index = tricon2index(tricon); 
      Fn(Tinit arr) : __arr(arr) {};
      const T value = __arr[index]; 
    };
    Tinit _arr;
    MakeMatrixForContraction(Tinit arr) : _arr(arr) {};
    using Twrapper = NDArray<Arr, Fn, Type, IndexSeq<Nans,Nrhs>, Tinit>;
    //! The matrix type is Arr<Arr<Type, Nrhs>,Nans>
    using type = typename Twrapper::type; 
    const type value = Twrapper{_arr}.value;
  };
/*   //! NOTE: I don't think it's often the case that we need to perform contraction of all ranks 
 *   //so initializing this matrices in construction is not neccessary -- if you need it you can 
 *   //uncommented this out
 *   template <size_t i> using MakeMatrixForContraction_refwrap = MakeMatrixForContraction<i, Trefelem>;
 *   using Tmaker_contraction = Tuple_Maker<std::tuple, MakeMatrixForContraction_refwrap, container_index<R+1>, Tinit>;
 *   //! A tuple of matrices with the elements being the reference to _data. 
 *   //The i'th matrix has C(R-i+2,2) rows and C(i+2,2) columns so that
 *   //a tensor of rank i, expressed in a an array of C(i+2,2), multipling
 *   //this matrix will result in a tensor of rank R-i
 *   const typename Tmaker_contraction::type MatricesForContraction = Tmaker_contraction{_data}.value;
 */

  
  //! contraction with an array (target operator)
  template < size_t Nrhs >
  typename std::enable_if<
    (N >= Nrhs) && 
    isTensorSize<Nrhs>::value,
    ThisType_tr<R-isTensorSize<Nrhs>::Rank, (R - isTensorSize<Nrhs>::Rank >= 2) && Detraced >
  >::type
  operator*(const Arr<Elem,Nrhs>& rhs) const {
    static constexpr size_t Rrhs = isTensorSize<Nrhs>::Rank;
    static constexpr size_t Rans = R - Rrhs;
    static constexpr bool ifDetraced = (Rans >= 2) && Detraced;
    const auto M = MakeMatrixForContraction<Rrhs,Elem>{_data}.value; 
    const auto v = rhs & ThisType_igtr<Rrhs>::_g; //just element-wise product, resulting in a vector
    //TODO: we could optimize this matrix-vector product
    return ThisType_tr<Rans, ifDetraced>{M * v};
  };

  //! contraction with another object whose rank is smaller than this one's (delegating operator)
  template < size_t Rrhs, bool ifDetraced >
  auto operator*(const ThisType_tr<Rrhs,ifDetraced>& rhs) const ->
  typename std::enable_if<
    R >= Rrhs,
    decltype((*this)*(rhs._data))
  >::type
  const {
    return this->operator*(rhs._data);
  };

  //! contraction with another object whose rank is larger than this one's (delegating operator)
  template < size_t Rrhs, bool ifDetraced >
  auto operator*(const ThisType_tr<Rrhs, ifDetraced>& rhs) const ->
  typename std::enable_if<
    (R < Rrhs),
    decltype(rhs*_data)
  >::type
  const {
    return rhs.operator*(_data);
  };
 
  //! Maximal rank of trace
  static constexpr size_t Maxlm = Detraced ? 0 : ( R % 2 ? (R-1)/2 : R/2 );
  //! Min lc
  static constexpr size_t Minlc = 0;
  //! Max lc
  static constexpr size_t Maxlc = R;
  //! helper class to initialize MatricesForTrace
  template <size_t lm, class Type> struct MakeMatrixForTrace {
    static constexpr size_t Rans = R - 2*lm;
    using Tans = ThisType<Rans>;
    //we are using lm not 2*lm because the number of non-zero elments 
    //in a 2*lm delta tensor is the number of elements in a tensor of rank lm
    using Tdelta = ThisType<lm>;
    static constexpr size_t Nans = Tans::N;
    //Ndelta is basically how many non-zero elments in the delta tensor of rank 2*lm
    static constexpr size_t Ndelta = Tdelta::N;
    Tinit _arr;
    MakeMatrixForTrace(Tinit arr) : _arr(arr) {};
    template < class T, size_t idelta, size_t ians > struct Fn {
      static constexpr auto tricondelta = Tdelta::index2tricon(idelta);
      static constexpr auto triconans = Tans::index2tricon(ians);
      //This is just equation 2.4 in Applequist, J., J. of Phys. A: Math. Gen., 1989, 22, 4303  
      static constexpr Tricon tricon = {triconans[0]+2*tricondelta[0],triconans[1]+2*tricondelta[1],triconans[2]+2*tricondelta[2]};
      static constexpr size_t index = tricon2index(tricon); 
      Tinit __arr;
      Fn(Tinit arr) : __arr(arr) {};
      const T value = __arr[index];
    };
    using Twrapper = NDArray<Arr, Fn, Type, IndexSeq<Nans,Ndelta>, Tinit>;
    using type = typename Twrapper::type; 
    const type value = Twrapper{_arr}.value;
  };

/*   //! NOTE: I don't think it's often the case that we need to perform all the traces 
 *   //so initializing this matrices in construction is not neccessary -- if you need it you can 
 *   //uncommented this out
 *   template<size_t lm> using MakeMatrixForTrace_refwrap = MakeMatrixForTrace<lm,Trefelem>;
 *   using Tmaker_trace = Tuple_Maker<std::tuple, MakeMatrixForTrace_refwrap, container_index<Maxlm+1>, Tinit>;
 *   //! A tuple of matrices for performing traces
 *   //The i'th matrix has C(R-2*lm+2,2) rows and C(lm+2,2) columns so 
 *   //that the sum over the elements in one row gives the corresponding 
 *   //component the resultant tensor -- which is equivalent to take a 
 *   //dot product between that row and a delta tensor of rank 2*lm
 *   const typename Tmaker_trace::type MatricesForTrace = Tmaker_trace{_data}.value;
 */
  
  //! lm-fold trace
  template < size_t lm >
  typename std::enable_if<lm <= Maxlm, ThisType<R-2*lm>>::type
  trace() const {
    // just a matrix-vector product
    return ThisType<R-2*lm>{MakeMatrixForTrace<lm,Elem>{_data}.value * ThisType<lm>::_g};
  }

  //! map the i'th tricon into a set of index of the full tensor
  template < class Tarr,  size_t i > struct tricon2indexset {
    static constexpr auto& tricon = ThisType<R>::_tricon[i];
    template < class T, size_t ia > struct F {
      static constexpr size_t value = ia < tricon[0] ? 0 : ( ia < tricon[0]+tricon[1] ? 1 : 2 ); 
    };
    using Tmaker = NDArray<Arr, F, size_t, IndexSeq<R>>;
    static constexpr Tarr value = Tmaker::value;
  };
  using Tmaker_indexset = NDArray<Arr, tricon2indexset, Arr<size_t, R>, IndexSeq<N>>;
  using Tindexset = typename Tmaker_indexset::type;
  static constexpr Tindexset indexset = Tmaker_indexset::value; 

  //! Given a 3x3 rotation matrix, calculate its R-fold tensor product in compressed notation
  //as described in eq. (37) in J. Math. Phys. 24 (4), April 1983
  //TODO: we can generate the index pair (ai[k],bi[k]) at compile time
  TeMatrix<N,N> expandrotmatrix(const TeMatrix<3,3>& RM) const {
    TeMatrix<N,N> ans{};
    for(size_t i = 0; i < N; ++i) {
      const auto& ai = indexset[i];
      for(size_t j = 0; j < N; ++j) {
	auto bj = indexset[j];
	do {
	  Elem cum{1};
	  for(size_t k = 0; k < R; ++k) {
	    cum *= RM[ ai[k] ][ bj[k] ];
	  }
	  ans[i][j] += cum;
	} while(std::next_permutation(bj.begin(), bj.end()));
      }
    }
    return ans;
  }

  //! rotate the tensor given a quaternion
  ThisType<R> xform(const Arr<Elem,4>& quat) const {
    const auto RM = quat_to_mat(quat);
    //GM is the matrix described in eq. (37) in J. Math. Phys. 24 (4), April 1983
    TeMatrix<N,N> GM = expandrotmatrix(RM);
    return ThisType<R>{ GM * _data };
  }

  //! Initialize MatricesForBilinearForm 
  template<size_t lc, size_t qm = 0> struct MakeMatricesForBilinearForm {
    static_assert(lc + qm <= R, "Fold of contraction exceeds tensor rank");
    //maxlm is the maximum of lm given lc
    static constexpr size_t maxlm = (Detraced || R < 2) ? 0 : ((R-qm-lc) % 2 ? (R-qm-lc-1) / 2 : (R-qm-lc)/2);
    static constexpr size_t minlm = (Detraced || R < 2) ? 0 : ( maxlm > 0 ? 1 : 0 );
    static constexpr size_t Nx = maxlm + 1;
    static constexpr size_t Ny = binomial(lc+qm+2u, 2u);
    using Twrapper = NDArray<Arr, init_to_zeros, Elem, IndexSeq<Nx, Ny>>;
    using type = typename Twrapper::type; 
    static constexpr type value = Twrapper::value;
  };
  template<size_t lc> using MakeMatricesForBilinearForm0 = MakeMatricesForBilinearForm<lc,0>;
  template<size_t lc> using MakeMatricesForBilinearForm_qm = MakeMatricesForBilinearForm<lc,1>;
  using Tmaker_bilinear0 = Tuple_Maker<std::tuple, MakeMatricesForBilinearForm0, container_index<Maxlc+1>>;
  //! Matrices holding EMP^(m:lm) .Om. r^(Om), where Om = m-lm-lc, for each lc
  //(first dimension) and lm (second dimension) -- these are used in bilinear
  //form in EMP interaction
  using TMatricesForBilinearForm0 = typename Tmaker_bilinear0::type;

  using Tmaker_bilinear_qm = Tuple_Maker<std::tuple, MakeMatricesForBilinearForm_qm, container_index<Maxlc>>;
  using TMatricesForBilinearForm_qm = typename Tmaker_bilinear_qm::type;

  /**
  * @brief Directly populate the matrix for bilinear form when lc = [maxlc_w_lm-1, maxlc_w_lm] && lm = [1, maxlm]
  *
  * @tparam lc lc 
  * @tparam lm lm
  * @tparam maxlc_w_lm the largest lc where lm >= 1 is needed
  * @param Matrices matrix
  *
  * @return matrix
  */
  template <size_t lc, size_t lm, size_t maxlc_w_lm>
  typename std::enable_if<
    lm >= 1 && lm <= MakeMatricesForBilinearForm<lc>::maxlm &&
    lc >= maxlc_w_lm-1 && lc <= maxlc_w_lm,
    TMatricesForBilinearForm0&
  >::type
  GenMatricesForBilinearForm_lm(TMatricesForBilinearForm0& Matrices) const {
    static constexpr size_t RankSource = lc + 2*lm;
    const auto& v = std::get<RankSource>(Matrices)[0];
    using TMake_matrix = typename ThisType<RankSource>::template MakeMatrixForTrace<lm,Elem>;
    std::get<lc>(Matrices)[lm] = TMake_matrix{v}.value * ThisType_igtr<lm>::_g;
    return GenMatricesForBilinearForm_lm<lc, lm+1, maxlc_w_lm>(Matrices);
  }

  /**
  * @brief Recursively populate the matrix for bilinear form when lc = [0, maxlc_w_lm-1) && lm = [1, maxlm]
  *
  * @tparam lc lc 
  * @tparam lm lm
  * @tparam maxlc_w_lm the largest lc where lm >= 1 is needed
  * @param Matrices matrix
  *
  * @return matrix
  */
  template <size_t lc, size_t lm, size_t maxlc_w_lm>
  typename std::enable_if<
    lm >= 1 && lm <= MakeMatricesForBilinearForm<lc>::maxlm &&
    lc < maxlc_w_lm - 1,
    TMatricesForBilinearForm0&
  >::type
  GenMatricesForBilinearForm_lm(TMatricesForBilinearForm0& Matrices) const {
    //TODO: we can parallelize for all lm in this recursion since 
    //all the traces are taken from tensor of the same rank (lc+2)
    const auto& v = std::get<lc+2>(Matrices)[lm-1];
    //!Since v is itself a symmetric cartesian tensor, we can build a matrix for trace from it
    //and then we do the matrix-vector product to get the trace
    using TMake_matrix = typename ThisType<lc+2>::template MakeMatrixForTrace<1,Elem>;
    //! this matrix-vector product is just the trace
    std::get<lc>(Matrices)[lm] = TMake_matrix{v}.value * ThisType_igtr<1>::_g;
    return GenMatricesForBilinearForm_lm<lc, lm+1, maxlc_w_lm>(Matrices);
  }

  /**
  * @brief termination of the recursion in the case where lc > maxlc_w_lm or maxlm == 0
  * @details lc > maxlc_w_lm will be triggered when maxlc_w_lm < Maxlc; maxlm == 0 will
  * be triggered if Detraced == true or lc >= Maxlc - 1; lm > maxlm will be triggered
  * in the upward recursion about lm
  * @tparam lc lc 
  * @tparam lm lm
  * @tparam maxlc_w_lm the largest lc where lm >= 1 is needed
  * @param Matrices matrix
  *
  * @return matrix
  */
  template <size_t lc, size_t lm, size_t maxlc_w_lm>
  typename std::enable_if<
    (lm == 0 && MakeMatricesForBilinearForm<lc>::maxlm == 0) ||
    (lm > MakeMatricesForBilinearForm<lc>::maxlm) ||
    (lc > maxlc_w_lm),
    TMatricesForBilinearForm0&
  >::type
  GenMatricesForBilinearForm_lm(TMatricesForBilinearForm0& Matrices) const {
    return Matrices;
  }

  /**
  * @brief Termination of the recursion when lc == Maxlc
  *
  * @tparam lc lc
  * @tparam maxlc_w_lm the largest lc where lm >= 1 is needed
  * @param r distance vector
  * @param Matrices matrix
  *
  * @return matrix
  */
  template <size_t lc, size_t maxlc_w_lm>
  typename std::enable_if<
    lc == Maxlc,
    TMatricesForBilinearForm0&
  >::type
  GenMatricesForBilinearForm(const Arr<Elem,3>& r, TMatricesForBilinearForm0& Matrices) const {
    std::get<lc>(Matrices)[0] = _data;
    return Matrices;
  }

  /**
  * @brief Recursively populate the matrix for bilinear form (propagation)
  * @details The recursion is two dimensional. 1st dimension is upward recursion 
  * on lc and 2nd dimension is upward recursion about lm
  * @tparam lc lc
  * @tparam maxlc_w_lm the largest lc where lm >= 1 is needed
  * @param r distance vector
  * @param Matrices matrix
  *
  * @return matrix
  */
  template <size_t lc, size_t maxlc_w_lm>
  typename std::enable_if<
    lc < Maxlc,
    TMatricesForBilinearForm0&
  >::type
  GenMatricesForBilinearForm(const Arr<Elem,3>& r, TMatricesForBilinearForm0& Matrices) const {
    static constexpr size_t Rup = lc + 1;
    //TODO: we don't actually need to recurse upto lc == Maxlc since we might not actually need
    //all the results -- all we need is lc = [maxlc_w_lm, maxlc_w_lm + 2*MakeMatricesForBilinearForm<lc>::maxlm];
    //however, this range usually covers all lc = [0,Maxlc] for small R so we don't have to specialize 
    //it
    const auto& v = std::get<Rup>(GenMatricesForBilinearForm<Rup, maxlc_w_lm>(r,Matrices))[0];
    using TMake_matrix = typename ThisType_igtr<Rup>::template MakeMatrixForContraction<1,Elem>;
    std::get<lc>(Matrices)[0] = TMake_matrix{v}.value * r;
    GenMatricesForBilinearForm_lm<lc, MakeMatricesForBilinearForm<lc>::minlm, maxlc_w_lm>(Matrices);
    return Matrices;
  }

  /**
  * @brief Recursively populate the matrix for bilinear form (initiation)
  * @details The recursion is two dimensional. 1st dimension is upward recursion 
  * on lc and 2nd dimension is upward recursion about lm
  * @tparam lc lc
  * @tparam maxlc_w_lm the largest lc where lm >= 1 is needed
  * @param r distance vector
  * @param Matrices matrix
  *
  * @return matrix
  */
  template <size_t maxlc_w_lm = Maxlc>
  TMatricesForBilinearForm0 GenMatricesForBilinearForm(const Arr<Elem,3>& r) const {
    TMatricesForBilinearForm0 Matrices = Tmaker_bilinear0::value;
    return GenMatricesForBilinearForm<0,maxlc_w_lm>(r, Matrices);
  }

  template < size_t ... i > 
  TMatricesForBilinearForm_qm
  GenMatricesForBilinearForm_qm_impl(const TMatricesForBilinearForm0& M, const IndexSeq<i...>&) const {
    return TMatricesForBilinearForm_qm{ std::get<i+1>(M)... };
  }

  TMatricesForBilinearForm_qm
  GenMatricesForBilinearForm_qm(const TMatricesForBilinearForm0& M) const {
    static constexpr size_t N = std::tuple_size<TMatricesForBilinearForm0>::value - 1;
    return GenMatricesForBilinearForm_qm_impl(M, container_index<N>{});
  }

};

template < template<class, size_t> class Arr, class Elem, size_t R, bool Detraced >
constexpr size_t SymmetricCartesianTensor<Arr, Elem, R, Detraced>::N; 

template < template<class, size_t> class Arr, class Elem, size_t R, bool Detraced >
constexpr typename SymmetricCartesianTensor<Arr, Elem, R, Detraced>::Tricons SymmetricCartesianTensor<Arr, Elem, R, Detraced>::_tricon; 

template < template<class, size_t> class Arr, class Elem, size_t R, bool Detraced >
constexpr typename SymmetricCartesianTensor<Arr, Elem, R, Detraced>::Tdegen SymmetricCartesianTensor<Arr, Elem, R, Detraced>::_g; 

template < template<class, size_t> class Arr, class Elem, size_t R, bool Detraced >
template <size_t lc, size_t qm>
constexpr size_t SymmetricCartesianTensor<Arr, Elem, R, Detraced>::MakeMatricesForBilinearForm<lc,qm>::minlm;

template < template<class, size_t> class Arr, class Elem, size_t R, bool Detraced >
template <size_t lc, size_t qm>
constexpr size_t SymmetricCartesianTensor<Arr, Elem, R, Detraced>::MakeMatricesForBilinearForm<lc,qm>::maxlm;

template < template<class, size_t> class Arr, class Elem, size_t R, bool Detraced >
constexpr typename SymmetricCartesianTensor<Arr, Elem, R, Detraced>::Tindexset
SymmetricCartesianTensor<Arr, Elem, R, Detraced>::indexset;

/**
* @brief Handle cartesian tensor bilinear contraction  
* (result in a rank-0 tensor or a scalar)
*
* @tparam T1 Left operand symmetric cartesian tensor class
* @tparam T2 Right operand symmetric cartesian tensor class
* @tparam Tc The central tensor 
*/
template < class T1, class T2 > struct SymmetricCartesianTensorBilinearForm;

/**
* @brief Partial specialized for contraction of two symmetric 
* cartesian tensor with the gradient of a function of r
*
* @tparam Arr Array template holding the tensor components
* @tparam Elem Element type of the tensor
* @tparam R1 rank of the left-hand-side tensor
* @tparam R2 rank of the right-hand-side tensor
* @tparam Detraced1 if the left-hand-side tensor is detraced
* @tparam Detraced2 if the right-hand-side tensor is detraced
* @tparam C template for the central tensor functor
*/
template < template<class, size_t> class Arr, class Elem, size_t R1, size_t R2, bool Detraced1, bool Detraced2 >
struct SymmetricCartesianTensorBilinearForm
<
  SymmetricCartesianTensor<Arr,Elem,R1,Detraced1>, 
  SymmetricCartesianTensor<Arr,Elem,R2,Detraced2> 
>
{
  template <size_t R, bool Detraced> using Ttensor = SymmetricCartesianTensor<Arr,Elem,R,Detraced>;
  template <size_t Nrow, size_t Ncol> using TeMatrix = Arr<Arr<Elem,Ncol>,Nrow>;
  template <size_t Ra, size_t Rb, bool Detraceda, bool Detracedb> using ThisType_tr = 
  SymmetricCartesianTensorBilinearForm<
    SymmetricCartesianTensor<Arr,Elem,Ra,Detraceda>, 
    SymmetricCartesianTensor<Arr,Elem,Rb,Detracedb>
  >;
  template <size_t Ra, size_t Rb> using ThisType = ThisType_tr<Ra,Rb,Detraced1,Detraced2>;

  using T1 = Ttensor<R1,Detraced1>;
  using T2 = Ttensor<R2,Detraced2>;
  using TM1 = typename T1::TMatricesForBilinearForm0;
  using TM2 = typename T2::TMatricesForBilinearForm0;
  static constexpr size_t N1 = T1::N;
  static constexpr size_t N2 = T2::N;
  static constexpr size_t Maxlc = Min(R1,R2);
  static constexpr size_t Maxlc_qr = Maxlc;
  static constexpr size_t Maxlc_qm = R1 == 0 ? 0 : Min(R1-1,R2);
  static constexpr size_t Maxlc_qn = R2 == 0 ? 0 : Min(R1,R2-1);
  static constexpr size_t Maxlc_tr = Maxlc_qm;
  static constexpr size_t Maxlc_tn = (R2 == 0 || R1 == 0) ? 0 : Min(R1-1,R2-1);
  //NOTE that tm = 1 term is not needed for EMP torque calculation
  //due to the antisymmetry of the Levi civita symbol
  //static constexpr size_t Maxlc_tm = R2 < 2 ? 0 : min(R1-2,R2);
  static constexpr size_t Rtot = R1 + R2;

  template <size_t lc, size_t qr = 0, size_t qm = 0, size_t qn = 0> struct MakeCoefficientMatrices {
    static constexpr size_t q = qr + qm + qn;
    static constexpr size_t qmn = qm + qn;
    static_assert(q <= 1, "qr + qm + qn must <= 1");
    using MatricesMaker1 = typename T1::template MakeMatricesForBilinearForm<lc, qm>;
    using MatricesMaker2 = typename T2::template MakeMatricesForBilinearForm<lc, qn>;
    static constexpr size_t maxlm = MatricesMaker1::maxlm;
    static constexpr size_t maxln = MatricesMaker2::maxlm;
    template < class T, size_t ln, size_t lm > struct Fn {
      static constexpr size_t k = lc + lm + ln + qmn;
      static constexpr T value = 
        binomial(R1,qm) * binomial(R2,qn) * 
	(R1 == 0 ? 1 : NpairK(R1-qm,lm)) * (R2 == 0 ? 1 : NpairK(R2-qn,ln)) * 
	(R1 == 0 ? 1 : binomial(R1-qm-2*lm,lc)) * (R2 == 0 ? 1 : binomial(R2-qn-2*ln,lc)) * multifac(lc,1) * 
	(size_t{1} << (Rtot+q-k)) ; 
    };
    using Twrapper = NDArray<Arr, Fn, size_t, IndexSeq<maxlm+1, maxln+1>>;
    using type = typename Twrapper::type;
    static constexpr type value = Twrapper::value;
  };
  template <size_t lc> using MakeCoefficientMatrices0 = MakeCoefficientMatrices<lc,0,0,0>;
  template <size_t lc> using MakeCoefficientMatrices_qr = MakeCoefficientMatrices<lc,1,0,0>;
  template <size_t lc> using MakeCoefficientMatrices_qm = MakeCoefficientMatrices<lc,0,1,0>;
  template <size_t lc> using MakeCoefficientMatrices_qn = MakeCoefficientMatrices<lc,0,0,1>;

  using ThisType_for_T = ThisType_tr<R1==0 ? 0 : R1-1, R2, Detraced1 && R1 >= 3, Detraced2>;
  template <size_t lc> using MakeCoefficientMatrices_tr = typename ThisType_for_T::template MakeCoefficientMatrices_qr<lc>;
  template <size_t lc> using MakeCoefficientMatrices_tn = typename ThisType_for_T::template MakeCoefficientMatrices_qn<lc>;
  using Tmaker_coeff0 = Tuple_Maker<std::tuple, MakeCoefficientMatrices0, container_index<Maxlc+1>>;
  using Tmaker_coeff_qr = Tuple_Maker<std::tuple, MakeCoefficientMatrices_qr, container_index<Maxlc_qr+1>>;
  using Tmaker_coeff_qm = Tuple_Maker<std::tuple, MakeCoefficientMatrices_qm, container_index<R1 == 0 ? 0 : Maxlc_qm+1>>;
  using Tmaker_coeff_qn = Tuple_Maker<std::tuple, MakeCoefficientMatrices_qn, container_index<R2 == 0 ? 0 : Maxlc_qn+1>>;
  using Tmaker_coeff_tr = typename ThisType_for_T::Tmaker_coeff_qr; 
  using Tmaker_coeff_tn = typename ThisType_for_T::Tmaker_coeff_qn; 

  using CoefficientMatricesType0 = typename Tmaker_coeff0::type;
  static constexpr CoefficientMatricesType0 CoefficientMatrices = Tmaker_coeff0::value;

  using CoefficientMatricesType_qr = typename Tmaker_coeff_qr::type;
  static constexpr CoefficientMatricesType_qr CoefficientMatrices_qr = Tmaker_coeff_qr::value;

  using CoefficientMatricesType_qm = typename Tmaker_coeff_qm::type;
  static constexpr CoefficientMatricesType_qm CoefficientMatrices_qm = Tmaker_coeff_qm::value;

  using CoefficientMatricesType_qn = typename Tmaker_coeff_qn::type;
  static constexpr CoefficientMatricesType_qn CoefficientMatrices_qn = Tmaker_coeff_qn::value;

  using CoefficientMatricesType_tr = typename Tmaker_coeff_tr::type;
  static constexpr CoefficientMatricesType_tr CoefficientMatrices_tr = ThisType_for_T::CoefficientMatrices_qr;

  using CoefficientMatricesType_tn = typename Tmaker_coeff_tn::type;
  static constexpr CoefficientMatricesType_tn CoefficientMatrices_tn = ThisType_for_T::CoefficientMatrices_qn;

  //!Multiply each of the coefficients with the appropriate F component
  template < size_t qr = 0, size_t qm = 0, size_t qn = 0 >
  struct MultiplyCoefficientMatrixWithDeriv {
    static constexpr size_t q = qr + qm + qn;
    static constexpr size_t qmn = qm + qn;
    static constexpr size_t maxgrad = Rtot+q;
    MultiplyCoefficientMatrixWithDeriv(const Arr<Elem,maxgrad+1>& F) : _F(F) {}; 
    const Arr<Elem,maxgrad+1>& _F;
    //note that the input matrix is matrix of size_t while the output is matrix of Elem 
    template<size_t Nrow, size_t Ncol> using TiMatrix = Arr<Arr<size_t,Ncol>,Nrow>;
    template<size_t Nrow, size_t Ncol> TeMatrix<Nrow, Ncol> operator()(size_t lc, const TiMatrix<Nrow,Ncol>& m) {
      TeMatrix<Nrow, Ncol> ans;
      for(size_t lm = 0; lm < Nrow; ++lm) {
        for(size_t ln = 0; ln < Ncol; ++ln) {
    	  ans[lm][ln] = m[lm][ln] * _F[maxgrad-qmn-lm-ln-lc];
        }
      }
      return ans;
    }
  };

  //!Given two matrices, each representing an array of tensor, 
  //perform element-wise contraction and 
  //then sum the results together
  template <size_t qr = 0, size_t qm = 0, size_t qn = 0>
  struct ContractTwoTensorArrays {
    /**
    * @brief Contracted two arrays of tensors 
    * @details The rank of the tensors in the first array must be qm 
    * higher or qn lower than the rank of those in the second array
    * TODO: This operator should be avoided -- we can simply overload the array 
    * multiplication operator for the tensor class and make the matrix for bilinear form
    * in the tensor class array of tensors -- that way, simply M1*M2 should be enough 
    * to give the result
    *
    * @tparam Nrow number of tensors in each array
    * @tparam Ncol1 number of elements for tensors in the 1st array
    * @tparam Ncol2 number of elements for tensors in the 2nd array
    * @param M1 first array
    * @param M2 second array
    *
    * @return sum of the element-wise contraction between the two arrays
    */
    template < size_t Nrow, size_t Ncol1, size_t Ncol2>
    auto operator()(const TeMatrix<Nrow, Ncol1>& M1, const TeMatrix<Nrow, Ncol2>& M2) ->
    typename std::enable_if<
      T1::template isTensorSize<Ncol1>::value &&
      T2::template isTensorSize<Ncol2>::value &&
      ((T1::template isTensorSize<Ncol1>::Rank == T2::template isTensorSize<Ncol2>::Rank && qm + qn == 0) || 
       (T1::template isTensorSize<Ncol1>::Rank != T2::template isTensorSize<Ncol2>::Rank && qm + qn == 1)) &&
      Nrow == MakeCoefficientMatrices<T1::template isTensorSize<Ncol1>::Rank-qm, qr, qm, qn>::maxlm + 1,
      decltype( (Ttensor<T1::template isTensorSize<Ncol1>::Rank,Detraced1 && T1::template isTensorSize<Ncol1>::Rank >= 2>{M1[0]} * 
                 Ttensor<T2::template isTensorSize<Ncol2>::Rank,Detraced2 && T2::template isTensorSize<Ncol2>::Rank >= 2>{M2[0]})._data )
    >::type
    {
      static constexpr size_t rank1 = T1::template isTensorSize<Ncol1>::Rank;
      static constexpr size_t rank2 = T1::template isTensorSize<Ncol2>::Rank;
      using Tt1 = Ttensor<rank1,Detraced1 && rank1 >= 2>;
      using Tt2 = Ttensor<rank2,Detraced1 && rank2 >= 2>;
      using Tans = decltype( (Tt1{M1[0]} * Tt2{M2[0]})._data );
      Tans ans{};
      for(size_t i = 0; i < Nrow; ++i) {
	ans += (Tt1{M1[i]} * Tt2{M2[i]})._data;
      }
      return ans;
    }

    /**
    * @brief Contracted two arrays of tensors and then contract the resultant tensor with the Levi Civita symbol 
    * @details The rank of the tensors in the first array must be the same as
    * the rank of those in the second array. TODO: The underlying operation is 
    * such that the two tensor remain 1 rank after the contraction and then the 
    * two resulting vector is passed to a cross product -- this may be implemented 
    * as member operator of the tensor class
    *
    * @tparam Nrow number of tensors in each array
    * @tparam Ncol1 number of elements for tensors in the 1st array
    * @tparam Ncol2 number of elements for tensors in the 2nd array
    * @param M1 first array
    * @param M2 second array
    *
    * @return sum of the element-wise contraction between the two arrays resulting 
    * in a rank-2 tensor, which is then contracted with the Levi Civita symbol
    * 
    */
    template < size_t Nrow, size_t Ncol1, size_t Ncol2>
    auto operator()(const TeMatrix<Nrow, Ncol1>& M1, const TeMatrix<Nrow, Ncol2>& M2) ->
    typename std::enable_if<
      T1::template isTensorSize<Ncol1>::value &&
      T2::template isTensorSize<Ncol2>::value &&
      T1::template isTensorSize<Ncol1>::Rank != 0 && 
      T2::template isTensorSize<Ncol2>::Rank != 0 &&
      T1::template isTensorSize<Ncol1>::Rank == T2::template isTensorSize<Ncol2>::Rank &&
      qr == 0 && qm == 1 && qn == 0 &&
      Nrow == MakeCoefficientMatrices<T1::template isTensorSize<Ncol1>::Rank-qm, qr, qm, qn>::maxlm + 1,
      decltype( std::declval<Ttensor<1,false>>()._data )
    >::type
    {
      static constexpr size_t rank1 = T1::template isTensorSize<Ncol1>::Rank;
      static constexpr size_t rank2 = T1::template isTensorSize<Ncol2>::Rank;
      static constexpr size_t lc = rank1 - 1;
      using Tt1 = Ttensor<rank1,Detraced1 && rank1 >= 2>;
      using Tt2 = Ttensor<rank2,Detraced1 && rank2 >= 2>;
      using Tlc = Ttensor<lc,false>;
      using Tans = decltype( std::declval<Ttensor<1,false>>()._data );
      Tans ans{};
      for(size_t i = 0; i < Nrow; ++i) {
        const auto v1 = typename Tt1::template MakeMatrixForContraction<lc,Elem>{M1[i]}.value;
        const auto v2 = typename Tt2::template MakeMatrixForContraction<lc,Elem>{M2[i]}.value;
	ans[0] += (Tlc{v2[1]} * Tlc{v1[2]})._data[0] - (Tlc{v2[2]} * Tlc{v1[1]})._data[0];
	ans[1] += (Tlc{v2[2]} * Tlc{v1[0]})._data[0] - (Tlc{v2[0]} * Tlc{v1[2]})._data[0];
	ans[2] += (Tlc{v2[0]} * Tlc{v1[1]})._data[0] - (Tlc{v2[1]} * Tlc{v1[0]})._data[0];
      }
      return ans;
    }
  };

  /**
  * @brief Evaluate the bilinear form (user interface)
  * @details The central tensor is the same as 
  * Rtot, which results in a rank-0 tensor
  * @tparam F Function of r that returns 0 to Rtot fold gradient
  * @param t1 left-hand operand
  * @param t2 right-hand operand
  * @param r distance vector
  *
  * @return The image of the bilinear mapping
  */
  template <Arr<Elem,Rtot+1> F(const Elem&)> 
  Elem BF0(T1& t1, T2& t2, const Arr<Elem,3>& r) const {
    //! populate the matrices for bilinear form
    //NOTE that maxlc_w_lm should be Maxlc+2 because the traces at lc == maxlc_w_lm
    //has to come from lc = maxlc_w_lm + 2, where the traces are generated directly
    const auto M1 = t1.template GenMatricesForBilinearForm<Maxlc+2>(r);
    const auto M2 = t2.template GenMatricesForBilinearForm<Maxlc+2>(r);
    //! calculate the derivatives of the function
    const Elem rv = sqrt(r*r);
    const auto Fr = F(rv);
    //! perform the bilinear contraction
    const auto CWD = map_tuple(MultiplyCoefficientMatrixWithDeriv<0,0,0>{Fr}, CoefficientMatrices); 
    const auto CB = multiply_2tuples(CWD, M2);
    const auto ACB = map_2tuples(ContractTwoTensorArrays<0,0,0>{}, M1, CB); //ACB is an array of array<Elem,1> 

    return sumACB(ACB);
  }

  /**
  * @brief Evaluate the bilinear form (user interface)
  * @details The central tensor is 1-rank higher than 
  * Rtot, which results in a rank-1 tensor
  * @tparam F Function of r that returns 0 to Rtot fold gradient
  * @param t1 left-hand operand
  * @param t2 right-hand operand
  * @param r distance vector
  *
  * @return The image of the bilinear mapping
  */
  template <Arr<Elem,Rtot+2> F(const Elem&)> 
  Arr<Elem,3> BF1(T1& t1, T2& t2, const Arr<Elem,3>& r) const {
    //! populate the matrices for bilinear form for qr = 1
    //NOTE the maxlc_w_lm should be Maxlc + 1 + 2, where the '+1' is needed 
    //because we need to shift the index in MatricesForBilinearForm down by 1
    //in the case of qm = 1; '+2' is needed because the traces at lc == Maxlc+1
    //has to come from lc == Maxlc+1+2, where the traces are generated directly
    const auto M1_qr = t1.template GenMatricesForBilinearForm<Maxlc+1+2>(r);
    const auto M2_qr = t2.template GenMatricesForBilinearForm<Maxlc+1+2>(r);
    //! populate the matrices for bilinear form for qm = 1 or qn = 1
    const auto M1_qm = t1.GenMatricesForBilinearForm_qm(M1_qr);
    const auto M2_qn = t2.GenMatricesForBilinearForm_qm(M2_qr);
    //! calculate the derivatives of the function
    const Elem rv = sqrt(r*r);
    const auto Fr = F(rv);
    //! Multiply the coefficient matrix by the gradients
    const auto CWD_qr = map_tuple(MultiplyCoefficientMatrixWithDeriv<1,0,0>{Fr}, CoefficientMatrices_qr); 
    const auto CWD_qm = map_tuple(MultiplyCoefficientMatrixWithDeriv<0,1,0>{Fr}, CoefficientMatrices_qm); 
    const auto CWD_qn = map_tuple(MultiplyCoefficientMatrixWithDeriv<0,0,1>{Fr}, CoefficientMatrices_qn);
    //! perform the bilinear contraction for qr = 1
    const auto CB_qr = multiply_2tuples(CWD_qr, M2_qr);
    const auto ACB_qr = map_2tuples(ContractTwoTensorArrays<1,0,0>{}, M1_qr, CB_qr); 
    //! perform the bilinear contraction for qm = 1
    const auto CB_qm = multiply_2tuples(CWD_qm, M2_qr);
    const auto ACB_qm = map_2tuples(ContractTwoTensorArrays<0,1,0>{}, M1_qm, CB_qm); 
    //! perform the bilinear contraction for qn = 1
    const auto CB_qn = multiply_2tuples(CWD_qn, M2_qn);
    const auto ACB_qn = map_2tuples(ContractTwoTensorArrays<0,0,1>{}, M1_qr, CB_qn);

    return sumACB(ACB_qr, ACB_qm, ACB_qn, r);
  }

  /**
  * @brief Evaluate the bilinear form for calculating torque on the first EMP tensor (user interface)
  * @details The central tensor is the same as 
  * Rtot but contracts R1-1 with T1, which results in a rank-2 tensor
  * -- This is intended for calculating EMP torque and is not 
  *  for general purpose so the tm = 1 term is skipped since it 
  *  vanishes when contracted with the Levi Civita symbol as in 
  *  EMP torque calculation. NOTE: This will only be enable if 
  *  R1 >= 1 since monopole can't experience torque
  * @tparam F Function of r that returns 0 to Rtot fold gradient
  * @param t1 left-hand operand where the torque is exerted 
  * @param t2 right-hand operand
  * @param r distance vector
  *
  * @return The image of the bilinear mapping
  */
  template <Arr<Elem,Rtot+1> F(const Elem&)>
  typename std::conditional<R1 != 0, Arr<Elem,3>, void>::type
  BFT(T1& t1, T2& t2, const Arr<Elem,3>& r) const {
    //! populate the matrices for bilinear form  
    //NOTE the maxlc_w_lm should be Maxlc + 1 + 2, where the '+1' is needed 
    //because we need to shift the index in MatricesForBilinearForm down by 1
    //in the case of qm = 1; '+2' is needed because the traces at lc == Maxlc+1
    //has to come from lc == Maxlc+1+2, where the traces are generated directly
    const auto M1_t0 = t1.template GenMatricesForBilinearForm<Maxlc+1+2>(r);
    const auto M2_t0 = t2.template GenMatricesForBilinearForm<Maxlc+1+2>(r);
    //! populate the matrices for bilinear form for tr = 1 or tn = 1
    //NOTE that the matrices actually used are M1, M2_t0 and M2_tn
    const auto M1 = t1.GenMatricesForBilinearForm_qm(M1_t0);
    const auto M2_tn = t2.GenMatricesForBilinearForm_qm(M2_t0);
    //! calculate the derivatives of the function
    const Elem rv = sqrt(r*r);
    const auto Fr = F(rv);
    //! Multiply the coefficient matrix by the gradients
    const auto CWD_tr = map_tuple(typename ThisType_for_T::template MultiplyCoefficientMatrixWithDeriv<1,0,0>{Fr}, CoefficientMatrices_tr); 
    const auto CWD_tn = map_tuple(typename ThisType_for_T::template MultiplyCoefficientMatrixWithDeriv<0,0,1>{Fr}, CoefficientMatrices_tn);
    //! perform the bilinear contraction for tr = 1
    const auto CB_tr = multiply_2tuples(CWD_tr, M2_t0);
    const auto ACB_tr = map_2tuples(ContractTwoTensorArrays<0,1,0>{}, M1, CB_tr); 
    //! perform the bilinear contraction for tn = 1
    const auto CB_tn = multiply_2tuples(CWD_tn, M2_tn);
    const auto ACB_tn = map_2tuples(ContractTwoTensorArrays<0,1,0>{}, M1, CB_tn); 

    return sumACB_T(ACB_tr, ACB_tn, r); 
  }
  
  /**
  * @brief Sum ACB arrays from qr = 1
  *
  * @tparam Nqr number of elements in the array
  * @param ACB_qr the ACB array
  * @param r distance vector
  * 
  * @return the sum as a vector of 3
  */
  template <size_t Nqr>
  typename std::enable_if<
    Nqr == std::tuple_size<CoefficientMatricesType_qr>::value,
    Arr<Elem,3>
  >::type
  sumACB_qr(const Arr<Arr<Elem,1>,Nqr>& ACB_qr, const Arr<Elem,3>& r) const {
    return sum(ACB_qr) * r;
  }

  /**
  * @brief Sum ACB arrays for qm = 1 or qn = 1
  *
  * @tparam Nqmn number of elements in the array
  * @param ACB_qmn the ACB array
  *
  * @return the sum as a vector of 3
  */
  template <size_t Nqmn>
  typename std::enable_if<
    Nqmn == std::tuple_size<CoefficientMatricesType_qm>::value ||
    Nqmn == std::tuple_size<CoefficientMatricesType_qn>::value,
    Arr<Elem,3>
  >::type
  sumACB_qmn(const Arr<Arr<Elem,3>,Nqmn>& ACB_qmn) const {
    return sumrows(ACB_qmn);
  }

  /**
  * @brief Sum ACB arrays for qm = 1 or qn =1 in case the array are empty (R1 = 0 or R2 = 0)
  * @details The empty arrays are represented by an emtpy tuple
  *
  * @param ACB_qmn the ACB array
  *
  * @return a vector of 3 zeros
  */
  Arr<Elem,3> sumACB_qmn(const std::tuple<>& ACB_qmn) const {
    return Arr<Elem,3>{};
  }
  
  /**
  * @brief Sum the ACB arrays for qr = 1, qm = 1 and qn = 1 and sum the results
  *
  * @tparam TACB_qr type of ACB_qr
  * @tparam TACB_qm type of ACB_qm
  * @tparam TACB_qn type of ACB_qn
  * @param ACB_qr ACB for qr = 1
  * @param ACB_qm ACB for qm = 1
  * @param ACB_qn ACB for qn = 1
  * @param r the distance vector
  *
  * @return a vector of 3
  */
  template <class TACB_qr, class TACB_qm, class TACB_qn>
  Arr<Elem,3> sumACB(TACB_qr&& ACB_qr, TACB_qm&& ACB_qm, TACB_qn&& ACB_qn, const Arr<Elem,3>& r) const
  {
    return sumACB_qr(std::forward<TACB_qr>(ACB_qr),r) 
          +sumACB_qmn(std::forward<TACB_qm>(ACB_qm))
          +sumACB_qmn(std::forward<TACB_qn>(ACB_qn));
  }

  /**
  * @brief Sum the ACB array for 0-rank gradient
  *
  * @tparam N number of elements of the ACB array
  * @param ACB the ACB array
  *
  * @return the sum
  */
  template <size_t N>
  typename std::enable_if<
    N == std::tuple_size<CoefficientMatricesType0>::value,
    Elem
  >::type
  sumACB(const Arr<Arr<Elem,1>,N>& ACB) const {
    return sum(ACB);
  }

  /**
  * @brief Sum ACB arrays from tr = 1
  *
  * @tparam Ntr number of elements in the array
  * @param ACB_tr the ACB array
  * @param r distance vector
  * 
  * @return the sum as a vector of 3
  */
  template <size_t Ntr>
  typename std::enable_if<
    Ntr == std::tuple_size<CoefficientMatricesType_tr>::value,
    Arr<Elem,3>
  >::type
  sumACB_tr(const Arr<Arr<Elem,3>,Ntr>& ACB_tr, const Arr<Elem,3>& r) const {
    //In Einstein notation for tensor contraction,
    //this is actually LV .2. (sumrows(ACB_tr) .outer. r ),
    //or LV_iab sumrows(ACB_tr)_b r_a, where
    //.outer. is the tensor product, LV is the Levi-Civita symbol
    //and _iab, _b and _a are the indices of the tensor
    //We can simply view the contractions as a cross product, where 
    //we need to shift the order between ACB_tr and r
    return r ^ sumrows(ACB_tr);
  }

  /**
  * @brief Sum ACB arrays from tn = 1
  *
  * @tparam Ntn number of elements in the array
  * @param ACB_tn the ACB array
  * 
  * @return the sum as a vector of 3
  */
  template <size_t Ntn>
  typename std::enable_if<
    Ntn == std::tuple_size<CoefficientMatricesType_tn>::value,
    Arr<Elem,3>
  >::type
  sumACB_tn(const Arr<Arr<Elem,3>,Ntn>& ACB_tn) const {
    return sumrows(ACB_tn);
  }

  /**
  * @brief Sum ACB arrays from tn = 1 in case R2 == 0
  * @details The empty arrays are represented by an emtpy tuple
  *
  * @tparam Ntn number of elements in the array
  * @param ACB_tn the ACB array
  * 
  * @return the sum as a vector of 3
  */
  Arr<Elem,3> sumACB_tn(const std::tuple<>& ACB_tn) const {
    return Arr<Elem,3>{};
  }

  /**
  * @brief Sum the ACB arrays for EMP torque calculation 
  *
  * @tparam TACB_tr type of ACB_tr
  * @tparam TACB_tn type of ACB_tn
  * @param ACB_tr ACB for tr = 1
  * @param ACB_tn ACB for tn = 1
  * @param r the distance vector
  *
  * @return a vector of 3
  */
  template <class TACB_tr, class TACB_tn>
  Arr<Elem,3> sumACB_T(TACB_tr&& ACB_tr, TACB_tn&& ACB_tn, const Arr<Elem,3>& r) const
  {
    return sumACB_tr(std::forward<TACB_tr>(ACB_tr),r) 
          +sumACB_tn(std::forward<TACB_tn>(ACB_tn));
  }
};

template < template<class, size_t> class Arr, class Elem, size_t R1, size_t R2, bool Detraced1, bool Detraced2 >
constexpr typename 
SymmetricCartesianTensorBilinearForm<SymmetricCartesianTensor<Arr,Elem,R1,Detraced1>, SymmetricCartesianTensor<Arr,Elem,R2,Detraced2>>::CoefficientMatricesType0 
SymmetricCartesianTensorBilinearForm<SymmetricCartesianTensor<Arr,Elem,R1,Detraced1>, SymmetricCartesianTensor<Arr,Elem,R2,Detraced2>>::CoefficientMatrices; 

template < template<class, size_t> class Arr, class Elem, size_t R1, size_t R2, bool Detraced1, bool Detraced2 >
constexpr typename 
SymmetricCartesianTensorBilinearForm<SymmetricCartesianTensor<Arr,Elem,R1,Detraced1>, SymmetricCartesianTensor<Arr,Elem,R2,Detraced2>>::CoefficientMatricesType_qr 
SymmetricCartesianTensorBilinearForm<SymmetricCartesianTensor<Arr,Elem,R1,Detraced1>, SymmetricCartesianTensor<Arr,Elem,R2,Detraced2>>::CoefficientMatrices_qr; 

template < template<class, size_t> class Arr, class Elem, size_t R1, size_t R2, bool Detraced1, bool Detraced2 >
constexpr typename 
SymmetricCartesianTensorBilinearForm<SymmetricCartesianTensor<Arr,Elem,R1,Detraced1>, SymmetricCartesianTensor<Arr,Elem,R2,Detraced2>>::CoefficientMatricesType_qm 
SymmetricCartesianTensorBilinearForm<SymmetricCartesianTensor<Arr,Elem,R1,Detraced1>, SymmetricCartesianTensor<Arr,Elem,R2,Detraced2>>::CoefficientMatrices_qm; 

template < template<class, size_t> class Arr, class Elem, size_t R1, size_t R2, bool Detraced1, bool Detraced2 >
constexpr typename 
SymmetricCartesianTensorBilinearForm<SymmetricCartesianTensor<Arr,Elem,R1,Detraced1>, SymmetricCartesianTensor<Arr,Elem,R2,Detraced2>>::CoefficientMatricesType_qn 
SymmetricCartesianTensorBilinearForm<SymmetricCartesianTensor<Arr,Elem,R1,Detraced1>, SymmetricCartesianTensor<Arr,Elem,R2,Detraced2>>::CoefficientMatrices_qn; 

template < template<class, size_t> class Arr, class Elem, size_t R1, size_t R2, bool Detraced1, bool Detraced2 >
constexpr typename 
SymmetricCartesianTensorBilinearForm<SymmetricCartesianTensor<Arr,Elem,R1,Detraced1>, SymmetricCartesianTensor<Arr,Elem,R2,Detraced2>>::CoefficientMatricesType_tr 
SymmetricCartesianTensorBilinearForm<SymmetricCartesianTensor<Arr,Elem,R1,Detraced1>, SymmetricCartesianTensor<Arr,Elem,R2,Detraced2>>::CoefficientMatrices_tr; 

template < template<class, size_t> class Arr, class Elem, size_t R1, size_t R2, bool Detraced1, bool Detraced2 >
constexpr typename 
SymmetricCartesianTensorBilinearForm<SymmetricCartesianTensor<Arr,Elem,R1,Detraced1>, SymmetricCartesianTensor<Arr,Elem,R2,Detraced2>>::CoefficientMatricesType_tn 
SymmetricCartesianTensorBilinearForm<SymmetricCartesianTensor<Arr,Elem,R1,Detraced1>, SymmetricCartesianTensor<Arr,Elem,R2,Detraced2>>::CoefficientMatrices_tn; 
#endif   /* ----- #ifndef SymmetricCartesianTensor_INC  ----- */
