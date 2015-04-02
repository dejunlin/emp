#ifndef  Tuple_INC
#define  Tuple_INC

/*
 * =====================================================================================
 *
 *       Filename:  Tuple.hpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  03/01/2015 06:48:01 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */
#include <type_traits>
#include <tuple>
#include "indexseq_generator.hpp"
#include "minmax.hpp"
#include "Typetraits.hpp"

/**
* @brief Tuple maker
* @details See the specializations for more details
*
* @tparam Tuple tuple class
* @tparam Fn map a index sequence to a sequence of types and corresponding values
*/
template 
< 
  template <class...> class Tuple, 
  template <size_t> class Fn,
  class ... T
>
struct Tuple_Maker;

/**
* @brief Make a tuple out of a index sequence
* @details By expanding the input index sequence and map them to a set of types and values
*
* @tparam Tuple tuple class
* @tparam Fn map a index sequence to a sequence of types and corresponding values
* @tparam I index sequence wrapper
* @tparam i... indices
*/
template 
< 
  template <class...> class Tuple, 
  template <size_t> class Fn,
  template <size_t...> class I,
  size_t ... i 
>
struct Tuple_Maker<Tuple, Fn, I<i...>> {
  using type = Tuple<typename Fn<i>::type...>;
  static_assert(std::is_literal_type<type>::value, "Types inside the tuple must be literal");
  static constexpr type value{Fn<i>::value...};
};

/**
* @brief Partial specialization for non literal type
* @details The idea is similar to the specialization for literal type except that
* we need an initializer to do the math at run time
*
* @tparam Tuple tuple class
* @tparam Fn map a index sequence to a sequence of types and corresponding values
* @tparam I index sequence wrapper
* @tparam i... indices
* @tparam Tinit Initializer class for Fn
*/
template 
< 
  template <class...> class Tuple, 
  template <size_t> class Fn,
  template <size_t...> class I,
  size_t ... i,
  class Tinit 
>
struct Tuple_Maker<Tuple, Fn, I<i...>, Tinit> {
  using type = Tuple<typename Fn<i>::type...>;
  Tinit _init;
  Tuple_Maker(const Tinit& init) : _init(init) {}; 
  const type value{(Fn<i>){_init}.value...};
};

template < class Tuple, class Ret, template<size_t...> class I, size_t i>
Ret sum_tuple_impl(Tuple&& t, const Ret&, const I<i>&) {
  return std::get<i>(std::forward<Tuple>(t)); 
}

template < class Tuple, class Ret, template<size_t...> class I, size_t i1, size_t ... i>
typename std::enable_if
<
  sizeof...(i) != 0,
  Ret
>::type
sum_tuple_impl(Tuple&& t, const Ret& ret, const I<i1, i...>&) {
  return std::get<i1>(std::forward<Tuple>(t)) + sum_tuple_impl(t, ret, I<i...>{}); 
}

template < class ... T >
typename std::common_type<T...>::type
sum_tuple(const std::tuple<T...>& t) {
  using Ret = typename std::common_type<T...>::type;
  return sum_tuple_impl(t, Ret{}, container_index<sizeof...(T)>{});
}

/**
* @brief Helper for map_tuple 
* @details F is a functor taking the index of the tuple element as the 1st argument and
* returns the same type for each element of the input tuples
*
* @tparam F Functor type
* @tparam Tuple tuple type
* @tparam I Index sequence wrapper
* @tparam i... index for the tuples
* @param f Functor
* @param t Tuple 
* @param I index sequence
*
* @return a std::array of the returns of F
*/
template < class F, class Tuple, template <size_t...> class I, size_t ... i, class ... Args >
auto map_tuple_impl(F&& f, Tuple&& t, const I<i...>& Dummy, Args&&... args) ->
typename std::enable_if
< 
  check_if<bare<Tuple>, can_be_get>::value &&
  check_if<bare<F>, is_functor_with_args<size_t, decltype(std::get<0>(std::forward<Tuple>(t))), Args...>>::value && 
  are_same<decltype(f(i,std::get<i>(std::forward<Tuple>(t)),std::forward<Args>(args)...))...>::value,
  std::array<decltype(f(0,std::get<0>(std::forward<Tuple>(t)),std::forward<Args>(args)...)),sizeof...(i)>
>::type 
{
  using T = typename std::common_type<decltype(f(i,std::get<i>(std::forward<Tuple>(t)),std::forward<Args>(args)...))...>::type;
  constexpr size_t N = sizeof...(i);
  return std::array<T,N>{ f(i,std::get<i>(std::forward<Tuple>(t)),std::forward<Args>(args)...)... };
}

/**
* @brief Helper for map_tuple 
* @details F is a functor taking a tuple element as the 1st argument and
* returns the same type for each element of the input tuples
*
* @tparam F Functor type
* @tparam Tuple tuple type
* @tparam I Index sequence wrapper
* @tparam i... index for the tuples
* @param f Functor
* @param t Tuple 
* @param I index sequence
*
* @return a std::array of the returns of F
*/
template < class F, class Tuple, template <size_t...> class I, size_t ... i, class ... Args >
auto map_tuple_impl(F&& f, Tuple&& t, const I<i...>& Dummy, Args&&... args) ->
typename std::enable_if
< 
  check_if<bare<Tuple>, can_be_get>::value &&
  check_if<bare<F>, is_functor_with_args<decltype(std::get<0>(std::forward<Tuple>(t))), Args...>>::value && 
  are_same<decltype(f(std::get<i>(std::forward<Tuple>(t)),std::forward<Args>(args)...))...>::value,
  std::array<decltype(f(std::get<0>(std::forward<Tuple>(t)),std::forward<Args>(args)...)),sizeof...(i)>
>::type 
{
  using T = typename std::common_type<decltype(f(std::get<i>(std::forward<Tuple>(t)),std::forward<Args>(args)...))...>::type;
  constexpr size_t N = sizeof...(i);
  return std::array<T,N>{ f(std::get<i>(std::forward<Tuple>(t)),std::forward<Args>(args)...)... };
}

/**
* @brief Helper for map_tuple
* @details F is a functor taking the index of the tuple element as the 1st argument and
* returns the different types for each element of the input tuples
*
* @tparam F Functor type
* @tparam Tuple tuple type
* @tparam I Index sequence wrapper
* @tparam i... index for the tuples
* @param f Functor
* @param t Tuple 
* @param I index sequence
*
* @return a new tuple with each elment being the return of F
*/
template < class F, class Tuple, template <size_t...> class I, size_t ... i, class ... Args >
auto map_tuple_impl(F&& f, Tuple&& t, const I<i...>& Dummy, Args&&... args) ->
typename std::enable_if
< 
  check_if<bare<Tuple>, can_be_get>::value &&
  check_if<bare<F>, is_functor_with_args<size_t, decltype(std::get<0>(std::forward<Tuple>(t))), Args...>>::value && 
  !are_same<decltype(f(i,std::get<i>(std::forward<Tuple>(t)),std::forward<Args>(args)...))...>::value,
  decltype(std::make_tuple( f(i,std::get<i>(std::forward<Tuple>(t)),std::forward<Args>(args)...)...))
>::type 
{
  return std::make_tuple( f(i,std::get<i>(std::forward<Tuple>(t)),std::forward<Args>(args)...)... );
}

/**
* @brief Helper for map_tuple
* @details F is a functor taking a tuple element as the 1st argument and
* returns the different types for each element of the input tuples
*
* @tparam F Functor type
* @tparam Tuple tuple type
* @tparam I Index sequence wrapper
* @tparam i... index for the tuples
* @param f Functor
* @param t Tuple 
* @param I index sequence
*
* @return a new tuple with each elment being the return of F
*/
template < class F, class Tuple, template <size_t...> class I, size_t ... i, class ... Args >
auto map_tuple_impl(F&& f, Tuple&& t, const I<i...>& Dummy, Args&&... args) ->
typename std::enable_if
< 
  check_if<bare<Tuple>, can_be_get>::value &&
  check_if<bare<F>, is_functor_with_args<decltype(std::get<0>(std::forward<Tuple>(t))), Args...>>::value && 
  !are_same<decltype(f(std::get<i>(std::forward<Tuple>(t)),std::forward<Args>(args)...))...>::value,
  decltype(std::make_tuple( f(std::get<i>(std::forward<Tuple>(t)),std::forward<Args>(args)...)...))
>::type 
{
  return std::make_tuple( f(std::get<i>(std::forward<Tuple>(t)),std::forward<Args>(args)...)... );
}

/**
* @brief Helper for map_tuple
* @details input tuple is an empty tuple, i.e. std::tuple_size<Tuple>::value == 0 
*
* @tparam F Functor type
* @tparam Tuple tuple type
* @tparam I Index sequence wrapper
* @tparam i... index for the tuples
* @param f Functor
* @param t Tuple 
* @param I index sequence
*
* @return a empty tuple 
*/
template < class F, class Tuple, template <size_t...> class I, class ... Args >
std::tuple<> map_tuple_impl(F&& f, Tuple&& t, const I<>& Dummy, Args&&... args)
{
  return std::tuple<>{};
}

/**
* @brief map the elements of a tuple to another tuple
*
* @tparam F Functor type that performs the mapping
* @tparam Tuple Tuple type
* @tparam Args... Other parameters for F
* @param f functor
* @param t Tuple 
* @param args... Other parameters for F
*
* @return a new tuple with each elment being the return of F
*/
template < class F, class Tuple, class ... Args >
auto map_tuple(F&& f, Tuple&& t, Args&&... args) -> 
decltype(map_tuple_impl(std::forward<F>(f), std::forward<Tuple>(t), container_index<std::tuple_size<bare<Tuple>>::value>{}, std::forward<Args>(args)...))
{
  constexpr size_t N = std::tuple_size<bare<Tuple>>::value;
  return map_tuple_impl(std::forward<F>(f), std::forward<Tuple>(t), container_index<N>{});
};


/**
* @brief Just an list-initializable empty class
*/
template < class ... T > struct Swallow { Swallow(T&& ... args) {}; };

/**
* @brief Helper for foreach_2tuples
*
* @tparam F Functor type
* @tparam Tuple1 tuple 1 type
* @tparam Tuple2 tuple 2 type
* @tparam I Index sequence wrapper
* @tparam i... index for the tuples
* @param f Functor
* @param t1 Tuple 1
* @param t2 Tuple 2
* @param I index sequence
*
* @return void
*/
template < class F, class Tuple1, class Tuple2, template <size_t...> class I, size_t ... i >
typename std::enable_if
< 
  check_if<bare<Tuple1>, can_be_get>::value && 
  check_if<bare<Tuple2>, can_be_get>::value 
>::type 
foreach_2tuples_impl(F&& f, Tuple1&& t1, Tuple2&& t2, const I<i...>&) {
  Swallow
  <
    decltype
    (
      f
      (
        typename std::tuple_element<i, bare<Tuple1>>::type{}, 
	typename std::tuple_element<i, bare<Tuple2>>::type{}
      )
    )...
  >
  { f(std::get<i>(std::forward<Tuple1>(t1)),std::get<i>(std::forward<Tuple2>(t2)))... };
}

/**
* @brief Iterate through two tuple of the same sizes and apply a functor on each pair of 
* underlying elements in the order where they are stored in the input tuples
*
* @tparam F Functor type that performs the mapping
* @tparam Tuple1 Tuple 1 type
* @tparam Tuple2 Tuple 2 type 
* @param f functor
* @param t1 Tuple 1
* @param t2 Tuple 2
*
* @return void 
*/
template < class F, class Tuple1, class Tuple2 >
typename std::enable_if
< 
  std::tuple_size<bare<Tuple1>>::value == std::tuple_size<bare<Tuple2>>::value
>::type 
foreach_2tuples(F&& f, Tuple1&& t1, Tuple2&& t2) {
  constexpr size_t N = std::tuple_size<bare<Tuple1>>::value;
  foreach_2tuples_impl(std::forward<F>(f), std::forward<Tuple1>(t1), std::forward<Tuple2>(t2), container_index<N>{});
};

/**
* @brief Helper for map_2tuples 
* @details F is a functor taking the index of tuple element as the 1st argument and
* returns the same type for each element of the input tuples
*
* @tparam F Functor type
* @tparam Tuple1 tuple 1 type
* @tparam Tuple2 tuple 2 type
* @tparam I Index sequence wrapper
* @tparam i... index for the tuples
* @param f Functor
* @param t1 Tuple 1
* @param t2 Tuple 2
* @param I index sequence
*
* @return a std::array the returns of F
*/
template < class F, class Tuple1, class Tuple2, template <size_t...> class I, size_t ... i, class ... Args >
auto map_2tuples_impl(F&& f, Tuple1&& t1, Tuple2&& t2, const I<i...>& Dummy, Args&&... args) ->
typename std::enable_if
< 
  check_if<bare<Tuple1>, can_be_get>::value && 
  check_if<bare<Tuple2>, can_be_get>::value &&
  check_if<bare<F>, is_functor_with_args<size_t, decltype(std::get<0>(std::forward<Tuple1>(t1))), decltype(std::get<0>(std::forward<Tuple2>(t2))), Args...>>::value &&
  are_same<decltype(f(i,std::get<i>(std::forward<Tuple1>(t1)),std::get<i>(std::forward<Tuple2>(t2)),std::forward<Args>(args)...))...>::value,
  std::array<decltype(f(0,std::get<0>(std::forward<Tuple1>(t1)),std::get<0>(std::forward<Tuple2>(t2)),std::forward<Args>(args)...)), sizeof...(i)>
>::type 
{
  using T = typename std::common_type<decltype(f(i,std::get<i>(std::forward<Tuple1>(t1)),std::get<i>(std::forward<Tuple2>(t2)),std::forward<Args>(args)...))...>::type;
  constexpr size_t N = sizeof...(i); 
  return std::array<T,N>{ f(i,std::get<i>(std::forward<Tuple1>(t1)),std::get<i>(std::forward<Tuple2>(t2)),std::forward<Args>(args)...)... };
}

/**
* @brief Helper for map_2tuples 
* @details F is a functor taking just two tuple elements and
* returns the same type for each element of the input tuples
*
* @tparam F Functor type
* @tparam Tuple1 tuple 1 type
* @tparam Tuple2 tuple 2 type
* @tparam I Index sequence wrapper
* @tparam i... index for the tuples
* @param f Functor
* @param t1 Tuple 1
* @param t2 Tuple 2
* @param I index sequence
*
* @return a std::array the returns of F
*/
template < class F, class Tuple1, class Tuple2, template <size_t...> class I, size_t ... i, class ... Args >
auto map_2tuples_impl(F&& f, Tuple1&& t1, Tuple2&& t2, const I<i...>& Dummy, Args&&... args) ->
typename std::enable_if
< 
  check_if<bare<Tuple1>, can_be_get>::value && 
  check_if<bare<Tuple2>, can_be_get>::value &&
  check_if<bare<F>, is_functor_with_args<decltype(std::get<0>(std::forward<Tuple1>(t1))), decltype(std::get<0>(std::forward<Tuple2>(t2))), Args...>>::value &&
  are_same<decltype(f(std::get<i>(std::forward<Tuple1>(t1)),std::get<i>(std::forward<Tuple2>(t2)),std::forward<Args>(args)...))...>::value,
  std::array<decltype(f(std::get<0>(std::forward<Tuple1>(t1)),std::get<0>(std::forward<Tuple2>(t2)),std::forward<Args>(args)...)), sizeof...(i)>
>::type 
{
  using T = typename std::common_type<decltype(f(std::get<i>(std::forward<Tuple1>(t1)),std::get<i>(std::forward<Tuple2>(t2)),std::forward<Args>(args)...))...>::type;
  constexpr size_t N = sizeof...(i); 
  return std::array<T,N>{ f(std::get<i>(std::forward<Tuple1>(t1)),std::get<i>(std::forward<Tuple2>(t2)),std::forward<Args>(args)...)... };
}

/**
* @brief Helper for map_2tuples
* @details F is a functor taking the index of tuple element as the 1st argument and
* return different types depending on which tuple elments are provided as its arguments
* @tparam F Functor type
* @tparam Tuple1 tuple 1 type
* @tparam Tuple2 tuple 2 type
* @tparam I Index sequence wrapper
* @tparam i... index for the tuples
* @param f Functor
* @param t1 Tuple 1
* @param t2 Tuple 2
* @param I index sequence
*
* @return a new tuple with each elment being the return of F
*/
template < class F, class Tuple1, class Tuple2, template <size_t...> class I, size_t ... i, class ... Args >
auto map_2tuples_impl(F&& f, Tuple1&& t1, Tuple2&& t2, const I<i...>& Dummy, Args&&... args) ->
typename std::enable_if
< 
  check_if<bare<Tuple1>, can_be_get>::value && 
  check_if<bare<Tuple2>, can_be_get>::value &&
  check_if<bare<F>, is_functor_with_args<size_t, decltype(std::get<0>(std::forward<Tuple1>(t1))), decltype(std::get<0>(std::forward<Tuple2>(t2))), Args...>>::value &&
  !are_same<decltype(f(i,std::get<i>(std::forward<Tuple1>(t1)),std::get<i>(std::forward<Tuple2>(t2)),std::forward<Args>(args)...))...>::value,
  decltype(std::make_tuple( f(i,std::get<i>(std::forward<Tuple1>(t1)),std::get<i>(std::forward<Tuple2>(t2)),std::forward<Args>(args)...)... ))
>::type 
{
  return std::make_tuple( f(i,std::get<i>(std::forward<Tuple1>(t1)),std::get<i>(std::forward<Tuple2>(t2)),std::forward<Args>(args)...)... );
}

/**
* @brief Helper for map_2tuples
* @details F is a functor taking two tuple elements as the 1st two arguments and
* return different types depending on which tuple elments are provided as its arguments
* @tparam F Functor type
* @tparam Tuple1 tuple 1 type
* @tparam Tuple2 tuple 2 type
* @tparam I Index sequence wrapper
* @tparam i... index for the tuples
* @param f Functor
* @param t1 Tuple 1
* @param t2 Tuple 2
* @param I index sequence
*
* @return a new tuple with each elment being the return of F
*/
template < class F, class Tuple1, class Tuple2, template <size_t...> class I, size_t ... i, class ... Args >
auto map_2tuples_impl(F&& f, Tuple1&& t1, Tuple2&& t2, const I<i...>& Dummy, Args&&... args) ->
typename std::enable_if
< 
  check_if<bare<Tuple1>, can_be_get>::value && 
  check_if<bare<Tuple2>, can_be_get>::value &&
  check_if<bare<F>, is_functor_with_args<decltype(std::get<0>(std::forward<Tuple1>(t1))), decltype(std::get<0>(std::forward<Tuple2>(t2))), Args...>>::value &&
  !are_same<decltype(f(std::get<i>(std::forward<Tuple1>(t1)),std::get<i>(std::forward<Tuple2>(t2)),std::forward<Args>(args)...))...>::value,
  decltype(std::make_tuple( f(std::get<i>(std::forward<Tuple1>(t1)),std::get<i>(std::forward<Tuple2>(t2)),std::forward<Args>(args)...)... ))
>::type 
{
  return std::make_tuple( f(std::get<i>(std::forward<Tuple1>(t1)),std::get<i>(std::forward<Tuple2>(t2)),std::forward<Args>(args)...)... );
}

/**
* @brief Helper for map_2tuples 
* @details input tuples are empty tuples, i.e. std::tuple_size<Tuple1>::value == 0 && std::tuple_size<Tuple2>::value == 0 
*
* @tparam F Functor type
* @tparam Tuple1 tuple 1 type
* @tparam Tuple2 tuple 2 type
* @tparam I Index sequence wrapper
* @tparam i... index for the tuples
* @param f Functor
* @param t1 Tuple 1
* @param t2 Tuple 2
* @param I index sequence
*
* @return an empty tuple 
*/
template < class F, class Tuple1, class Tuple2, template <size_t...> class I, class ... Args >
std::tuple<> map_2tuples_impl(F&& f, Tuple1&& t1, Tuple2&& t2, const I<>& Dummy, Args&&... args)
{
  return std::tuple<>{};
}

/**
* @brief map two tuple of the same sizes to another tuple 
*
* @tparam F Functor type
* @tparam Tuple1 Tuple 1 type
* @tparam Tuple2 Tuple 2 type 
* @param f functor
* @param t1 Tuple 1
* @param t2 Tuple 2
*
* @return a new tuple with each elment being the return of F
*/
template < class F, class Tuple1, class Tuple2, class ... Args >
auto map_2tuples(F&& f, Tuple1&& t1, Tuple2&& t2, Args&&... args) -> 
decltype
(
  map_2tuples_impl
  (
    std::forward<F>(f), 
    std::forward<Tuple1>(t1), 
    std::forward<Tuple2>(t2), 
    container_index<Min(std::tuple_size<bare<Tuple1>>::value, std::tuple_size<bare<Tuple2>>::value)>{}, 
    std::forward<Args>(args)...
  )
)
{
  constexpr size_t N = Min(std::tuple_size<bare<Tuple1>>::value, std::tuple_size<bare<Tuple2>>::value);
  return map_2tuples_impl(std::forward<F>(f), std::forward<Tuple1>(t1), std::forward<Tuple2>(t2), container_index<N>{});
};

/**
* @brief Helper for multiply_2tuples 
* @details the product types are the same for all the 
* elements in the tuples
*
* @tparam Tuple1 tuple 1 type
* @tparam Tuple2 tuple 2 type
* @tparam I Index sequence wrapper
* @tparam i... index for the tuples
* @param f Functor
* @param t1 Tuple 1
* @param t2 Tuple 2
* @param I index sequence
*
* @return a std::array of the elements 
*/
template < class Tuple1, class Tuple2, template <size_t...> class I, size_t ... i >
auto multiply_2tuples_impl(Tuple1&& t1, Tuple2&& t2, const I<i...>&) ->
typename std::enable_if
< 
  check_if<bare<Tuple1>, can_be_get>::value && 
  check_if<bare<Tuple2>, can_be_get>::value &&
  are_same<decltype(std::get<i>(std::forward<Tuple1>(t1))*std::get<i>(std::forward<Tuple2>(t2)))...>::value,
  std::array<decltype(std::get<0>(std::forward<Tuple1>(t1))*std::get<0>(std::forward<Tuple2>(t2))), sizeof...(i)>
>::type 
{
  using T = typename std::common_type<decltype(std::get<i>(std::forward<Tuple1>(t1))*std::get<i>(std::forward<Tuple2>(t2)))...>::type;
  constexpr size_t N = sizeof...(i);
  return std::array<T,N>{ std::get<i>(std::forward<Tuple1>(t1)) * std::get<i>(std::forward<Tuple2>(t2))... };
};

/**
* @brief Helper for multiply_2tuples
* @details the product types are different for all the 
* elements in the tuples
*
* @tparam Tuple1 tuple 1 type
* @tparam Tuple2 tuple 2 type
* @tparam I Index sequence wrapper
* @tparam i... index for the tuples
* @param f Functor
* @param t1 Tuple 1
* @param t2 Tuple 2
* @param I index sequence
*
* @return a new tuple with each elment being the return the corresponding multiplication 
*/
template < class Tuple1, class Tuple2, template <size_t...> class I, size_t ... i >
auto multiply_2tuples_impl(Tuple1&& t1, Tuple2&& t2, const I<i...>&) ->
typename std::enable_if
< 
  check_if<bare<Tuple1>, can_be_get>::value && 
  check_if<bare<Tuple2>, can_be_get>::value &&
  !are_same<decltype(std::get<i>(std::forward<Tuple1>(t1))*std::get<i>(std::forward<Tuple2>(t2)))...>::value,
  decltype(std::make_tuple( std::get<i>(std::forward<Tuple1>(t1))*std::get<i>(std::forward<Tuple2>(t2))...) )
>::type 
{
  return std::make_tuple( std::get<i>(std::forward<Tuple1>(t1)) * std::get<i>(std::forward<Tuple2>(t2))... );
};

/**
* @brief Helper for multiply_2tuples 
* @details input tuples are empty tuples, i.e. std::tuple_size<Tuple1>::value == 0 && std::tuple_size<Tuple2>::value == 0 
*
* @tparam Tuple1 tuple 1 type
* @tparam Tuple2 tuple 2 type
* @tparam I Index sequence wrapper
* @tparam i... index for the tuples
* @param f Functor
* @param t1 Tuple 1
* @param t2 Tuple 2
* @param I index sequence
*
* @return a empty tuple 
*/
template < class Tuple1, class Tuple2, template <size_t...> class I >
std::tuple<> multiply_2tuples_impl(Tuple1&& t1, Tuple2&& t2, const I<>&)
{
  return std::tuple<>{}; 
};

/**
* @brief map 2 tuples to a new one by multiplying the corresponding elements 
*
* @tparam Tuple1 Tuple 1 type
* @tparam Tuple2 Tuple 2 type 
* @param f functor
* @param t1 Tuple 1
* @param t2 Tuple 2
*
* @return a new tuple of the same size as the smaller of two input tuples
* with each elment being the return the corresponding multiplication 
*/
template < class Tuple1, class Tuple2 >
auto multiply_2tuples(Tuple1&& t1, Tuple2&& t2) -> 
decltype(
  multiply_2tuples_impl(
    std::forward<Tuple1>(t1), 
    std::forward<Tuple2>(t2), 
    container_index<Min(std::tuple_size<bare<Tuple1>>::value,std::tuple_size<bare<Tuple2>>::value)>{}
    )
  )
{
  static constexpr size_t N = Min(std::tuple_size<bare<Tuple1>>::value,std::tuple_size<bare<Tuple2>>::value);
  return multiply_2tuples_impl(std::forward<Tuple1>(t1), std::forward<Tuple2>(t2), container_index<N>{});
};

#endif   /* ----- #ifndef Tuple_INC  ----- */
