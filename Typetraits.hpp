#ifndef  Typetraits_INC
#define  Typetraits_INC
/*
 * =====================================================================================
 *
 *       Filename:  typetraits.hpp
 *
 *    Description:  Anything related to type traits
 *
 *        Version:  1.0
 *        Created:  02/16/2015 05:06:30 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */

#include <type_traits>
#include <utility>
#include <functional>

/*
 * =====================================================================================
 *        Class:  check_if
 *  Description:  check if a class T satisify a condition 
 *                The idea is to make the compiler decide which member function template
 *                'f' to overload based on whether the parameter 'Condition' has a 
 *                valid definition of nested template class 'type' -- if yes, then 
 *                std::true_type f(int*) is seen because interpreting '0' in f<T>(0) 
 *                as int* is more specific than the vararg list '...'.
 *                The implementation of class Condition for checking if the class T has 
 *                a member function 'Mfn' could look like this:
 *                template <class Ret, class ... Arg>
 *                struct Condition {
 *                  template < class T, Ret (T::*)(Arg&&...) = &T::Mfn > struct type {};
 *                };
 * =====================================================================================
 */
template < class T, class Condition >
struct check_if
{
  template < class C > 
  static auto f(int*)
    ->
    decltype(
      //NOTE: the first expression (left operand in comma operator) is discarded
      std::declval<typename Condition::template type<C>>(), 
      //NOTE: std::true_type is used if the left operand is valid
      std::true_type() 
    ); 

  template < class C > static auto f(...) -> std::false_type;

  constexpr static bool value = decltype(f<T>(0))::value;
}; /* ----------  end of template class check_if  ---------- */

//Here are some examples of the Condition class
//! This class check if a class T has member function 
// T::size_type size() const 
struct has_memfn_size {
  template < class T, typename T::size_type (T::*)() const = &T::size >
  struct type {};
};

//! This class check if a class T has member function 
// void emplace_back(T::value_type&) 
struct has_memfn_emplace_back {
  template < class T, void (T::*)(typename T::value_type&) = &T::emplace_back >
  struct type {};
};

//! This class check if a class T has member function 
// T::type get() const 
struct has_memfn_get {
  template < class T, typename T::type& (T::*)() const = &T::get >
  struct type {};
};

//!This class check if a class T can be iterated using nested const_iterator 
struct can_be_const_iterate {
  template < class T, 
  class C = typename T::const_iterator,
  C (T::*)() const = &T::begin,
  C (T::*)() const = &T::end
           >
  struct type {};
};

//! This class check if a class T has unerlying elements
//that can be retrieved by std::get<0>
//NOTE: This depends on the STL header <functional> to work although 
//I'm not very sure why. My guess is somehow std::function help wrap
//std::get<0> in the default value of the second template parameter 
//of can_be_get::template <T, func> type
struct can_be_get {
  template < class T, decltype(std::get<0>(std::declval<T>()))& (*)(T&) = &std::get<0> >
  struct type {};
};

template <class S>
struct can_be_streamed {
  template < class T >
  using type = std::remove_reference<decltype(std::declval<S&>() << std::declval<T>())>;
};

template <class ... Args>
struct is_functor_with_args {
  template < class T, decltype(std::declval<T>().operator()(std::declval<Args>()...)) (T::*)(Args...) = &T::operator() >
  struct type{};
};

/*
 * =====================================================================================
 *        Class:  Type trait of std::string 
 *  Description:  
 * =====================================================================================
 */
#include <string>

template < class T > struct is_string : std::false_type {};
template <typename charT, typename traits, typename Alloc>
struct is_string<std::basic_string<charT, traits, Alloc> > : std::true_type {};

/*
 * =====================================================================================
 *        Class:  Type trait of the default floating point number, i.e., float, double or long double 
 *  Description:  
 * =====================================================================================
 */
template < class T > struct is_stdfloat : std::false_type {};
template <>
struct is_stdfloat<float> : std::true_type {};
template <>
struct is_stdfloat<double> : std::true_type {};
template <>
struct is_stdfloat<long double> : std::true_type {};

/*
 * =====================================================================================
 *        Class:  Type trait of shared_ptr 
 *  Description:  
 * =====================================================================================
 */
#include <memory>

namespace std {
template < class T >
struct is_pointer<std::shared_ptr<T>> : std::true_type {};
}

/*
 * =====================================================================================
 *        Class:  Are_same
 *  Description:  Check if all the types are the same
 * =====================================================================================
 */
template < class ... T > struct are_same;
template < class T > struct are_same<T> { static constexpr bool value = true; };
template < class Ti, class Tj, class ... Tk > 
struct are_same<Ti, Tj, Tk...> {
  static constexpr bool value = std::is_same<Ti, Tj>::value && are_same<Tj, Tk...>::value;
};

template < class Ti, class Tj >
struct are_same<Ti, Tj> {
  static constexpr bool value = std::is_same<Ti,Tj>::value;
};

/*
 * =====================================================================================
 *        Class:  bare
 *  Description:  remove cv-qualifiers and references from a class (similar to std::decay
 *  except that it doesn't handle function/pointer here)
 * =====================================================================================
 */
template < class T > using bare = typename std::remove_cv<typename std::remove_reference<T>::type>::type;

/*
 * =====================================================================================
 *        Class:  same_underlying_type
 *  Description:  same if the bare types (with cv-qualifiers and references removed) are the same  
 * =====================================================================================
 */
template < class ... T > using same_underlying_type = are_same<bare<T>...>;

/*
 * =====================================================================================
 *        Class:  common 
 *  Description:  short hand for std::common_type<T...>::type  
 * =====================================================================================
 */
template < class ... T > using common = typename std::common_type<T...>::type;

#endif   /* ----- #ifndef Typetraits_INC  ----- */


