#ifndef  NDArray_INC
#define  NDArray_INC

/*
 * =====================================================================================
 *
 *       Filename:  NDArray.hpp
 *
 *    Description:  Define compile-time multidimensional array
 *
 *        Version:  1.0
 *        Created:  02/25/2015 04:44:26 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */

#include <type_traits>
#include "Typetraits.hpp"
#include <utility>
#include <functional>
#include "indexseq_generator.hpp"

/**
* @brief Helper class to initialize NDArray (original template)
* @details NDArray instantiates a index sequence for the elements of current dimension
* and this helper class performs pack expansion of the index seqence and initialze the 
* elements. See the specializations for more details
* @tparam arr The requested array type
* @tparam NestedType The NDArray type in next level of recursion
* @tparam T... Specialization specific parameters
*/
template<class arr, template <class...> class NestedType, class ... T> struct NDArray_Helper;

/**
* @brief Specialize for dynamic array (user interface)
*
* @tparam arr The requested array type
* @tparam NestedType The NDArray type in next level of recursion
* @tparam I wrapper for the dimension specifiers
* @tparam T... dimension specifiers
* @tparam Ic wrapper for the indices of current dimension
* @tparam newindices... indices of current dimension
*/
template
<
  class arr, 
  template<class, class> class NestedType,
  template <class...> class I, 
  class ... T, 
  template <size_t...> class Ic, 
  size_t ... newindices
> 
struct NDArray_Helper<arr, NestedType, I<T...>, Ic<newindices...>> {
  static_assert(std::is_literal_type<arr>::value, "array type must be literal");
  static constexpr arr value = { NestedType<T, Ic<newindices>>::value... };
};

/**
* @brief Specialize for dynamic array (user interface)
*
* @tparam arr The requested array type
* @tparam NestedType The NDArray type in next level of recursion
* @tparam I wrapper for the dimension specifiers
* @tparam T... dimension specifiers
* @tparam Ic wrapper for the indices of current dimension
* @tparam newindices... indices of current dimension
* @tparam newindices... indices of upper dimensions
*/
template
<
  class arr,
  template<class, class> class NestedType,
  template <class...> class I,
  class ... T,
  template <size_t...> class Ic,
  size_t ... newindices,
  size_t ... oldindices
> 
struct NDArray_Helper<arr, NestedType, I<T...>, Ic<newindices...>, Ic<oldindices...>> {
  static_assert(std::is_literal_type<arr>::value, "array type must be literal");
  static constexpr arr value = { NestedType<T, Ic<newindices, oldindices...>>::value... };
};

/**
* @brief Helper for static array (for user-interface)
* @details The NDArray class creates a index sequence (newindices) to be expanded here
* @tparam arr The array type
* @tparam NestedType The type for next level of recursion
* @tparam I wrapper of the index sequence
* @tparam newindices The indices for current dimension of the array
*/
template
<
  class arr, 
  template<class> class NestedType,
  template <size_t...> class I, 
  size_t ... newindices
> 
struct NDArray_Helper<arr, NestedType, I<newindices...>> {
  static_assert(std::is_literal_type<arr>::value, "array type must be literal");
  static constexpr arr value = { NestedType<I<newindices>>::value... };
};

/**
* @brief Helper for static array (for propagation)
* @details The NDArray class creates a index sequence (newindices) to be expanded here and accumulated into oldindices
* @tparam arr The array type
* @tparam NestedType The type for next level of recursion
* @tparam I wrapper of the index sequence
* @tparam newindices The indices for current dimension of the array
* @tparam oldindices The indices for upper dimensions
*/
template
<
  class arr, 
  template<class> class NestedType,
  template <size_t...> class I, 
  size_t ... newindices,
  size_t ... oldindices
> 
struct NDArray_Helper<arr, NestedType, I<newindices...>, I<oldindices...>> {
  static_assert(std::is_literal_type<arr>::value, "array type must be literal");
  static constexpr arr value = { NestedType<I<newindices, oldindices...>>::value... };
};

/**
* @brief Helper for runtime static array (user interface)
* @details Basically the same as static array of literal type except that 
* we need an initializer to instantiate the value along the recursion
*
* @tparam arr The array type
* @tparam NestedType The type for next level of recursion
* @tparam I wrapper of the index sequence
* @tparam newindices The indices for current dimension of the array
* @tparam Tinit The initializer type that's passed along the recursion and 
* finally used to instantiate an object of Fn
*
*/
template
<
  class Trefarr, 
  template <class> class NestedType,
  template <size_t...> class I, 
  size_t ... newindices,
  class Tinit
> 
struct NDArray_Helper<Trefarr, NestedType, I<newindices...>, Tinit> {
  Tinit _init;
  NDArray_Helper(const Tinit& init) : _init(init) {}; 
  const Trefarr value = {{ NestedType<I<newindices>>(_init).value... }};
};

/**
* @brief Helper for runtime static array (propagation)
*
* @tparam arr The array type
* @tparam NestedType The type for next level of recursion
* @tparam I wrapper of the index sequence
* @tparam newindices The indices for current dimension of the array
* @tparam oldindices The indices for upper dimensions
* @tparam Tinit The initializer type that's passed along the recursion and 
* finally used to instantiate an object of Fn
*/
template
<
  class Trefarr, 
  template<class> class NestedType,
  template <size_t...> class I, 
  size_t ... newindices,
  size_t ... oldindices,
  class Tinit
> 
struct NDArray_Helper<Trefarr, NestedType, I<newindices...>, I<oldindices...>, Tinit> {
  Tinit _init;
  NDArray_Helper(const Tinit& init) : _init(init) {}; 
  const Trefarr value = {{ NestedType<I<newindices, oldindices...>>(_init).value... }};
};


/**
* @brief Original template for N-dimensional array initializer 
* @details See the partial specialization for more details
* @tparam Arr The template class that holds array
* @tparam Fn The function wrapper class that maps the indices of the elements in the array to a value
* @tparam T A list of types involved
*/
template 
< 
  template <class,size_t...> class Arr,
  template<class,size_t...> class Fn,
  class ... T
> 
struct NDArray;

/**
* @brief Partial specialize for dynamic array (user interface)
* @details Here we allow the sizes of all the subdimensions to vary, e.g., we
* can represent a triangular matrix by {{1,2}, {3}}. The user need to supply
* the sizes for all the subdimensions, e.g., I<In<2>, In<1>> for above example
*
* @tparam Arr The template class that holds array
* @tparam Fn The function wrapper class that maps the indices of the elements
* in the array to a value
* @tparam Elem The underlying element type of the array
* @tparam I Wrapper class that wraps a list of types
* @tparam Ii Wrapper class or index sequence specifying the sizes of the first
* subdimesnion
* @tparam Ij... Wrapper classes or index sequence specifying the sizes of the
* first subdimesnions
*/
template
< 
  template <class> class Arr,
  template<class,size_t...> class Fn,
  class Elem,
  template <class...> class I,
  class Ii,
  class ... Ij 
>
struct NDArray<Arr, Fn, Elem, I<Ii, Ij...>> {
  static_assert(std::is_literal_type<Elem>::value && std::is_literal_type<Arr<Elem>>::value, "Element and array types must both be literal");
  static constexpr size_t N = sizeof...(Ij) + 1;
  using IndicesOfThisDimension = container_index<N>;
  template<class Ipack, class Indices> using NestedType = NDArray<Arr, Fn, Elem, Ipack, Indices>; 
  using NestedArrType = typename NestedType<Ii, container_index<1>>::type; //TODO: This might be avoided by using decltype()
  using type = Arr<NestedArrType>;
  static constexpr type value = NDArray_Helper<type, NestedType, I<Ii, Ij...>, IndicesOfThisDimension>::value; 
};


/**
* @brief Partial specialize for dynamic array (propagation)
*
* @tparam Arr The template class that holds array
* @tparam Fn The function wrapper class that maps the indices of the elements
* in the array to a value
* @tparam Elem The underlying element type of the array
* @tparam I Wrapper class that wraps a list of types
* @tparam Ii Wrapper class or index sequence specifying the sizes of the first
* subdimesnion
* @tparam Ij... Wrapper classes or index sequence specifying the sizes of the
* first subdimesnions
* @tparam Ic Wrapper class for indidces
* @tparam i... The indices for the upper dimensions
*/
template 
< 
  template <class> class Arr,
  template<class,size_t...> class Fn,
  class Elem,
  template <class...> class I,
  class Ii,
  class ... Ij,
  template <size_t...> class Ic,
  size_t ... i 
>
struct NDArray<Arr, Fn, Elem, I<Ii, Ij...>, Ic<i...>> {
  static_assert(std::is_literal_type<Elem>::value && std::is_literal_type<Arr<Elem>>::value, "Element and array types must both be literal");
  static constexpr size_t N = sizeof...(Ij) + 1;
  using IndicesOfThisDimension = container_index<N>;
  template<class Ipack, class Indices> using NestedType = NDArray<Arr, Fn, Elem, Ipack, Indices>;
  using NestedArrType = typename NestedType<Ii, Ic<0, i...>>::type;
  using type = Arr<NestedArrType>;
  static constexpr type value = NDArray_Helper<type, NestedType, I<Ii, Ij...>, IndicesOfThisDimension, Ic<i...>>::value;
};

/**
* @brief Partial specialize for multidimensional dynamic array (pre-termination)
*
* @tparam Arr The template class that holds array  
* @tparam Fn The function wrapper class that maps the indices of the elements
* in the array to a value
* @tparam Elem The underlying element type of the array
* @tparam I A wrapper of the sizes for each dimension of the array. The number
* from left to right denotes the sizes for the outter to inner dimension
* @tparam Ni The size of the current most dimension
* @tparam i... The indices for the upper dimensions
*/
template
< 
  template <class> class Arr,
  template<class,size_t...> class Fn,
  class Elem,
  template <size_t...> class I,
  size_t Ni,
  size_t ... i
>
struct NDArray<Arr, Fn, Elem, I<Ni>, I<i...>> {
  static_assert(std::is_literal_type<Elem>::value && std::is_literal_type<Arr<Elem>>::value, "Element and array types must both be literal");
  using IndicesOfThisDimension = container_index<Ni>;
  template<class Indices> using NestedType = NDArray<Arr, Fn, Elem, I<>, Indices>; 
  using NestedArrType = typename NestedType<I<0, i...>>::type;
  using type = Arr<NestedArrType>;
  static constexpr type value = NDArray_Helper<type, NestedType, IndicesOfThisDimension, I<i...>>::value; 
};

/**
* @brief Partial specialize for multidimensional static array (user interface)
* @details The user need to supply the array type, the function that maps the
* indices to the value of each element and the sizes of each dimension
*
* @tparam Arr The template class that holds array (which takes 2 template
* parameters itself: a type of the elements and the number of elements 
* @tparam Fn The function wrapper class that maps the indices of the elements
* in the array to a value
* @tparam Elem The underlying element type of the array
* @tparam I A wrapper of the sizes for each dimension of the array. The number
* from left to right denotes the sizes for the outter to inner dimension
* @tparam Ni The size of the outter most dimension
* @tparam Nj... The sizes of the rest of the dimensions
*/
template
< 
  template <class,size_t> class Arr,
  template<class,size_t...> class Fn,
  class Elem,
  template <size_t...> class I,
  size_t Ni,
  size_t ... Nj
>
struct NDArray<Arr, Fn, Elem, I<Ni,Nj...>> {
  static_assert(std::is_literal_type<Elem>::value && std::is_literal_type<Arr<Elem,1>>::value, "Element and array types must both be literal");
  using IndicesOfThisDimension = container_index<Ni>;
  template<class Indices> using NestedType = NDArray<Arr, Fn, Elem, I<Nj...>, Indices>; 
  using NestedArrType = typename NestedType<I<0>>::type;
  using type = Arr<NestedArrType, Ni>;
  static constexpr type value = NDArray_Helper<type, NestedType, IndicesOfThisDimension>::value; 
};

/**
* @brief Partial specialize for multidimensional static array (propagation)
*
* @tparam Arr The template class that holds array (which takes 2 template
* parameters itself: a type of the elements and the number of elements 
* @tparam Fn The function wrapper class that maps the indices of the elements
* in the array to a value
* @tparam Elem The underlying element type of the array
* @tparam I A wrapper of the sizes for each dimension of the array. The number
* from left to right denotes the sizes for the outter to inner dimension
* @tparam Ni The size of the outter most dimension
* @tparam Nj... The sizes of the rest of the dimensions
* @tparam i... The indices for the upper dimensions
*/
template
< 
  template <class,size_t> class Arr,
  template<class,size_t...> class Fn,
  class Elem,
  template <size_t...> class I,
  size_t Ni,
  size_t ... Nj,
  size_t ... i
>
struct NDArray<Arr, Fn, Elem, I<Ni,Nj...>, I<i...>> {
  static_assert(std::is_literal_type<Elem>::value && std::is_literal_type<Arr<Elem,1>>::value, "Element and array types must both be literal");
  using IndicesOfThisDimension = container_index<Ni>;
  template<class Indices> using NestedType = NDArray<Arr, Fn, Elem, I<Nj...>, Indices>; 
  using NestedArrType = typename NestedType<I<0, i...>>::type;
  using type = Arr<NestedArrType, Ni>;
  static constexpr type value = NDArray_Helper<type, NestedType, IndicesOfThisDimension, I<i...>>::value; 
};

/**
* @brief Partial specialize for multidimensional static and dynamic array (termination)
* @details At this point, the index sequence for the sizes has been exhausted
* and we end up with a sequnce of indices for the current element of the array
* and we map the indices to the value of the element using Fn
*
* @tparam Arr The template class that holds array (which takes 2 template
* parameters itself: a type of the elements and the number of elements 
* @tparam Fn The function wrapper class that maps the indices of the elements
* in the array to a value
* @tparam Elem The underlying element type of the array
* @tparam I A wrapper of the sizes for each dimension of the array. The number
* from left to right denotes the sizes for the outter to inner dimension
* @tparam i... The indices for the current element (going from left to right 
* means going from deeper to upper dimensions)
*/
template 
< 
  template <class,size_t...> class Arr,
  template<class,size_t...> class Fn,
  class Elem,
  template <size_t...> class I,
  size_t ... i 
>
struct NDArray<Arr, Fn, Elem, I<>, I<i...>> {
  static_assert(std::is_literal_type<Elem>::value && std::is_literal_type<Arr<Elem,1>>::value, "Element and array types must both be literal");
  using type = Elem;
  static constexpr type value = Fn<Elem, i...>::value;  
};

/**
* @brief Partial specialization for runtime static array (user interface)
* @details The idea is basically the same as in the case of literal-type specialization. 
* The only difference is that the element type is initialized at runtime, which could be non-literal, 
* the user also has to supply an initializer type (Tinit) and one of its instance (init) 
* in order to initialize the array
*
* @tparam Arr The template class that holds array (which takes 2 template
* parameters itself: a type of the elements and the number of elements 
* @tparam Fn The function wrapper class that maps the indices of the elements
* in the array to a value
* @tparam Elem The underlying element type of the array
* @tparam I A wrapper of the sizes for each dimension of the array. The number
* from left to right denotes the sizes for the outter to inner dimension
* @tparam Ni The size of the outter most dimension
* @tparam Nj... The sizes of the rest of the dimensions
* @tparam Tinit The initializer type that's passed along the recursion and 
* finally used to instantiate an object of Fn
*
*/
template
< 
  template <class,size_t> class Arr,
  template<class,size_t...> class Fn,
  class Elem,
  template <size_t...> class I,
  size_t Ni,
  size_t ... Nj,
  class Tinit
>
struct NDArray<Arr, Fn, Elem, I<Ni,Nj...>, Tinit> {
  using IndicesOfThisDimension = container_index<Ni>;
  template<class Indices> using NestedType = NDArray<Arr, Fn, Elem, I<Nj...>, Indices, Tinit>; 
  using NestedArrType = typename NestedType<I<0>>::type;
  using type = Arr<NestedArrType, Ni>;
  
  Tinit _init; 
  /**
  * @brief Constructor
  *
  * @param init instance of the initializer
  */
  NDArray(const Tinit& init) : _init(init) {};
  using Thelper = NDArray_Helper<type, NestedType, IndicesOfThisDimension, Tinit>;
  const type value = (Thelper(this->_init)).value; 
};

/**
* @brief Partial specialization for runtime static array (propagation)
*
* @tparam Arr The template class that holds array (which takes 2 template
* parameters itself: a type of the elements and the number of elements 
* @tparam Fn The function wrapper class that maps the indices of the elements
* in the array to a value
* @tparam Elem The underlying element type of the array
* @tparam I A wrapper of the sizes for each dimension of the array. The number
* from left to right denotes the sizes for the outter to inner dimension
* @tparam Ni The size of the outter most dimension
* @tparam Nj... The sizes of the rest of the dimensions
* @tparam i... Indices of the outter most dimension
* @tparam Tinit The initializer type that's passed along the recursion and 
* finally used to instantiate an object of Fn
*/
template
< 
  template <class,size_t> class Arr,
  template<class,size_t...> class Fn,
  class Elem,
  template <size_t...> class I,
  size_t Ni,
  size_t ... Nj,
  size_t ... i,
  class Tinit
>
struct NDArray<Arr, Fn, Elem, I<Ni,Nj...>, I<i...>, Tinit> {
  using IndicesOfThisDimension = container_index<Ni>;
  template<class Indices> using NestedType = NDArray<Arr, Fn, Elem, I<Nj...>, Indices, Tinit>; 
  using NestedArrType = typename NestedType<I<0, i...>>::type;
  using type = Arr<NestedArrType, Ni>;

  Tinit _init; 
  /**
  * @brief Constructor
  *
  * @param init instance of the initializer
  */
  NDArray(const Tinit& init) : _init(init) {};
  using Thelper = NDArray_Helper<type, NestedType, IndicesOfThisDimension, I<i...>, Tinit>;
  const type value = Thelper(_init).value; 
};

/**
* @brief Partial specialization for runtime static array (termination)
*
* @tparam Arr The template class that holds array (which takes 2 template
* parameters itself: a type of the elements and the number of elements 
* @tparam Fn The function wrapper class that maps the indices of the elements
* in the array to a value
* @tparam Elem The underlying element type of the array
* @tparam I A wrapper of the sizes for each dimension of the array. The number
* from left to right denotes the sizes for the outter to inner dimension
* @tparam i... Indices of the element 
* @tparam Tinit The initializer type that's passed along the recursion and 
* finally used to instantiate an object of Fn
*/
template 
< 
  template <class,size_t...> class Arr,
  template<class,size_t...> class Fn,
  class Elem,
  template <size_t...> class I,
  size_t ... i,
  class Tinit
>
struct NDArray<Arr, Fn, Elem, I<>, I<i...>, Tinit> {
  using type = Elem;
    
  Tinit _init; 
  /**
  * @brief Constructor
  *
  * @param init instance of the initializer
  */
  NDArray(const Tinit& init) : _init(init) {};
  using Tfn = Fn<Elem, i...>;
  const type value = Tfn{_init}.value;  
};


/**
* @brief initialize all elements of an NDArray to zero if used as NDArray parameter
*
* @tparam T element type
* @tparam i dummy parameters
*/
template <class T, size_t ... i> struct init_to_zeros {
  static constexpr T value = T{0};
};

#endif   /* ----- #ifndef NDArray_INC  ----- */
