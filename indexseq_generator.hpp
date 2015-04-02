#ifndef  indexseq_generator_INC
#define  indexseq_generator_INC
/*
 * =====================================================================================
 *
 *       Filename:  seq_generator.hpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/16/2015 05:38:05 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */

/*
 * =====================================================================================
 *        Class:  Indices generator for constexpr integer 
 *  Description:  These helper classes define a compile-time sequence 
 *                Indices<N1, N2, ..., > for a given N1, a propagator template P<Indices<N1,...> > 
 *                where Ni+1 == P<Indices<N1,N2,...,Ni> >::value and a termination 
 *                condition T<Indices<N1...> > that determine the end of the sequence.
 *                The idea is to recursively define a nested type in the helper template
 *                class make_index_seq where a new Ni is pushed back in Indices<N1..., Ni-1> 
 *                at each recursion step until the termination condistion is satisfied
 * =====================================================================================
 */
//! These are the propagators
//! Propagate the arithmetic sequence with a commone difference of dN and initial term Ni
template < class T, class Arg > struct ArithmeticSeq{};
template < size_t ... indices, template <size_t...> class I, size_t Ni, size_t dN>
struct ArithmeticSeq<I<indices...>, I<Ni, dN> > {
  typedef I<indices..., Ni+dN> type;
  typedef I<Ni+dN, dN> seed;
};
//! Specialize for index sequence
template < template <size_t...> class I>
struct ArithmeticSeq<I<>, I<> > {
  typedef I<0> type;
  typedef I<0, 1> seed;
};

//! Propagate the Fibonacci sequence
template < class T, class Arg > struct Fibonacci{};
template < size_t ... indices, template <size_t...> class I, size_t Ni, size_t Nj > 
struct Fibonacci<I<indices...>, I<Ni,Nj> > {
  typedef I<indices...> Indices;
  typedef I<indices..., Ni+Nj> type;
  typedef I<Nj, Ni+Nj> seed; 
};

template <template <size_t...> class I >
struct Fibonacci<I<>, I<> > {
  typedef typename Fibonacci<I<1, 1>, I<1,1> >::type type;
  typedef I<1, 2> seed;
};

template <template <size_t...> class I >
struct Fibonacci<I<1>, I<> > {
  typedef typename Fibonacci<I<1, 1>, I<1,1> >::type type;
  typedef I<1, 2> seed;
};

//! These are the termination condition
//! Terminate if a maximum number of elements is reached in the sequence
template < class T, class Arg > struct StopAtMaxN {};
template < size_t ... indices, template <size_t...> class I, size_t Max> 
struct StopAtMaxN<I<indices...>, I<Max> > {
  constexpr static bool decision = sizeof...(indices) >= Max ? true : false;
};

//! This is the target type to be generated
template < size_t ... N > struct IndexSeq {};

//! Primary template
template < bool Terminate, class Indices, class Propagator, class Terminator > 
struct make_index_seq{};
//! Propagate if Terminate == false
template < class Indices, 
	   template <class, class> class Propagator,
	   class PropagatorArg, 
	   template <class, class> class Terminator, 
	   class TerminatorArg 
	 >
struct make_index_seq< false, Indices, Propagator<Indices, PropagatorArg>, Terminator<Indices, TerminatorArg> > {
  typedef typename Propagator<Indices, PropagatorArg>::type newIndices;
  typedef typename Propagator<Indices, PropagatorArg>::seed newPropagatorArg;
  typedef Propagator<newIndices,newPropagatorArg> newPropagator;
  typedef Terminator<newIndices,TerminatorArg> newTerminator;
  constexpr static bool newdecision = newTerminator::decision;
  typedef typename make_index_seq<newdecision, newIndices, newPropagator, newTerminator>::type type; 
};

//! Stop if Terminate == true 
template < class Indices, class Propagator, class Terminator > 
struct make_index_seq<true, Indices, Propagator, Terminator> {
  typedef Indices type;
};

template < size_t N >
using make_container_index = 
  make_index_seq< StopAtMaxN<IndexSeq<>, IndexSeq<N>>::decision,
                  IndexSeq<>,
		  ArithmeticSeq<IndexSeq<>, IndexSeq<>>,
		  StopAtMaxN<IndexSeq<>, IndexSeq<N>> >;

template < size_t N >
using container_index = typename make_container_index<N>::type;

#endif   /* ----- #ifndef indexseq_generator_INC  ----- */
