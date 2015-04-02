#include "SymmetricCartesianTensor.hpp"
#include "minmax.hpp"
#include "Typetraits.hpp"
#include <iostream>
#include <array>

using namespace std;

using real = double;

//! The gradients (up to rank N) of
//f(r^2) = (r^2)^(-1/2) as a function of r^2 = rx^2 + ry^2 + rz^2
//as expected for emp direct-space interaction  
template < size_t N >
array<real, N+1> Frecp (const real& r) {
  array<real, N+1> ans;
  const real r2 = r*r;
  ans[0] = 1/r;
  for(size_t i = 1; i <= N; ++i) {
    ans[i] = -((2*i-1)/(2*r2)) * ans[i-1];
  }
  return ans;
}

template < class EMP1, class EMP2 >
real EMPenergy (EMP1&& m, EMP2&& n, const array<real,3>& r) {
  using Tb = SymmetricCartesianTensorBilinearForm<bare<EMP1>,bare<EMP2>>;
  auto blf = Tb{};
  static constexpr size_t R1 = bare<EMP1>::_rank;
  static constexpr size_t R2 = bare<EMP2>::_rank;
  static constexpr bool D1 = bare<EMP1>::_detraced;
  static constexpr bool D2 = bare<EMP2>::_detraced;
  static constexpr real sign = R2 % 2 ? -1 : 1;
  static constexpr real factor1 = D1 ? 1/real(multifac(2*R1-1,2)) : 1;
  static constexpr real factor2 = D2 ? 1/real(multifac(2*R2-1,2)) : 1;
  auto e = sign * factor1 * factor2 * 
           blf.template BF0<Frecp<R1+R2>>(std::forward<EMP1>(m), std::forward<EMP2>(n),r);
  return e;
}

template < class EMP1, class EMP2 >
array<real,3> EMPforce (EMP1&& m, EMP2&& n, const array<real,3>& r) {
  using Tb = SymmetricCartesianTensorBilinearForm<bare<EMP1>,bare<EMP2>>;
  auto blf = Tb{};
  static constexpr size_t R1 = bare<EMP1>::_rank;
  static constexpr size_t R2 = bare<EMP2>::_rank;
  static constexpr bool D1 = bare<EMP1>::_detraced;
  static constexpr bool D2 = bare<EMP2>::_detraced;
  static constexpr real sign = (R2+1) % 2 ? -1 : 1;
  static constexpr real factor1 = D1 ? 1/real(multifac(2*R1-1,2)) : 1;
  static constexpr real factor2 = D2 ? 1/real(multifac(2*R2-1,2)) : 1;
  auto F = (sign * factor1 * factor2) *
           blf.template BF1<Frecp<R1+R2+1>>(std::forward<EMP1>(m), std::forward<EMP2>(n),r);
  return F;
}

template < class EMP1, class EMP2 >
array<real,3> EMPtorqm (EMP1&& m, EMP2&& n, const array<real,3>& r) {
  using Tb = SymmetricCartesianTensorBilinearForm<bare<EMP1>,bare<EMP2>>;
  auto blf = Tb{};
  static constexpr size_t R1 = bare<EMP1>::_rank;
  static constexpr size_t R2 = bare<EMP2>::_rank;
  static constexpr bool D1 = bare<EMP1>::_detraced;
  static constexpr bool D2 = bare<EMP2>::_detraced;
  const auto bft = blf.template BFT<Frecp<R1+R2>>(std::forward<EMP1>(m), std::forward<EMP2>(n),r);
  static constexpr real factor1 = D1 ? 1/real(multifac(2*R1-1,2)) : 1;
  static constexpr real factor2 = D2 ? 1/real(multifac(2*R2-1,2)) : 1;
  static constexpr real sign = R2 % 2 ? -1 : 1;
  return real(R1) * sign * factor1 * factor2 * bft;
}

template < class EMP1, class EMP2 >
array<real,3> EMPnumforce (EMP1&& m, EMP2&& n, const array<real,3>& r, const real& eps) {
  using Tb = SymmetricCartesianTensorBilinearForm<bare<EMP1>,bare<EMP2>>;
  auto blf = Tb{};
  static constexpr size_t R1 = bare<EMP1>::_rank;
  static constexpr size_t R2 = bare<EMP2>::_rank;
  static constexpr bool D1 = bare<EMP1>::_detraced;
  static constexpr bool D2 = bare<EMP2>::_detraced;
  const array<real, 3> dx = {eps/2, 0, 0};
  const array<real, 3> dy = {0, eps/2, 0};
  const array<real, 3> dz = {0, 0, eps/2};

  const auto rpdx = r + dx;
  const auto rmdx = r - dx;
  const auto Epdx =  EMPenergy(m,n,rpdx);
  const auto Emdx =  EMPenergy(m,n,rmdx);
  const auto Fx = (Emdx - Epdx)/eps;

  const auto rpdy = r + dy;
  const auto rmdy = r - dy;
  const auto Epdy =  EMPenergy(m,n,rpdy);
  const auto Emdy =  EMPenergy(m,n,rmdy);
  const auto Fy = (Emdy - Epdy)/eps;

  const auto rpdz = r + dz;
  const auto rmdz = r - dz;
  const auto Epdz =  EMPenergy(m,n,rpdz);
  const auto Emdz =  EMPenergy(m,n,rmdz);
  const auto Fz = (Emdz - Epdz)/eps;

  return array<real,3>{Fx, Fy, Fz};
}

template < class EMP1, class EMP2 >
array<real,3> EMPnumtorqm (EMP1&& m, EMP2&& n, const array<real,3>& r, const array<real,4>& quatm, const real& eps) {
  const array<real,4> epsw = {eps/2,0,0,0};
  const array<real,4> epsx = {0,eps/2,0,0};
  const array<real,4> epsy = {0,0,eps/2,0};
  const array<real,4> epsz = {0,0,0,eps/2};
  
  const auto qwp = quatm + epsw;
  const auto qwm = quatm - epsw;
  auto mxqwp = m.xform(qwp);
  auto mxqwm = m.xform(qwm);
  const auto ewp = EMPenergy( mxqwp, n, r );
  const auto ewm = EMPenergy( mxqwm, n, r );
  const real Tw = (ewm - ewp)/eps;

  const auto qxp = quatm + epsx;
  const auto qxm = quatm - epsx;
  auto mxqxp = m.xform(qxp);
  auto mxqxm = m.xform(qxm);
  const auto exp = EMPenergy( mxqxp, n, r );
  const auto exm = EMPenergy( mxqxm, n, r );
  const real Tx = (exm - exp)/eps;

  const auto qyp = quatm + epsy;
  const auto qym = quatm - epsy;
  auto mxqyp = m.xform(qyp);
  auto mxqym = m.xform(qym);
  const auto eyp = EMPenergy( mxqyp, n, r );
  const auto eym = EMPenergy( mxqym, n, r );
  const real Ty = (eym - eyp)/eps;

  const auto qzp = quatm + epsz;
  const auto qzm = quatm - epsz;
  auto mxqzp = m.xform(qzp);
  auto mxqzm = m.xform(qzm);
  const auto ezp = EMPenergy( mxqzp, n, r );
  const auto ezm = EMPenergy( mxqzm, n, r );
  const real Tz = (ezm - ezp)/eps;

  const array<real,4> Tg = { Tw, Tx, Ty, Tz }; //torque in global quaternion frame

  //The following transformation is found in
  //T. F. Miller, J. Chem. Phys. 116, 8649 (2002)
  //The jacobian to transform torque from space-fixed quaternion frame
  //to space-fixed Cartesian frame is:
  // 0.5*rlg*S'
  // where S' is the 4*4 matrix representation of quaternion conjugate q* = conjugate(q)
  // where q describes the orientation of the local frame wrsp. to global frame
  // rlg is the 4*4 matrix with 1st row and 1st column are zeroes and  the lower right 
  // 3*3 block is the rotation matrix that transforms vectors from local to 
  // global Cartesian space
  
  const array<real,4> quatm_conjugate = { quatm[0], -quatm[1], -quatm[2], -quatm[3] };

  const auto Tl = quatquat(quatm_conjugate, Tg)*0.5; //torque in local quaternion frame
  const array<real, 3> Tlxyz = {Tl[1], Tl[2], Tl[3]};
  const auto MLG = quat_to_mat(quatm);
  return MLG * Tlxyz;
}

int main(int argc, char* argv[]) {
  const real eps = atof(argv[1]);
  
  const array<real,3> r = {1,2,3};
  const array<real,4> quatm = {1,0,0,0};
  SymmetricCartesianTensor<array, real, 4, false> m = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
  SymmetricCartesianTensor<array, real, 4, false> n = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
  cout << "Energy is " << EMPenergy(m, n, r) << endl;
  const auto F = EMPforce(m, n, r);
  const auto Fnum = EMPnumforce(m, n, r, eps);
  decltype(F) Fperdiff = {(F[0]-Fnum[0])/Fnum[0], (F[1]-Fnum[1])/Fnum[1], (F[2]-Fnum[2])/Fnum[2]};
  const auto Tm = EMPtorqm(m,n,r);
  const auto Tnum = EMPnumtorqm(m, n, r, quatm, eps);
  decltype(Tm) Tperdiff = {(Tm[0]-Tnum[0])/Tnum[0], (Tm[1]-Tnum[1])/Tnum[1], (Tm[2]-Tnum[2])/Tnum[2]};
  cout << "Force on m is " << F[0] << " " << F[1] << " " << F[2] << endl;
  cout << "Difference between Force and Num. Force is " << Fperdiff[0] << " " << Fperdiff[1] << " " << Fperdiff[2] << endl; 
  cout << "Torque on m is " << Tm[0] << " " << Tm[1] << " " << Tm[2] << endl;
  cout << "Difference between Torque and Num. Torque is " << Tperdiff[0] << " " << Tperdiff[1] << " " << Tperdiff[2] << endl; 

}
