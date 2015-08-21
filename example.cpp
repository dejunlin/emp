/*
 * =====================================================================================
 *
 *       Filename:  example.cpp
 *
 *    Description:  generate 100 random monopole, dipole, quadrupole, octupol
 *    and hexadecapole at random positions and compute the numerical forces and
 *    torques and compare them to the analytic ones for each pair of multipoles
 *    generated
 *
 *        Version:  1.0
 *        Created:  08/20/2015 17:06:51 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dejun Lin (DL), dejun.lin@gmail.com
 *   Organization:  Department of Biochemistry and Biophysics, Medical Center, University of Rochester
 *
 * =====================================================================================
 */

#include "SymmetricCartesianTensor.hpp"
#include "minmax.hpp"
#include "Typetraits.hpp"
#include <iostream>
#include <array>
#include <cmath>
#include "NDArray.hpp"
#include <chrono>
#include <random>
#include <iomanip>

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
  static constexpr real factor1 = D1 ? 1/real(multifac(2u*R1-1u,2u)) : 1;
  static constexpr real factor2 = D2 ? 1/real(multifac(2u*R2-1u,2u)) : 1;
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
  static constexpr real factor1 = D1 ? 1/real(multifac(2u*R1-1u,2u)) : 1;
  static constexpr real factor2 = D2 ? 1/real(multifac(2u*R2-1u,2u)) : 1;
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
  static constexpr real factor1 = D1 ? 1/real(multifac(2u*R1-1u,2u)) : 1;
  static constexpr real factor2 = D2 ? 1/real(multifac(2u*R2-1u,2u)) : 1;
  static constexpr real sign = R2 % 2 ? -1 : 1;
  return real(R1) * sign * factor1 * factor2 * bft;
}

template < class EMP1, class EMP2 >
array<real,3> EMPnumforce (EMP1&& m, EMP2&& n, const array<real,3>& r, const real& eps) {
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

//! check the consistency between the analytic and numerical force between a pair of EMP
template < class EMP1, class EMP2 > 
array<real,3> chkEMPforce(EMP1&& m, EMP2&& n, const array<real,3>& r, const real& eps) {
      const auto F = EMPforce(forward<EMP1>(m), forward<EMP2>(n), r);
      const auto Fnum = EMPnumforce(forward<EMP1>(m), forward<EMP2>(n), r, eps);
      return {(F[0]-Fnum[0])/Fnum[0], (F[1]-Fnum[1])/Fnum[1], (F[2]-Fnum[2])/Fnum[2]};
}

//! check the consistency between the analytic and numerical torque (on m) between a pair of EMP
template < class EMP1, class EMP2 > 
array<real,3> chkEMPtorqm(EMP1&& m, EMP2&& n, const array<real,3>& r, const array<real,4>& quat, const real& eps) {
      const auto Tm = EMPtorqm(forward<EMP1>(m), forward<EMP2>(n), r);
      const auto Tmnum = EMPnumtorqm(forward<EMP1>(m), forward<EMP2>(n), r, quat, eps);
      return {(Tm[0]-Tmnum[0])/Tmnum[0], (Tm[1]-Tmnum[1])/Tmnum[1], (Tm[2]-Tmnum[2])/Tmnum[2]};
}


int main(int argc, char* argv[]) {
  const real epsf = atof(argv[1]);
  const real epst = atof(argv[2]);
  
  using hexdecp = SymmetricCartesianTensor<array, real, 4, false>;
  using octp = SymmetricCartesianTensor<array, real, 3, false>;
  using quadp = SymmetricCartesianTensor<array, real, 2, false>;
  using dp = SymmetricCartesianTensor<array, real, 1, false>;
  using mp = SymmetricCartesianTensor<array, real, 0, false>;

  constexpr size_t N = 1000;
  array<hexdecp, N> hps;
  array<octp, N> ops;
  array<quadp, N> qps;
  array<dp, N> dps;
  array<mp, N> mps;
  array<array<real,3>, N*5> pos;
  const array<real,4> quatm = {1,0,0,0}; 

  //randomly pick the components from uniform distribution in [0,1]
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> dis(0, 0.05);
  for(size_t i = 0; i < N; ++i) {
    for(size_t j = 0; j < hexdecp::N; ++j) {
      hps[i][j] = dis(gen);
    }
    for(size_t j = 0; j < octp::N; ++j) {
      ops[i][j] = dis(gen);
    }
    for(size_t j = 0; j < quadp::N; ++j) {
      qps[i][j] = dis(gen);
    }
    for(size_t j = 0; j < dp::N; ++j) {
      dps[i][j] = dis(gen);
    }
    for(size_t j = 0; j < mp::N; ++j) {
      mps[i][j] = dis(gen);
    }
  }
  //randomly pick the position from uniform distribution in [-5,5] in each dimension for each site
  random_device rdpos;
  mt19937 genpos(rdpos());
  uniform_real_distribution<> dispos(-10, 10);
  for(size_t i = 0; i < pos.size(); ++i) {
    for(size_t j = 0; j < 3; ++j) {
      pos[i][j] = dispos(genpos);
    }
  }
  
  //first check forces
  cout << setprecision(13) << scientific;
  cout << "# Difference and -log10(difference) between analytic and numerical forces (eps = "<<epsf<<") \n";
  //monopole
  for(size_t i = 0; i < N; ++i) {
    //check mp-mp
    for(size_t j = i+1; j < N; ++j) {
      const auto diff = chkEMPforce(mps[i], mps[j], pos[i]-pos[j], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check mp-dp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPforce(mps[i], dps[j], pos[i]-pos[j+N], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check mp-quadp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPforce(mps[i], qps[j], pos[i]-pos[j+2*N], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check mp-octp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPforce(mps[i], ops[j], pos[i]-pos[j+3*N], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check mp-hexdecp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPforce(mps[i], hps[j], pos[i]-pos[j+4*N], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
  }

  //dipole
  for(size_t i = 0; i < N; ++i) {
    //check dp-dp
    for(size_t j = i+1; j < N; ++j) {
      const auto diff = chkEMPforce(dps[i], dps[j], pos[i+N]-pos[j+N], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check dp-mp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPforce(dps[i], mps[j], pos[i+N]-pos[j], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check dp-quadp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPforce(dps[i], qps[j], pos[i+N]-pos[j+2*N], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check dp-octp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPforce(dps[i], ops[j], pos[i+N]-pos[j+3*N], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check dp-hexdecp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPforce(dps[i], hps[j], pos[i+N]-pos[j+4*N], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
  }

  //quadrupole
  for(size_t i = 0; i < N; ++i) {
    //check qp-qp
    for(size_t j = i+1; j < N; ++j) {
      const auto diff = chkEMPforce(qps[i], qps[j], pos[i+2*N]-pos[j+2*N], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check qp-mp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPforce(qps[i], mps[j], pos[i+2*N]-pos[j], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check qp-dp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPforce(qps[i], dps[j], pos[i+2*N]-pos[j+N], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check qp-octp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPforce(qps[i], ops[j], pos[i+2*N]-pos[j+3*N], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check qp-hexdecp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPforce(qps[i], hps[j], pos[i+2*N]-pos[j+4*N], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
  }

  //octupole
  for(size_t i = 0; i < N; ++i) {
    //check op-op
    for(size_t j = i+1; j < N; ++j) {
      const auto diff = chkEMPforce(ops[i], ops[j], pos[i+3*N]-pos[j+3*N], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check op-mp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPforce(ops[i], mps[j], pos[i+3*N]-pos[j], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check op-dp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPforce(ops[i], dps[j], pos[i+3*N]-pos[j+N], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check op-quaop
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPforce(ops[i], qps[j], pos[i+3*N]-pos[j+2*N], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check op-hexdecp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPforce(ops[i], hps[j], pos[i+3*N]-pos[j+4*N], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
  }

  //hexadecapole
  for(size_t i = 0; i < N; ++i) {
    //check hp-hp
    for(size_t j = i+1; j < N; ++j) {
      const auto diff = chkEMPforce(hps[i], hps[j], pos[i+4*N]-pos[j+4*N], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check hp-mp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPforce(hps[i], mps[j], pos[i+4*N]-pos[j], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check hp-dp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPforce(hps[i], dps[j], pos[i+4*N]-pos[j+N], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check hp-quadp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPforce(hps[i], qps[j], pos[i+4*N]-pos[j+2*N], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check hp-octp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPforce(hps[i], ops[j], pos[i+4*N]-pos[j+3*N], epsf);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
  }

  //then check torques
  cout << "# Difference and -log10(difference) between analytic and numerical torques (eps = "<<epst<<") \n";
  //monopole doesn't experience torque

  //dipole
  for(size_t i = 0; i < N; ++i) {
    //check dp-dp
    for(size_t j = i+1; j < N; ++j) {
      const auto diff = chkEMPtorqm(dps[i], dps[j], pos[i+N]-pos[j+N], quatm, epst);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check dp-mp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPtorqm(dps[i], mps[j], pos[i+N]-pos[j], quatm, epst);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check dp-quadp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPtorqm(dps[i], qps[j], pos[i+N]-pos[j+2*N], quatm, epst);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check dp-octp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPtorqm(dps[i], ops[j], pos[i+N]-pos[j+3*N], quatm, epst);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check dp-hexdecp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPtorqm(dps[i], hps[j], pos[i+N]-pos[j+4*N], quatm, epst);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
  }

  //quadrupole
  for(size_t i = 0; i < N; ++i) {
    //check qp-qp
    for(size_t j = i+1; j < N; ++j) {
      const auto diff = chkEMPtorqm(qps[i], qps[j], pos[i+2*N]-pos[j+2*N], quatm, epst);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check qp-mp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPtorqm(qps[i], mps[j], pos[i+2*N]-pos[j], quatm, epst);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check qp-dp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPtorqm(qps[i], dps[j], pos[i+2*N]-pos[j+N], quatm, epst);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check qp-octp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPtorqm(qps[i], ops[j], pos[i+2*N]-pos[j+3*N], quatm, epst);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check qp-hexdecp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPtorqm(qps[i], hps[j], pos[i+2*N]-pos[j+4*N], quatm, epst);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
  }

  //octupole
  for(size_t i = 0; i < N; ++i) {
    //check op-op
    for(size_t j = i+1; j < N; ++j) {
      const auto diff = chkEMPtorqm(ops[i], ops[j], pos[i+3*N]-pos[j+3*N], quatm, epst);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check op-mp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPtorqm(ops[i], mps[j], pos[i+3*N]-pos[j], quatm, epst);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check op-dp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPtorqm(ops[i], dps[j], pos[i+3*N]-pos[j+N], quatm, epst);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check op-quaop
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPtorqm(ops[i], qps[j], pos[i+3*N]-pos[j+2*N], quatm, epst);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check op-hexdecp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPtorqm(ops[i], hps[j], pos[i+3*N]-pos[j+4*N], quatm, epst);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
  }

  //hexadecapole
  for(size_t i = 0; i < N; ++i) {
    //check hp-hp
    for(size_t j = i+1; j < N; ++j) {
      const auto diff = chkEMPtorqm(hps[i], hps[j], pos[i+4*N]-pos[j+4*N], quatm, epst);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check hp-mp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPtorqm(hps[i], mps[j], pos[i+4*N]-pos[j], quatm, epst);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check hp-dp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPtorqm(hps[i], dps[j], pos[i+4*N]-pos[j+N], quatm, epst);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check hp-quadp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPtorqm(hps[i], qps[j], pos[i+4*N]-pos[j+2*N], quatm, epst);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
    //check hp-octp
    for(size_t j = 0; j < N; ++j) {
      const auto diff = chkEMPtorqm(hps[i], ops[j], pos[i+4*N]-pos[j+3*N], quatm, epst);
      const real rmsd = sqrt((diff*diff)/3);
      cout << rmsd << " " << -log10(rmsd) << endl;
    }
  }
}
