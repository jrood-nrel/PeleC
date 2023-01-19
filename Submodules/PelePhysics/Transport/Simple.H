#ifndef SIMPLETRANSPORT_H
#define SIMPLETRANSPORT_H

#include "TransportParams.H"
#include "PhysicsConstants.H"

namespace pele::physics::transport {

template <typename EOSType>
struct NonIdealChungCorrections
{
  template <class... Args>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void operator()(Args... /*unused*/)
  {
  }
};

template <typename EOSType>
struct BinaryDiff
{
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void operator()(
    const amrex::Real* Xloc,
    const amrex::Real* Yloc,
    const amrex::Real* logT,
    const amrex::Real /*rholoc*/,
    const amrex::Real Tloc,
    amrex::Real* Ddiag,
    TransParm<EOSType, SimpleTransport> const* tparm)
  {
    const amrex::Real scale = Constants::PATM / (Constants::RU * Tloc);
    for (int i = 0; i < NUM_SPECIES; ++i) {
      amrex::Real term1 = 0.0;
      amrex::Real term2 = 0.0;
      for (int j = 0; j < NUM_SPECIES; ++j) {
        // cppcheck-suppress knownConditionTrueFalse
        if (i != j) {
          const int four_idx_ij = 4 * (i + NUM_SPECIES * j);
          const amrex::Real dbintemp =
            tparm->fitdbin[four_idx_ij] +
            tparm->fitdbin[1 + four_idx_ij] * logT[0] +
            tparm->fitdbin[2 + four_idx_ij] * logT[1] +
            tparm->fitdbin[3 + four_idx_ij] * logT[2];
          term1 += Yloc[j];
          term2 += Xloc[j] * std::exp(-dbintemp);
        }
      }
      Ddiag[i] = tparm->wt[i] * term1 / term2 * scale;
    }
  }
};

struct SimpleTransport
{
  using transport_type = SimpleTransport;

  static std::string identifier() { return "SimpleTransport"; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE static void comp_pure_bulk(
    amrex::Real Tloc,
    const amrex::Real* muloc,
    amrex::Real* xiloc,
    TransParm<EosType, transport_type> const* tparm)
  {
    amrex::Real cvk[NUM_SPECIES] = {0.0};
    amrex::Real cvkint[NUM_SPECIES] = {0.0};
    amrex::Real cvkrot[NUM_SPECIES] = {0.0};
    amrex::Real FofT[NUM_SPECIES] = {0.0};
    amrex::Real Fnorm[NUM_SPECIES] = {0.0};
    amrex::Real pi3_2 = std::pow(M_PI, 1.5);

    CKCVMS(&Tloc, cvk);

    for (int i = 0; i < NUM_SPECIES; ++i) {
      if (tparm->nlin[i] == 0) {
        cvkint[i] = 0.0;
        cvkrot[i] = 0.0;
      } else if (tparm->nlin[i] == 1) {
        cvkint[i] = cvk[i] * tparm->wt[i] / Constants::RU - 1.50;
        cvkrot[i] = 1.0;
      } else {
        cvkint[i] = cvk[i] * tparm->wt[i] / Constants::RU - 1.50;
        cvkrot[i] = 1.5;
      }
    }

    for (int i = 0; i < NUM_SPECIES; ++i) {
      amrex::Real epskoverTstd = tparm->eps[i] / 298.0;
      amrex::Real epskoverT = tparm->eps[i] / Tloc;
      Fnorm[i] = 1.0 + 0.50 * pi3_2 * std::sqrt(epskoverTstd) +
                 (2.0 + .50 * M_PI * M_PI) * epskoverTstd +
                 pi3_2 * std::sqrt(epskoverTstd) * epskoverTstd;

      FofT[i] = 1.0 + 0.50 * pi3_2 * std::sqrt(epskoverT) +
                (2.0 + .50 * M_PI * M_PI) * epskoverT +
                pi3_2 * std::sqrt(epskoverT) * epskoverT;
    }

    for (int i = 0; i < NUM_SPECIES; ++i) {
      if (tparm->nlin[i] == 0) {
        xiloc[i] = 0.0;
      } else {
        //   zrot/crot approximately zint / cint by assuming vibrational
        //   internal energy is small cvkrot is scaled by wk / Ru = mk / kb
        //   relative to standard specific cv
        xiloc[i] =
          0.250 * M_PI *
          (cvkint[i] / (cvkint[i] + 1.50) * cvkint[i] / (cvkint[i] + 1.50)) *
          tparm->zrot[i] / cvkrot[i] * Fnorm[i] / FofT[i] * muloc[i];
      }
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE static void transport(
    const bool wtr_get_xi,
    const bool wtr_get_mu,
    const bool wtr_get_lam,
    const bool wtr_get_Ddiag,
    const amrex::Real Tloc,
    const amrex::Real rholoc,
    amrex::Real* Yloc,
    amrex::Real* Ddiag,
    amrex::Real& mu,
    amrex::Real& xi,
    amrex::Real& lam,
    TransParm<EosType, transport_type> const* tparm)
  {
    amrex::Real trace = 1.e-15;
    amrex::Real Xloc[NUM_SPECIES] = {0.0};
    amrex::Real muloc[NUM_SPECIES] = {0.0};
    amrex::Real xiloc[NUM_SPECIES] = {0.0};
    amrex::Real logT[NUM_FIT - 1] = {0.0};

    logT[0] = std::log(Tloc);
    logT[1] = logT[0] * logT[0];
    logT[2] = logT[0] * logT[1];

    amrex::Real sum = 0.0;

    for (int i = 0; i < NUM_SPECIES; ++i) {
      sum += Yloc[i];
    }

    for (int i = 0; i < NUM_SPECIES; ++i) {
      Yloc[i] += trace * (sum / NUM_SPECIES - Yloc[i]);
    }

    amrex::Real wbar = 0.0;
    for (int i = 0; i < NUM_SPECIES; ++i) {
      wbar += Yloc[i] * tparm->iwt[i];
    }
    wbar = 1.0 / wbar;
    for (int i = 0; i < NUM_SPECIES; ++i) {
      Xloc[i] = Yloc[i] * wbar * tparm->iwt[i];
    }
    if (wtr_get_mu) {
      for (int i = 0; i < NUM_SPECIES; ++i) {
        muloc[i] = tparm->fitmu[4 * i] + tparm->fitmu[1 + 4 * i] * logT[0] +
                   tparm->fitmu[2 + 4 * i] * logT[1] +
                   tparm->fitmu[3 + 4 * i] * logT[2];
        muloc[i] = std::exp(muloc[i]);
      }

      mu = 0.0;

      for (int i = 0; i < NUM_SPECIES; ++i) {
        mu += Xloc[i] *
              (muloc[i] * muloc[i] * muloc[i] * muloc[i] * muloc[i] * muloc[i]);
      }
      mu = std::cbrt(std::sqrt(mu));

      //  assumption that we only get bulk viscosity if we are already getting
      //  shear viscosity

      if (wtr_get_xi) {
        comp_pure_bulk(Tloc, muloc, xiloc, tparm);
        xi = 0.0;
        for (int i = 0; i < NUM_SPECIES; ++i) {
          xi +=
            Xloc[i] * (std::sqrt(std::sqrt(xiloc[i] * xiloc[i] * xiloc[i])));
        }
        xi = std::cbrt(xi * xi * xi * xi);
      }
    }

    if (wtr_get_lam) {
      amrex::Real lamloc[NUM_SPECIES] = {0.0};
      for (int i = 0; i < NUM_SPECIES; ++i) {
        lamloc[i] = tparm->fitlam[4 * i] + tparm->fitlam[1 + 4 * i] * logT[0] +
                    tparm->fitlam[2 + 4 * i] * logT[1] +
                    tparm->fitlam[3 + 4 * i] * logT[2];
        lamloc[i] = std::exp(lamloc[i]);
      }

      lam = 0.0;
      for (int i = 0; i < NUM_SPECIES; ++i) {
        lam += Xloc[i] * std::sqrt(std::sqrt(lamloc[i]));
      }
      lam = lam * lam * lam * lam;
    }
    // Add nonideal corrections if necessary
    NonIdealChungCorrections<EosType>()(
      wtr_get_mu, wtr_get_lam, Tloc, Xloc, rholoc, wbar, mu, lam, tparm);

    if (wtr_get_Ddiag) {
      BinaryDiff<EosType>()(Xloc, Yloc, logT, rholoc, Tloc, Ddiag, tparm);
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE static void get_transport_coeffs(
    amrex::Box const& bx,
    amrex::Array4<const amrex::Real> const& Y_in,
    amrex::Array4<const amrex::Real> const& T_in,
    amrex::Array4<const amrex::Real> const& Rho_in,
    amrex::Array4<amrex::Real> const& D_out,
    amrex::Array4<amrex::Real> const& mu_out,
    amrex::Array4<amrex::Real> const& xi_out,
    amrex::Array4<amrex::Real> const& lam_out,
    TransParm<EosType, transport_type> const* tparm)
  {
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    const bool wtr_get_xi = true;
    const bool wtr_get_mu = true;
    const bool wtr_get_lam = true;
    const bool wtr_get_Ddiag = true;

    for (int k = lo.z; k <= hi.z; ++k) {
      for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {

          const amrex::Real T = T_in(i, j, k);
          const amrex::Real rho = Rho_in(i, j, k);
          amrex::Real massloc[NUM_SPECIES] = {0.0};
          for (int n = 0; n < NUM_SPECIES; ++n) {
            massloc[n] = Y_in(i, j, k, n);
          }

          amrex::Real muloc, xiloc, lamloc;
          amrex::Real Ddiag[NUM_SPECIES] = {0.0};
          transport(
            wtr_get_xi, wtr_get_mu, wtr_get_lam, wtr_get_Ddiag, T, rho, massloc,
            Ddiag, muloc, xiloc, lamloc, tparm);

          //   mu, xi and lambda are stored after D in the diffusion multifab
          for (int n = 0; n < NUM_SPECIES; ++n) {
            D_out(i, j, k, n) = Ddiag[n];
          }

          mu_out(i, j, k) = muloc;
          xi_out(i, j, k) = xiloc;
          lam_out(i, j, k) = lamloc;
        }
      }
    }
  }
  template <class... Args>
  AMREX_GPU_HOST_DEVICE explicit SimpleTransport(Args... /*unused*/)
  {
  }
};
} // namespace pele::physics::transport
#endif
