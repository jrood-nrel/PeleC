#ifndef TRANSPORT_PARAMS_H
#define TRANSPORT_PARAMS_H

#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include "TransportTypes.H"
#include "EOS.H"

namespace pele::physics::transport {

struct ConstTransport;
struct SimpleTransport;

template <typename EOSType, typename TransportType>
struct TransParm
{
};

template <typename EOSType, typename TransportType>
struct InitTransParm
{
  void operator()(TransParm<EOSType, TransportType> /*tparm*/) {}
};

template <typename EOSType>
struct TransParm<EOSType, ConstTransport>
{
  amrex::Real const_viscosity{0.0};
  amrex::Real const_bulk_viscosity{0.0};
  amrex::Real const_diffusivity{0.0};
  amrex::Real const_conductivity{0.0};
};

template <typename EOSType>
struct InitTransParm<EOSType, ConstTransport>
{
  void operator()(TransParm<EOSType, ConstTransport>* tparm)
  {
    amrex::ParmParse pp("transport");
    pp.query("const_viscosity", tparm->const_viscosity);
    pp.query("const_bulk_viscosity", tparm->const_bulk_viscosity);
    pp.query("const_conductivity", tparm->const_conductivity);
    pp.query("const_diffusivity", tparm->const_diffusivity);
    std::string units = "CGS";
    pp.query("units", units);
    if (units == "MKS") {
      // If input in MKS, convert to CGS
      tparm->const_viscosity *= amrex::Real(1.0e1);
      tparm->const_bulk_viscosity *= amrex::Real(1.0e1);
      tparm->const_diffusivity *= amrex::Real(1.0e1);
      tparm->const_conductivity *= amrex::Real(1.0e5);
    } else if (units != "CGS") {
      amrex::Abort(
        "Unable to initialize transport using unknown units " + units +
        ". Use CGS (def) or MKS");
    }
  }
};

template <typename EOSType>
struct TransParm<EOSType, SimpleTransport>
{
  amrex::GpuArray<amrex::Real, NUM_SPECIES> wt = {0.0};
  amrex::GpuArray<amrex::Real, NUM_SPECIES> iwt = {0.0};
  amrex::GpuArray<amrex::Real, NUM_SPECIES> eps = {0.0};
  amrex::GpuArray<amrex::Real, NUM_SPECIES> sig = {0.0};
  amrex::GpuArray<amrex::Real, NUM_SPECIES> dip = {0.0};
  amrex::GpuArray<amrex::Real, NUM_SPECIES> pol = {0.0};
  amrex::GpuArray<amrex::Real, NUM_SPECIES> zrot = {0.0};
  amrex::GpuArray<amrex::Real, NUM_SPECIES* NUM_FIT> fitmu = {0.0};
  amrex::GpuArray<amrex::Real, NUM_SPECIES* NUM_FIT> fitlam = {0.0};
  amrex::GpuArray<amrex::Real, NUM_SPECIES* NUM_SPECIES* NUM_FIT> fitdbin = {
    0.0};
  amrex::GpuArray<int, NUM_SPECIES> nlin = {0};
};

template <typename EOSType>
struct InitTransParm<EOSType, SimpleTransport>
{
  void operator()(TransParm<EOSType, SimpleTransport>* tparm)
  {
    egtransetWT(tparm->wt.data());
    egtransetEPS(tparm->eps.data());
    egtransetSIG(tparm->sig.data());
    egtransetDIP(tparm->dip.data());
    egtransetPOL(tparm->pol.data());
    egtransetZROT(tparm->zrot.data());
    egtransetNLIN(tparm->nlin.data());
    egtransetCOFETA(tparm->fitmu.data());
    egtransetCOFLAM(tparm->fitlam.data());
    egtransetCOFD(tparm->fitdbin.data());

    for (int i = 0; i < NUM_SPECIES; ++i) {
      tparm->iwt[i] = 1. / tparm->wt[i];
    }
  }
};

template <typename TransportType>
class TransportParams
{
public:
  TransportParams() = default;

  ~TransportParams() = default;

  void allocate()
  {
    InitTransParm<EosType, TransportType>()(&m_h_trans_parm);
    if (!m_device_allocated) {
      m_d_trans_parm =
        (TransParm<EosType, TransportType>*)amrex::The_Device_Arena()->alloc(
          sizeof(m_h_trans_parm));
      m_device_allocated = true;
      sync_to_device();
    }
  }

  void deallocate()
  {
    if (m_device_allocated) {
      amrex::The_Device_Arena()->free(m_d_trans_parm);
    }
  }

  void sync_to_device()
  {
    if (!m_device_allocated) {
      amrex::Abort("Device params not allocated yet");
    } else {
      amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, &m_h_trans_parm, &m_h_trans_parm + 1,
        m_d_trans_parm);
    }
  }

  TransParm<EosType, TransportType>& host_trans_parm()
  {
    return m_h_trans_parm;
  }
  const TransParm<EosType, TransportType>* device_trans_parm()
  {
    return m_d_trans_parm;
  }

private:
  TransParm<EosType, TransportType> m_h_trans_parm;
  TransParm<EosType, TransportType>* m_d_trans_parm{nullptr};
  bool m_device_allocated{false};
};
} // namespace pele::physics::transport
#endif
