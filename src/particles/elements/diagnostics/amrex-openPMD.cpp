/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */

/// #include "openPMD.H"
#include "ImpactXVersion.H"
#include "particles/ImpactXParticleContainer.H"

#include <ablastr/particles/IndexHandling.H>

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>

#include <utility>

#include "amrex-openPMD.H"

namespace impactx::diagnostics
{

  class StepMgr
{
public:
  StepMgr(int step, AMReXWithOpenPMD* owner)
    :m_Step(step),
     m_Owner(owner)
  {
    m_Owner = owner;
    m_Owner->m_UserHandler->m_Writer->SetStep(m_Step);
  }
  ~StepMgr()
  {
    m_Owner->m_UserHandler->m_Writer->CloseStep(m_Step);
  }

private:
  int m_Step;
  AMReXWithOpenPMD* m_Owner;
};


  AMReXWithOpenPMD::AMReXWithOpenPMD()
  {
    // warpx has multiple diags, each should maintain its own handler
    m_UserHandler = amrex::openpmd_api::InitUserHandler(m_Prefix);
  }

  void AMReXWithOpenPMD::SetWriter(amrex::openpmd_api::AMReX_openPMDWriter* w)
  {
    BL_ASSERT ( m_UserHandler != nullptr );
    BL_ASSERT ( w != nullptr );

    // so the openpmd filepath assigned from input file is still in use
    w->m_openPMDPrefix = m_UserHandler->m_Writer->m_openPMDPrefix;
    w->m_openPMDEncoding = m_UserHandler->m_Writer->m_openPMDEncoding;
    w->m_openPMDFileType = m_UserHandler->m_Writer->m_openPMDFileType;
    w->m_openPMDSeriesOptions = m_UserHandler->m_Writer->m_openPMDSeriesOptions;

    m_UserHandler->m_Writer.reset(w);
  }


  AMReXWithOpenPMD::~AMReXWithOpenPMD()
{
  amrex::Print()<<" openpmd_api::close handler "<<m_Prefix<<" \n";
  amrex::openpmd_api::CloseUserHandler(m_UserHandler);
}

bool AMReXWithOpenPMD::InitLocalHandler(const std::string& prefix)
{
  if (m_Prefix.compare(prefix) == 0)
    return false;

  amrex::Print()<<" openpmd_api::Init handler "<<m_Prefix<<" \n";
  m_Prefix = prefix;
  m_UserHandler = amrex::openpmd_api::InitUserHandler(prefix);
  return true;
}



  //
  // BeamMonitor, with I/O handled by amrex::openpmd-api
  //
    void BeamMonitor::finalize ()
    {
      if (m_plotWriter != NULL) {
	auto m_Writer = (AMReXWithOpenPMD*)(m_plotWriter);
	delete m_Writer;
      }
      amrex::Print()<<" => [check]: is things in BeamMonitor::finalize() addressed? \n";
      /*
        // close shared series alias
        if (m_series.has_value())
        {
            auto series = std::any_cast<io::Series>(m_series);
            series.close();
            m_series.reset();
        }

        // remove from unique series map
        if (m_unique_series.count(m_series_name) != 0u)
            m_unique_series.erase(m_series_name);
      */
    }

    BeamMonitor::BeamMonitor (std::string series_name, std::string backend, std::string encoding)
    {
      auto m_Writer =  new AMReXWithOpenPMD();
#ifdef ImpactX_USE_OPENPMD
        // encoding of iterations in the series
        openPMD::IterationEncoding series_encoding = openPMD::IterationEncoding::groupBased;
        if ( "v" == encoding )
            series_encoding = openPMD::IterationEncoding::variableBased;
        else if ( "g" == encoding )
            series_encoding = openPMD::IterationEncoding::groupBased;
        else if ( "f" == encoding )
            series_encoding = openPMD::IterationEncoding::fileBased;

	if ( m_Writer->InitLocalHandler(series_name) )
	  {
	    AMReX_impactxWriter* testWriter = new AMReX_impactxWriter(series_encoding);
	    m_Writer->SetWriter(testWriter);
	  }
#else
        amrex::AllPrint() << "Warning: openPMD output requested but not compiled for series=" << m_series_name << "\n";
#endif
	m_plotWriter = m_Writer;
    }

    void BeamMonitor::prepare (
        PinnedContainer & pc,
        RefPart const & ref_part,
        int step
    ) {
#ifdef ImpactX_USE_OPENPMD

	/* should be covered by amrex-openpmd-io
        // SoA: Real
        {
            std::vector<std::string> real_soa_names(RealSoA::names_s.size());
            std::copy(RealSoA::names_s.begin(), RealSoA::names_s.end(), real_soa_names.begin());
            for (auto real_idx = 0; real_idx < RealSoA::nattribs; real_idx++) {
                auto const component_name = real_soa_names.at(real_idx);
                getComponentRecord(component_name).resetDataset(d_fl);
            }
        }
        // SoA: Int
        static_assert(IntSoA::nattribs == 0); // not yet used
	*/
#else
        amrex::ignore_unused(pc, step);
#endif
    }

    void
    BeamMonitor::operator() (
        ImpactXParticleContainer & pc,
        int step
    )
    {
        // preparing to access reference particle data: RefPart
        RefPart & ref_part = pc.GetRefParticle();

        // pinned memory copy
        PinnedContainer pinned_pc = pc.make_alike<amrex::PinnedArenaAllocator>();
        pinned_pc.copyParticles(pc, true);  // no filtering

        // TODO: filtering
        /*
        using SrcData = WarpXParticleContainer::ParticleTileType::ConstParticleTileDataType;
        tmp.copyParticles(*pc,
                          [=] AMREX_GPU_HOST_DEVICE (const SrcData& src, int ip, const amrex::RandomEngine& engine)
                          {
                              const SuperParticleType& p = src.getSuperParticle(ip);
                              return random_filter(p, engine) * uniform_filter(p, engine)
                                     * parser_filter(p, engine) * geometry_filter(p, engine);
                          }, true);
        */

	//auto ixWriter = std::any_cast<AMReXWithOpenPMD>(m_Writer);
	auto m_Writer = (AMReXWithOpenPMD*)(m_plotWriter);
	//StepMgr sm(step, m_Writer->get());
	StepMgr sm(step, m_Writer);
	pinned_pc.CountParticles();

	//auto hi = pc.m_PtlCounter.m_Total;

      AMReX_impactxWriter* impactxWriter =  (AMReX_impactxWriter*) (m_Writer->m_UserHandler->m_Writer.get());
      amrex::Vector<std::string> real_names;
      amrex::Vector<std::string> int_names;
      amrex::Vector<int> int_flags;
      amrex::Vector<int> real_flags;
      impactxWriter->GetNames(real_names, int_names, int_flags, real_flags);
      m_Writer->m_UserHandler->m_Writer->DumpParticles(pinned_pc,
						       "beam",
						       real_flags,
						       int_flags,
						       real_names,
						       int_names,

						       [=] ([[maybe_unused]] auto& ppc, openPMD::ParticleSpecies& currSpecies, [[maybe_unused]]  unsigned long long localTotal)
						       {
							 //impactxWriter->SetConstantRefPart(currSpecies, localTotal, ref_part);
							      impactxWriter->SetConstantRefPart(currSpecies,  ref_part);
						       },
						       [=] (auto& pti, openPMD::ParticleSpecies& currSpecies, unsigned long long offset)
						       {
							 // use the default
							 impactxWriter->Save_impactx_PosID(pti, currSpecies, offset);
						       });


      // prepare element access
      //this->prepare(pinned_pc, ref_part, step);
    }

/*

    void
    BeamMonitor::operator() (
        PinnedContainer::ParIterType & pti,
        RefPart const & ref_part
    )
    {
#ifdef ImpactX_USE_OPENPMD
        int const currentLevel = pti.GetLevel();

        auto & offset = m_offset.at(currentLevel); // ...

        // preparing access to particle data: AoS
        auto& aos = pti.GetArrayOfStructs();

        // series & iteration
        auto series = std::any_cast<io::Series>(m_series);
        io::WriteIterations iterations = series.writeIterations();
        io::Iteration iteration = iterations[m_step];

        // writing
        io::ParticleSpecies beam = iteration.particles["beam"];

        auto const numParticleOnTile = pti.numParticles();
        uint64_t const numParticleOnTile64 = static_cast<uint64_t>( numParticleOnTile );

        // Do not call storeChunk() with zero-sized particle tiles:
        //   https://github.com/openPMD/openPMD-api/issues/1147
        //if (numParticleOnTile == 0) continue;

        auto const scalar = openPMD::RecordComponent::SCALAR;
        auto const getComponentRecord = [&beam](std::string const comp_name) {
            return detail::get_component_record(beam, comp_name);
        };

        // AoS: position and particle ID
        {
            using vs = std::vector<std::string>;
            vs const positionComponents{"x", "y", "t"}; // TODO: generalize
            for (auto currDim = 0; currDim < AMREX_SPACEDIM; currDim++) {
                std::shared_ptr<amrex::ParticleReal> const curr(
                    new amrex::ParticleReal[numParticleOnTile],
                    [](amrex::ParticleReal const *p) { delete[] p; }
                );
                for (auto i = 0; i < numParticleOnTile; i++) {
                    curr.get()[i] = aos[i].pos(currDim);
                }
                std::string const positionComponent = positionComponents[currDim];
                beam["position"][positionComponent].storeChunk(curr, {offset},
                                                               {numParticleOnTile64});
            }

            // save particle ID after converting it to a globally unique ID
            std::shared_ptr<uint64_t> const ids(
                new uint64_t[numParticleOnTile],
                [](uint64_t const *p) { delete[] p; }
            );
            for (auto i = 0; i < numParticleOnTile; i++) {
                ids.get()[i] = ablastr::particles::localIDtoGlobal(aos[i].id(), aos[i].cpu());
            }
            beam["id"][scalar].storeChunk(ids, {offset}, {numParticleOnTile64});
        }

        // SoA: everything else
        auto const& soa = pti.GetStructOfArrays();
        //   SoA floating point (ParticleReal) properties
        {
            std::vector<std::string> real_soa_names(RealSoA::names_s.size());
            std::copy(RealSoA::names_s.begin(), RealSoA::names_s.end(), real_soa_names.begin());

            for (auto real_idx=0; real_idx < RealSoA::nattribs; real_idx++) {
                auto const component_name = real_soa_names.at(real_idx);
                getComponentRecord(component_name).storeChunkRaw(
                soa.GetRealData(real_idx).data(), {offset}, {numParticleOnTile64});
            }
        }
        //   SoA integer (int) properties (not yet used)
        {
            static_assert(IntSoA::nattribs == 0); // not yet used
        }

        // TODO
        amrex::ignore_unused(ref_part);

        // needs to be higher for next pti; must be reset for next step via prepare
        offset += numParticleOnTile64;

        // TODO could be done once after all pti are processed
        // TODO at that point, we could also close the iteration/step
        series.flush();
#else
        amrex::ignore_unused(pti, ref_part);
#endif
    }
*/
} // namespace impactx::diagnostics
