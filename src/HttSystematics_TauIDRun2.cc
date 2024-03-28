#include "CombineHarvester/TauIDSFMeasurement/interface/HttSystematics_TauIDRun2.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include <string>
#include <vector>

using namespace std;

namespace ch {

using ch::syst::SystMap;
using ch::syst::SystMapAsymm;
using ch::syst::era;
using ch::syst::channel;
using ch::syst::bin_id;
using ch::syst::process;
using ch::syst::bin;
using ch::JoinStr;

void AddTauIDRun2Systematics(CombineHarvester &cb, bool jetfakes, bool embedding, bool regional_jec, int era, bool tauES) {

  // ##########################################################################
  // Define groups of processes
  // ##########################################################################

  std::vector<std::string> mc_processes;
  // Signal processes
  if (embedding){
    //   std::vector<std::string> signals = {"EMB"};
      std::vector<std::string> signals = {"EMB_Pt20to25", "EMB_Pt25to30", "EMB_Pt30to35", "EMB_Pt35to40", "EMB_PtGt40", "EMB_DM0" , "EMB_DM1", "EMB_DM10_11", "EMB_Inclusive" };

      mc_processes =
                  {"TT", "TTL", "TTJ", "W", "ZJ", "ZL", "VV", "VVL", "VVJ", "ST"};
  }
  else {
      std::vector<std::string> signals = {"ZTT"};
      mc_processes =
                  {"ZTT", "TT", "TTT", "TTL", "TTJ", "W", "ZJ", "ZL", "VV", "VVT", "VVL", "VVJ", "ST"};
  }

  // ##########################################################################
  // Uncertainty: Lumi
  // References:
  // - "CMS Luminosity Measurements for the 2016 Data Taking Period"
  //   (PAS, https://cds.cern.ch/record/2257069)
  // - Recommendation twiki
  //    https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#LumiComb
  // Notes:
  // - FIXME: Adapt for fake factor and embedding
  // ##########################################################################

  float lumi_unc = 1.0;
  float lumi_unc_corr = 1.0;
  float lumi_unc_1718 = 1.0;
  if (era == 2016) {
      lumi_unc = 1.010;
      lumi_unc_corr = 1.006;
  } else if (era == 2017) {
      lumi_unc = 1.020;
      lumi_unc_corr = 1.009;
      lumi_unc_1718 = 1.006;
  } else if (era == 2018) {
      lumi_unc = 1.015;
      lumi_unc_corr = 1.020;
      lumi_unc_1718 = 1.002;
  }
  cb.cp()
      .channel({"et", "mt", "tt", "em", "mm"})
      .process(mc_processes)
      .AddSyst(cb, "lumi_13TeV_$ERA", "lnN", SystMap<>::init(lumi_unc));
  cb.cp()
      .channel({"et", "mt", "tt", "em", "mm"})
      .process(mc_processes)
      .AddSyst(cb, "lumi_13TeV_correlated", "lnN", SystMap<>::init(lumi_unc_corr));
  cb.cp()
      .channel({"et", "mt", "tt", "em", "mm"})
      .process(mc_processes)
      .AddSyst(cb, "lumi_13TeV_1718", "lnN", SystMap<>::init(lumi_unc_1718));

  // ##########################################################################
  // Uncertainty: Prefiring
  // References:
  // - "https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe"
  // Notes:
  // - FIXME: assumed as uncorrelated accross the years for now, what is the recommendation?
  // ##########################################################################
  if (era != 2018) {
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_prefiring", "shape", SystMap<>::init(1.00));
  }


  // ##########################################################################
  // Uncertainty: Trigger efficiency
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################
// TODO add trigger sf shifts
  cb.cp()
      .channel({"mt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_trigger_mt_$ERA", "shape", SystMap<>::init(1.00));

  // 100% uncorrelated for embedded
  cb.cp()
      .channel({"mt", "mm"})
    //   .process({"EMB", "MUEMB"})
    //   .process({"EMB", "MUEMB"})
      .process({"EMB_Pt20to25", "EMB_Pt25to30", "EMB_Pt30to35", "EMB_Pt35to40", "EMB_PtGt40", "EMB_DM0" , "EMB_DM1", "EMB_DM10_11", "EMB_Inclusive", "MUEMB"})
      .AddSyst(cb, "CMS_eff_trigger_emb_mt_$ERA", "shape", SystMap<>::init(1.00));

    std::cout << "Before Muon ID systematics" << std::endl;
    cb.cp().PrintAll();

//  Muon ID 
    cb.cp()
        .channel({"mt"})
        .process({"EMB_Pt20to25", "EMB_Pt25to30", "EMB_Pt30to35", "EMB_Pt35to40", "EMB_PtGt40", "EMB_DM0", "EMB_DM1", "EMB_DM10_11", "EMB_Inclusive", "MUEMB"})
        .AddSyst(cb, "CMS_eff_m_emb", "lnN", SystMap<>::init(1.02));

    cb.cp()
        .channel({"mm"})
        .process({"EMB_Pt20to25", "EMB_Pt25to30", "EMB_Pt30to35", "EMB_Pt35to40", "EMB_PtGt40", "EMB_DM0", "EMB_DM1", "EMB_DM10_11", "EMB_Inclusive", "MUEMB"})
        .AddSyst(cb, "CMS_eff_m_emb", "lnN", SystMap<>::init(1.04));
    cb.cp()
        .channel({"mt"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_eff_m_mc", "lnN", SystMap<>::init(1.02));
    cb.cp()
        .channel({"mm"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_eff_m_mc", "lnN", SystMap<>::init(1.04));

    std::cout << "After Muon ID systematics" << std::endl;

    cb.cp().PrintAll();


    // cb.cp()
    //     .channel({"mt"})
    //     .process({"EMB_Pt20to25", "EMB_Pt25to30", "EMB_Pt30to35", "EMB_Pt35to40", "EMB_PtGt40", "EMB_DM0" , "EMB_DM1", "EMB_DM10_11", "EMB_Inclusive", "MUEMB"})
    //     .AddSyst(cb, "CMS_eff_m_emb", "lnN", SystMap<>::init(0.98));



  // ##########################################################################
  // Uncertainty: b-tag and mistag efficiency
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

  // cb.cp()
  //     .channel({"mt"})
  //     .process(mc_processes)
  //     .AddSyst(cb, "CMS_htt_eff_b_$ERA", "shape", SystMap<>::init(1.00));

  // cb.cp()
  //     .channel({"mt"})
  //     .process(mc_processes)
  //     .AddSyst(cb, "CMS_htt_mistag_b_$ERA", "shape", SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: Tau energy scale
  // References:
  // Notes:
  // - Tau energy scale is splitted by decay mode.
  // - FIXME: References?
  // ##########################################################################

  //// correlated part between ERAs
  // MC uncorrelated uncertainty
// TODO set to correlated again once embedding is available (0.71)
  if (tauES) {
        cb.cp()
            .channel({"mt"})
            .process({{"ZTT", "TTT", "TTL", "VVT", "VVL"}})
            // .process({{"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes", "QCD"}}) // TODO: Check if this should be propagated to QCD
            .AddSyst(cb, "CMS_scale_t_1prong_$ERA", "shape", SystMap<>::init(1.00));

        cb.cp()
            .channel({"mt"})
            .process({"ZTT", "TTT", "TTL", "VVT", "VVL"})
            .AddSyst(cb, "CMS_scale_t_1prong1pizero_$ERA", "shape",
                    SystMap<>::init(1.00));

        cb.cp()
            .channel({"mt"})
            .process({"ZTT", "TTT", "TTL", "VVT", "VVL"})
            .AddSyst(cb, "CMS_scale_t_3prong_$ERA", "shape", SystMap<>::init(1.00));

        cb.cp()
            .channel({"mt"})
            .process({"ZTT", "TTT", "TTL", "VVT", "VVL"})
            .AddSyst(cb, "CMS_scale_t_3prong1pizero_$ERA", "shape", SystMap<>::init(1.00));
    }
    else{
        cb.cp()
        .channel({"mt"})
        .process({{"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes"}})
        // .process({{"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes", "QCD"}}) // TODO: Check if this should be propagated to QCD
        .AddSyst(cb, "CMS_scale_t_1prong_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"mt"})
        .process({"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes"})
        .AddSyst(cb, "CMS_scale_t_1prong1pizero_$ERA", "shape",
                SystMap<>::init(1.00));

    cb.cp()
        .channel({"mt"})
        .process({"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes"})
        .AddSyst(cb, "CMS_scale_t_3prong_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"mt"})
        .process({"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes"})
        .AddSyst(cb, "CMS_scale_t_3prong1pizero_$ERA", "shape", SystMap<>::init(1.00));
    }

  // Embedded uncorrelated uncertainty
// TODO check EMB ES
if (not tauES) {
  cb.cp()
      .channel({"mt"})
    //   .process({"EMB", "jetFakes", "QCDEMB"})
      .process({"EMB_Pt20to25", "EMB_Pt25to30", "EMB_Pt30to35", "EMB_Pt35to40", "EMB_PtGt40", "EMB_DM0" , "EMB_DM1", "EMB_DM10_11", "EMB_Inclusive", "jetFakes", "QCDEMB"})
      .AddSyst(cb, "CMS_scale_t_emb_1prong_$ERA", "shape", SystMap<>::init(0.71));

  cb.cp()
      .channel({"mt"})
    //   .process({"EMB", "jetFakes"})
      .process({"EMB_Pt20to25", "EMB_Pt25to30", "EMB_Pt30to35", "EMB_Pt35to40", "EMB_PtGt40", "EMB_DM0" , "EMB_DM1", "EMB_DM10_11", "EMB_Inclusive", "jetFakes"})
      .AddSyst(cb, "CMS_scale_t_emb_1prong1pizero_$ERA", "shape", SystMap<>::init(0.71));

  cb.cp()
      .channel({"mt"})
    //   .process({"EMB", "jetFakes"})
      .process({"EMB_Pt20to25", "EMB_Pt25to30", "EMB_Pt30to35", "EMB_Pt35to40", "EMB_PtGt40", "EMB_DM0" , "EMB_DM1", "EMB_DM10_11", "EMB_Inclusive", "jetFakes"})
      .AddSyst(cb, "CMS_scale_t_emb_3prong_$ERA", "shape", SystMap<>::init(0.71));
  cb.cp()
      .channel({"mt"})
    //   .process({"EMB", "jetFakes"})
      .process({"EMB_Pt20to25", "EMB_Pt25to30", "EMB_Pt30to35", "EMB_Pt35to40", "EMB_PtGt40", "EMB_DM0" , "EMB_DM1", "EMB_DM10_11", "EMB_Inclusive", "jetFakes"})
      .AddSyst(cb, "CMS_scale_t_emb_3prong1pizero_$ERA", "shape", SystMap<>::init(0.71));

//   MC + embedded correlated uncertainty

  cb.cp()
      .channel({"mt"})
    //   .process(JoinStr({{"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes"}, {"EMB"}}))
       .process(JoinStr({{"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes"}, {"EMB_Pt20to25", "EMB_Pt25to30", "EMB_Pt30to35", "EMB_Pt35to40", "EMB_PtGt40", "EMB_DM0" , "EMB_DM1", "EMB_DM10_11", "EMB_Inclusive"}}))
      // .process(JoinStr({{"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes", "QCD"}, {"EMB", "QCDEMB"}}))
      .AddSyst(cb, "CMS_scale_t_1prong_$ERA", "shape", SystMap<>::init(0.71));

  cb.cp()
      .channel({"mt"})
    //   .process(JoinStr({{"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes"}, {"EMB"}}))
      .process(JoinStr({{"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes"}, {"EMB_Pt20to25", "EMB_Pt25to30", "EMB_Pt30to35", "EMB_Pt35to40", "EMB_PtGt40", "EMB_DM0" , "EMB_DM1", "EMB_DM10_11", "EMB_Inclusive"}}))
      .AddSyst(cb, "CMS_scale_t_1prong1pizero_$ERA", "shape", SystMap<>::init(0.71));

  cb.cp()
      .channel({"mt"})
    //   .process(JoinStr({{"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes"}, {"EMB"}}))
      .process(JoinStr({{"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes"}, {"EMB_Pt20to25", "EMB_Pt25to30", "EMB_Pt30to35", "EMB_Pt35to40", "EMB_PtGt40", "EMB_DM0" , "EMB_DM1", "EMB_DM10_11", "EMB_Inclusive"}}))
      .AddSyst(cb, "CMS_scale_t_3prong_$ERA", "shape", SystMap<>::init(0.71));

  cb.cp()
      .channel({"mt"})
    //   .process(JoinStr({{"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes"}, {"EMB"}}))
      .process(JoinStr({{"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes"}, {"EMB_Pt20to25", "EMB_Pt25to30", "EMB_Pt30to35", "EMB_Pt35to40", "EMB_PtGt40", "EMB_DM0" , "EMB_DM1", "EMB_DM10_11", "EMB_Inclusive"}}))
      .AddSyst(cb, "CMS_scale_t_3prong1pizero_$ERA", "shape", SystMap<>::init(0.71));
}

  // ##########################################################################
  // Uncertainty: Jet energy scale
  // References:
  // - Talk in CMS Htt meeting by Daniel Winterbottom about regional JES splits:
  //   https://indico.cern.ch/event/740094/contributions/3055870/
  // Notes:
  // ##########################################################################

  // Regional JES
    // uncorrelated between eras
    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_Absolute_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_BBEC1_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_EC2_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_HF_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_RelativeSample_$ERA", "shape", SystMap<>::init(1.00));
    // correlated between eras
    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_Absolute", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_BBEC1", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_EC2", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_HF", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_FlavorQCD", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_RelativeBal", "shape", SystMap<>::init(1.00));

  // JER
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_res_j_$ERA", "shape", SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: MET energy scale and Recoil
  // References:
  // Notes:
  // - FIXME: Clustered vs unclustered MET? Inclusion of JES splitting?
  // - FIXME: References?
  // ##########################################################################

  cb.cp()
      .channel({"mt"})
      .process({"ZTT", "TT", "TTT", "TTL", "TTJ", "W", "ZJ", "ZL", "VV", "VVT", "VVL", "VVJ", "ST"})  //Z and W processes are only included due to the EWK fraction. Make sure that there is no contribution to the shift from the DY or Wjets samples.
      .AddSyst(cb, "CMS_scale_met_unclustered_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"ZTT", "ZL", "ZJ", "W"})
      .AddSyst(cb, "CMS_htt_boson_scale_met_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
     .process({"ZTT", "ZL", "ZJ", "W"})
      .AddSyst(cb, "CMS_htt_boson_res_met_$ERA", "shape", SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: Background normalizations
  // References:
  // Notes:
  // - FIXME: Remeasure QCD extrapolation factors for SS and ABCD methods?
  //          Current values are measured by KIT.
  // - FIXME: Adapt for fake factor and embedding
  // - FIXME: W uncertainties: Do we need lnN uncertainties based on the Ersatz
  //          study in Run1 (found in HIG-16043 uncertainty model)
  // - FIXME: References?
  // ##########################################################################

  // VV
  cb.cp()
      .channel({"mt", "mm"})
      .process({"VVT", "VVJ", "VVL", "VV"})
      .AddSyst(cb, "CMS_htt_vvXsec", "lnN", SystMap<>::init(1.06));

  // TT
  cb.cp()
      .channel({"mt", "mm"})
      .process({"TTT", "TTL", "TTJ", "TT"})
      .AddSyst(cb, "CMS_htt_tjXsec", "lnN", SystMap<>::init(1.055));

  cb.cp()
      .channel({"mt", "mm"})
      .process({"ST"})
      .AddSyst(cb, "CMS_htt_stXsec", "lnN", SystMap<>::init(1.055));

  // W
  cb.cp()
      .channel({"mt", "mm"})
      .process({"W"})
      .AddSyst(cb, "CMS_htt_wjXsec", "lnN", SystMap<>::init(1.05));

  // Z
  if (embedding) {
      cb.cp()
          .channel({"mt", "mm"})
          .process({"ZL", "ZJ"})
          .AddSyst(cb, "CMS_htt_zjXsec", "lnN", SystMap<>::init(1.04));
  }
  else {
      cb.cp()
          .channel({"mt"})
          .process({"ZTT", "ZL", "ZJ"})
          .AddSyst(cb, "CMS_htt_zjXsec", "rateParam", SystMap<>::init(1.00));
      cb.cp()
          .channel({"mm"})
          .process({"ZLL"})
          .AddSyst(cb, "CMS_htt_zjXsec", "rateParam", SystMap<>::init(1.00));
  }

  // QCD
  cb.cp()
      .channel({"mt", "mm"})
      .process({"QCD", "QCDEMB"})
      .AddSyst(cb, "CMS_ExtrapSSOS_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.30));

  // ##########################################################################
  // Uncertainty: Drell-Yan LO->NLO reweighting
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

  cb.cp()
      .channel({"mt"})
      .process({"ZTT", "ZL", "ZJ"})
      .AddSyst(cb, "CMS_htt_dyShape_$ERA", "shape", SystMap<>::init(0.10));

  // ##########################################################################
  // Uncertainty: TT shape reweighting
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

  cb.cp()
      .channel({"mt"})
      .process({"TTT", "TTL", "TTJ", "TT"})
      .AddSyst(cb, "CMS_htt_ttbarShape", "shape", SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: Electron/muon to tau fakes and ZL energy scale
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

  // ZL energy scale splitted by decay mode
  cb.cp()
      .channel({"mt"})
      .process({"ZL"})
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_$ERA", "shape",
               SystMap<>::init(1.00));

  // Muon fakes
  cb.cp()
      .channel({"mt"})
      .process({"ZL"})
      .AddSyst(cb, "CMS_fake_m_WH1_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"ZL"})
      .AddSyst(cb, "CMS_fake_m_WH2_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"ZL"})
      .AddSyst(cb, "CMS_fake_m_WH3_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"ZL"})
      .AddSyst(cb, "CMS_fake_m_WH4_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"ZL"})
      .AddSyst(cb, "CMS_fake_m_WH5_$ERA", "shape",
               SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: Jet to tau fakes
  // References:
  // Notes:
  // - FIXME: Adapt for fake factor and embedding
  // - FIXME: References?
  // ##########################################################################

  cb.cp()
      .channel({"mt"})
      .process({"W", "TTJ", "ZJ", "VVJ"})
      .AddSyst(cb, "CMS_htt_fake_j_$ERA", "shape", SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: Embedded events
  // References:
  // - https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTauEmbeddingSamples2016
  // Notes:
  // ##########################################################################

  // Embedded Normalization: No Lumi, Zjxsec information used, instead derived from data using dimuon selection efficiency
  cb.cp()
      .channel({"mt", "mm"})
    //   .process({"EMB", "MUEMB"})
      .process({"EMB_Pt20to25", "EMB_Pt25to30", "EMB_Pt30to35", "EMB_Pt35to40", "EMB_PtGt40", "EMB_DM0" , "EMB_DM1", "EMB_DM10_11", "EMB_Inclusive", "MUEMB"})
      .AddSyst(cb, "CMS_htt_doublemutrg_$ERA", "lnN", SystMap<>::init(1.04));

  // TTbar contamination in embedded events: 10% shape uncertainty of assumed ttbar->tautau event shape
  cb.cp()
    .channel({"mt"})
    // .process({"EMB"})
    .process({"EMB_Pt20to25", "EMB_Pt25to30", "EMB_Pt30to35", "EMB_Pt35to40", "EMB_PtGt40", "EMB_DM0" , "EMB_DM1", "EMB_DM10_11", "EMB_Inclusive"})
    .AddSyst(cb, "CMS_htt_emb_ttbar_$ERA", "shape", SystMap<>::init(1.00));

//   cb.cp()
//       .channel({"mt"})
//       .process({"EMB"})
//       .AddSyst(cb, "CMS_htt_emb_tauID_$ERA", "lnN", SystMap<>::init(1.05));

  // Uncertainty of hadronic tau track efficiency correction
  // uncorrelated between eras
  // TODO add hadronic tau track efficiency correction uncertainty
//   cb.cp()
//     .channel({"mt"})
//     .process({"EMB"})
//     .AddSyst(cb, "CMS_3ProngEff", "shape", SystMap<>::init(1.00));

//   cb.cp()
//     .channel({"mt"})
//     .process({"EMB"})
//     .AddSyst(cb, "CMS_1ProngPi0Eff", "shape", SystMap<>::init(1.00));


  // ##########################################################################
  // Uncertainty: Jet fakes
  // References:
  // - https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauJet2TauFakes
  // Notes:
  // - FIXME: add 2017 norm uncertainties, and properly correlate across years
  // ##########################################################################

  // QCD shape stat.
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_dr0_njet0_morphed_stat_", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_dr0_njet1_morphed_stat_", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_dr0_njet2_morphed_stat_", "shape", SystMap<>::init(1.00));


  // W shape stat.
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_dr0_njet0_morphed_stat_", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_dr0_njet1_morphed_stat_", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_dr0_njet2_morphed_stat_", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_dr1_njet0_morphed_stat_", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_dr1_njet1_morphed_stat_", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_dr1_njet2_morphed_stat_", "shape", SystMap<>::init(1.00));



  // TT shape stat.
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_tt_dr0_njet0_morphed_stat_", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_tt_dr0_njet1_morphed_stat_", "shape", SystMap<>::init(1.00));
  
  // MC subtraction uncertainty
  // uncorrelated between eras
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_mc_", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_mc_", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_frac_w_", "shape", SystMap<>::init(1.0));

      
  // Shape syst. of different contributions (QCD/W/tt)
  // uncorrelated between eras
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_mvis_", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_mvis_osss_", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_corr_qcd_mvis_", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_corr_qcd_mvis_osss_", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_muiso_", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_corr_qcd_muiso_", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_tau2_pt_0jet_", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_corr_qcd_tau2_pt_0jet_", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_tau2_pt_1jet_", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_corr_qcd_tau2_pt_1jet_", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_syst_", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_tt_syst_", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_tt_morphed_", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_corr_tt_syst_", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_tt_sf_", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_lepPt_", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_corr_w_lepPt_", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_mt_", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_corr_w_mt_", "shape", SystMap<>::init(1.0));

  //below: jetFakes norm uncertainties. Current values are for 2016, which are probably a good approx. for 2017. To be updated.

//   // Stat. norm (uncorrelated across years)
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_norm_stat_$CHANNEL_$BIN_$ERA", "lnN", SystMap<channel, bin_id>::init
// 	       ({"mt"}, {11},  1.04) //w
// 	       ({"mt"}, {12},  1.052) //ztt
// 	       ({"mt"}, {13},  1.051) //tt
// 	       ({"mt"}, {14},  1.047) //ss
// 	       ({"mt"}, {15},  1.04) //zll
// 	       ({"mt"}, {16},  1.059) //misc
// 	       ({"mt"}, {20},  1.052) //emb
// 	       ({"mt"}, {21},  1.047) //ff
// 	       ({"mt"}, {300}, 1.037) //incl
// 	       ({"et"}, {11},  1.066) //w
// 	       ({"et"}, {12},  1.095) //ztt
// 	       ({"et"}, {13},  1.083) //tt
// 	       ({"et"}, {14},  1.054) //ss
// 	       ({"et"}, {15},  1.095) //zll
// 	       ({"et"}, {16},  1.107) //misc
// 	       ({"et"}, {20},  1.095) //emb
// 	       ({"et"}, {21},  1.066) //ff
// 	       ({"et"}, {300}, 1.065) //incl
// 	       ({"tt"}, {12},  1.049) //ztt
// 	       ({"tt"}, {16},  1.028) //misc
// 	       ({"tt"}, {17},  1.041) //noniso
// 	       ({"tt"}, {20},  1.049) //emb
// 	       ({"tt"}, {21},  1.041) //ff
// 	       ({"tt"}, {300}, 1.041) //incl
// 	       );

 //  // Syst. norm: Bin-correlated
 //  // uncorrelated between eras
 //  cb.cp()
 //      .channel({"et", "mt", "tt"})
 //      .process({"jetFakes"})
 //      .AddSyst(cb, "CMS_ff_norm_syst_$CHANNEL_$ERA", "lnN", SystMap<channel, bin_id>::init
 //               ({"mt"}, {1},   1.049) //ggh
 //               ({"mt"}, {2},   1.041) //qqh
 //               ({"mt"}, {11},  1.038) //w
 //               ({"mt"}, {12},  1.069) //ztt
 //               ({"mt"}, {13},  1.037) //tt
 //               ({"mt"}, {14},  1.064) //ss
 //               ({"mt"}, {15},  1.048) //zll
 //               ({"mt"}, {16},  1.064) //misc
 //               ({"mt"}, {100}, 1.042) //incl
 //               ({"et"}, {1},   1.042) //ggh
 //               ({"et"}, {2},   1.040) //qqh
 //               ({"et"}, {11},  1.037) //w
 //               ({"et"}, {12},  1.062) //ztt
 //               ({"et"}, {13},  1.040) //tt
 //               ({"et"}, {14},  1.045) //ss
 //               ({"et"}, {15},  1.051) //zll
 //               ({"et"}, {16},  1.041) //misc
 //               ({"et"}, {100}, 1.042) //incl
 //               ({"tt"}, {1},   1.068) //ggh
 //               ({"tt"}, {2},   1.067) //qqh
 //               ({"tt"}, {12},  1.067) //ztt
 //               ({"tt"}, {16},  1.078) //misc
 //               ({"tt"}, {17},  1.070) //noniso
 //               ({"tt"}, {100}, 1.067) //incl
 //               );
 //  // correlated between eras
 //  cb.cp()
 //      .channel({"et", "mt", "tt"})
 //      .process({"jetFakes"})
 //      .AddSyst(cb, "CMS_ff_norm_syst_$CHANNEL", "lnN", SystMap<channel, bin_id>::init
 //               ({"mt"}, {1},   1.049) //ggh
 //               ({"mt"}, {2},   1.041) //qqh
 //               ({"mt"}, {11},  1.038) //w
 //               ({"mt"}, {12},  1.069) //ztt
 //               ({"mt"}, {13},  1.037) //tt
 //               ({"mt"}, {14},  1.064) //ss
 //               ({"mt"}, {15},  1.048) //zll
 //               ({"mt"}, {16},  1.064) //misc
 //               ({"mt"}, {100}, 1.042) //incl
 //               ({"et"}, {1},   1.042) //ggh
 //               ({"et"}, {2},   1.040) //qqh
 //               ({"et"}, {11},  1.037) //w
 //               ({"et"}, {12},  1.062) //ztt
 //               ({"et"}, {13},  1.040) //tt
 //               ({"et"}, {14},  1.045) //ss
 //               ({"et"}, {15},  1.051) //zll
 //               ({"et"}, {16},  1.041) //misc
 //               ({"et"}, {100}, 1.042) //incl
 //               ({"tt"}, {1},   1.068) //ggh
 //               ({"tt"}, {2},   1.067) //qqh
 //               ({"tt"}, {12},  1.067) //ztt
 //               ({"tt"}, {16},  1.078) //misc
 //               ({"tt"}, {17},  1.070) //noniso
 //               ({"tt"}, {100}, 1.067) //incl
 //               );

 //  // Syst. norm: Bin-dependent, correlated across years
 //  // uncorrelated between eras
 //  cb.cp()
 //      .channel({"et", "mt", "tt"})
 //      .process({"jetFakes"})
 //      .AddSyst(cb, "CMS_ff_sub_syst_$CHANNEL_$BIN_$ERA", "lnN", SystMap<channel, bin_id>::init
 //               ({"mt"}, {1},   1.028) //ggh
 //               ({"mt"}, {2},   1.028) //qqh
 //               ({"mt"}, {11},  1.018) //w
 //               ({"mt"}, {12},  1.032) //ztt
 //               ({"mt"}, {13},  1.021) //tt
 //               ({"mt"}, {14},  1.014) //ss
 //               ({"mt"}, {15},  1.028) //zll
 //               ({"mt"}, {16},  1.025) //misc
 //               ({"mt"}, {100}, 1.025) //incl
 //               ({"et"}, {1},   1.028) //ggh
 //               ({"et"}, {2},   1.025) //qqh
 //               ({"et"}, {11},  1.014) //w
 //               ({"et"}, {12},  1.028) //ztt
 //               ({"et"}, {13},  1.021) //tt
 //               ({"et"}, {14},  1.014) //ss
 //               ({"et"}, {15},  1.028) //zll
 //               ({"et"}, {16},  1.025) //misc
 //               ({"et"}, {100}, 1.025) //incl
 //               ({"tt"}, {1},   1.021) //ggh
 //               ({"tt"}, {2},   1.021) //qqh
 //               ({"tt"}, {12},  1.025) //ztt
 //               ({"tt"}, {16},  1.021) //misc
 //               ({"tt"}, {17},  1.014) //noniso
 //               ({"tt"}, {100}, 1.021) //incl
 //               );
 //  // correlated between eras
 //  cb.cp()
 //      .channel({"et", "mt", "tt"})
 //      .process({"jetFakes"})
 //      .AddSyst(cb, "CMS_ff_sub_syst_$CHANNEL_$BIN", "lnN", SystMap<channel, bin_id>::init
 //               ({"mt"}, {1},   1.028) //ggh
 //               ({"mt"}, {2},   1.028) //qqh
 //               ({"mt"}, {11},  1.018) //w
 //               ({"mt"}, {12},  1.032) //ztt
 //               ({"mt"}, {13},  1.021) //tt
 //               ({"mt"}, {14},  1.014) //ss
 //               ({"mt"}, {15},  1.028) //zll
 //               ({"mt"}, {16},  1.025) //misc
 //               ({"mt"}, {100}, 1.025) //incl
 //               ({"et"}, {1},   1.028) //ggh
 //               ({"et"}, {2},   1.025) //qqh
 //               ({"et"}, {11},  1.014) //w
 //               ({"et"}, {12},  1.028) //ztt
 //               ({"et"}, {13},  1.021) //tt
 //               ({"et"}, {14},  1.014) //ss
 //               ({"et"}, {15},  1.028) //zll
 //               ({"et"}, {16},  1.025) //misc
 //               ({"et"}, {100}, 1.025) //incl
 //               ({"tt"}, {1},   1.021) //ggh
 //               ({"tt"}, {2},   1.021) //qqh
 //               ({"tt"}, {12},  1.025) //ztt
 //               ({"tt"}, {16},  1.021) //misc
 //               ({"tt"}, {17},  1.014) //noniso
 //               ({"tt"}, {100}, 1.021) //incl
 //               );
}
} // namespace ch
