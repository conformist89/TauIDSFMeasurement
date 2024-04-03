#include "CombineHarvester/CombinePdfs/interface/MorphFunctions.h"
#include "CombineHarvester/CombineTools/interface/Algorithm.h"
#include "CombineHarvester/CombineTools/interface/AutoRebin.h"
#include "CombineHarvester/CombineTools/interface/BinByBin.h"
#include "CombineHarvester/CombineTools/interface/CardWriter.h"
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/Observation.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include "CombineHarvester/SMRun2Legacy/interface/BinomialBinByBin.h"
#include "CombineHarvester/TauIDSFMeasurement/interface/HttSystematics_TauIDRun2.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TF1.h"
#include "TH2.h"
#include "boost/algorithm/string/predicate.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/regex.hpp"
#include <cstdlib>
#include <iostream>
#include <map>
#include <math.h>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <iomanip> // Add this header for std::setprecision

using namespace std;
using boost::starts_with;
namespace po = boost::program_options;

int main(int argc, char **argv) {
  typedef vector<string> VString;
  typedef vector<pair<int, string>> Categories;
  using ch::JoinStr;
  using ch::syst::bin_id;

  // Define program options
  string output_folder = "sm_run2";
  string base_path = string(getenv("CMSSW_BASE")) +
                     "/src/CombineHarvester/TauIDSFMeasurement/shapes";
  string input_folder_mt = "Vienna/";
  string input_folder_mm = "Vienna/";
  string postfix = "-ML";
  bool regional_jec = true;
  bool auto_rebin = false;
  bool rebin_categories = true;
  bool manual_rebin_for_yields = false;
  bool real_data = false;
  bool jetfakes = true;
  bool embedding = false;
  bool classic_bbb = false;
  bool binomial_bbb = false;
  bool use_control_region = true;
  bool verbose = false;
  string categories = "pt_binned"; // "pt_binned", "ptdm_binned"
  int era = 2016;                  // 2016 or 2017
  po::variables_map vm;
  po::options_description config("configuration");
  config.add_options()("base_path",
                       po::value<string>(&base_path)->default_value(base_path))(
      "input_folder_mt",
      po::value<string>(&input_folder_mt)->default_value(input_folder_mt))(
      "input_folder_mm",
      po::value<string>(&input_folder_mm)->default_value(input_folder_mm))(
      "postfix", po::value<string>(&postfix)->default_value(postfix))(
      "auto_rebin", po::value<bool>(&auto_rebin)->default_value(auto_rebin))(
      "rebin_categories",
      po::value<bool>(&rebin_categories)->default_value(rebin_categories))(
      "manual_rebin_for_yields", po::value<bool>(&manual_rebin_for_yields)
                                     ->default_value(manual_rebin_for_yields))(
      "regional_jec",
      po::value<bool>(&regional_jec)->default_value(regional_jec))(
      "real_data", po::value<bool>(&real_data)->default_value(real_data))(
      "verbose", po::value<bool>(&verbose)->default_value(verbose))(
      "output_folder",
      po::value<string>(&output_folder)->default_value(output_folder))(
      "categories", po::value<string>(&categories)->default_value(categories))(
      "jetfakes", po::value<bool>(&jetfakes)->default_value(jetfakes))(
      "embedding", po::value<bool>(&embedding)->default_value(embedding))(
      "classic_bbb", po::value<bool>(&classic_bbb)->default_value(classic_bbb))(
      "binomial_bbb",
      po::value<bool>(&binomial_bbb)->default_value(binomial_bbb))(
      "use_control_region",
      po::value<bool>(&use_control_region)->default_value(use_control_region))(
      "era", po::value<int>(&era)->default_value(era));
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);

  // Define the location of the "auxiliaries" directory where we can
  // source the input files containing the datacard shapes
  std::map<string, string> input_dir;
  input_dir["mt"] = base_path + "/" + input_folder_mt + "/";
  input_dir["mm"] = base_path + "/" + input_folder_mm + "/";

  // Define channels
  VString chns;
  chns.push_back("mt");
  if (use_control_region)
    chns.push_back("mm");

  // Define background processes
  map<string, VString> bkg_procs;
  VString bkgs, bkgs_mm;
  bkgs = {"W", "QCD", "ZL", "ZJ", "TTT", "TTL", "TTJ", "VVJ", "VVT", "VVL"};

  // bkgs_mm = {"W", "TT", "VV"};

  if (embedding) {
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "TTT"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "VVT"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "QCD"), bkgs.end());
    bkgs = JoinStr({bkgs, {"QCD"}});
    bkgs_mm.erase(std::remove(bkgs_mm.begin(), bkgs_mm.end(), "ZLL"),
                  bkgs_mm.end());
    bkgs_mm.erase(std::remove(bkgs_mm.begin(), bkgs_mm.end(), "TT"),
                  bkgs_mm.end());
    bkgs_mm.erase(std::remove(bkgs_mm.begin(), bkgs_mm.end(), "VV"),
                  bkgs_mm.end());
  }
  if (jetfakes) {
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "QCD"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "QCD"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "W"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "VVJ"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "TTJ"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "ZJ"), bkgs.end());
    bkgs = JoinStr({bkgs, {"jetFakes"}});
  }

  bkgs_mm = {"W", "TTL", "VVL", "QCD", "MUEMB"};
  std::cout << "[INFO] Considerung the following processes:\n";
  std::cout << "For mt channel : \n";
  for (unsigned int i = 0; i < bkgs.size(); i++)
    std::cout << bkgs[i] << std::endl;
  if (use_control_region) {
    std::cout << "For mm channel : \n";
    for (unsigned int i = 0; i < bkgs_mm.size(); i++)
      std::cout << bkgs_mm[i] << std::endl;
  }
  bkg_procs["mt"] = bkgs;
  bkg_procs["mm"] = bkgs_mm;

  // vector<string> es_shifts;
  // for (float i = -2.5; i <= 2.5; i += 0.1) {
  //   es_shifts.push_back(std::to_string(i));
  // }
vector<string> es_shifts;
es_shifts.reserve(51);
es_shifts.emplace_back("-2.5");
es_shifts.emplace_back("-2.4");
es_shifts.emplace_back("-2.3");
es_shifts.emplace_back("-2.2");
es_shifts.emplace_back("-2.1");
es_shifts.emplace_back("-2.0");
es_shifts.emplace_back("-1.9");
es_shifts.emplace_back("-1.8");
es_shifts.emplace_back("-1.7");
es_shifts.emplace_back("-1.6");
es_shifts.emplace_back("-1.5");
es_shifts.emplace_back("-1.4");
es_shifts.emplace_back("-1.3");
es_shifts.emplace_back("-1.2");
es_shifts.emplace_back("-1.1");
es_shifts.emplace_back("-1.0");
es_shifts.emplace_back("-0.9");
es_shifts.emplace_back("-0.8");
es_shifts.emplace_back("-0.7");
es_shifts.emplace_back("-0.6");
es_shifts.emplace_back("-0.5");
es_shifts.emplace_back("-0.4");
es_shifts.emplace_back("-0.3");
es_shifts.emplace_back("-0.2");
es_shifts.emplace_back("-0.1");
// es_shifts.emplace_back("0.0");
es_shifts.emplace_back("0.1");
es_shifts.emplace_back("0.2");
es_shifts.emplace_back("0.3");
es_shifts.emplace_back("0.4");
es_shifts.emplace_back("0.5");
es_shifts.emplace_back("0.6");
es_shifts.emplace_back("0.7");
es_shifts.emplace_back("0.8");
es_shifts.emplace_back("0.9");
es_shifts.emplace_back("1.0");
es_shifts.emplace_back("1.1");
es_shifts.emplace_back("1.2");
es_shifts.emplace_back("1.3");
es_shifts.emplace_back("1.4");
es_shifts.emplace_back("1.5");
es_shifts.emplace_back("1.6");
es_shifts.emplace_back("1.7");
es_shifts.emplace_back("1.8");
es_shifts.emplace_back("1.9");
es_shifts.emplace_back("2.0");
es_shifts.emplace_back("2.1");
es_shifts.emplace_back("2.2");
es_shifts.emplace_back("2.3");
es_shifts.emplace_back("2.4");
es_shifts.emplace_back("2.5");

  
int es_shits_size = es_shifts.size(); 

std::vector<std::string> es_shapes_cats(es_shits_size * 3);
for(int i = 0; i < es_shits_size; ++i) {
  es_shapes_cats[i] = "DM0_"+es_shifts[i];
  es_shapes_cats[i + es_shits_size] = "DM1_"+es_shifts[i];
  es_shapes_cats[i + 2 * es_shits_size] = "DM10_11_"+es_shifts[i];
}

std::cout << "Elements of es_shapes_cats vector:" << std::endl;
for (const auto& element : es_shapes_cats) {
  std::cout << element << std::endl;
}


int es_shapes_cats_size =  es_shapes_cats.size();

// copilot please create a map of type <string, Categories> like this one map<string, Categories> es_shifts_cats
map<string, Categories> es_shifts_cats_map;
for (int i = 0; i < es_shapes_cats_size; ++i) {
  es_shifts_cats_map[es_shapes_cats[i]] = {
    {10 + i, es_shapes_cats[i]}
  };
}

vector<string>es_signals = {"es_shifts"};

// Define MSSM model-independent mass parameter MH
RooRealVar ES("ES", "ES", 0.0, -2.5, 2.5);
ES.setConstant(true);



std::vector<std::string> es_shapes_cats_dm0(es_shits_size);
for(int i = 0; i < es_shits_size; ++i) {
  es_shapes_cats_dm0[i] = "DM0_"+es_shifts[i];
}

std::vector<std::string> es_shapes_cats_dm1(es_shits_size);
for(int i = 0; i < es_shits_size; ++i) {
  es_shapes_cats_dm1[i] = "DM1_"+es_shifts[i];
}

std::vector<std::string> es_shapes_cats_dm10_11(es_shits_size);
for(int i = 0; i < es_shits_size; ++i) {
  es_shapes_cats_dm10_11[i] = "DM10_11_"+es_shifts[i];
}


std::cout << "Elements of es_shapes_cats_dm0 vector:" << std::endl;
for (const auto& element_dm0 : es_shapes_cats_dm0) {
  std::cout << element_dm0 << std::endl;
}

std::cout << "Successfully executed" << std::endl;
  map<string, Categories> cats;
  // TODO: Introduce maps for decay mode dependent splitting.
  vector<string> sig_procs;
  if (categories == "all") {
    cats["mt"] = {
         {1, "mt_Pt20to25"}, {2, "mt_Pt25to30"}, {3, "mt_Pt30to35"},
         {4, "mt_Pt35to40"}, {5, "mt_PtGt40"},   {6, "mt_Inclusive"},
         {7, "mt_DM0"},      {8, "mt_DM1"},      {9, "mt_DM10_11"},
    };
  

  }else if (categories == "Pt20to25"){
        cats["mt"] = {
         {1, "mt_Pt20to25"}
    };

    sig_procs = {"EMB_Pt20to25_0.0"};

    

  }else if (categories == "Pt25to30"){
        cats["mt"] = {
         {2, "mt_Pt25to30"}
    };
    sig_procs = {"EMB_Pt25to30_0.0"};

  
  }else if (categories == "Pt30to35"){
        cats["mt"] = {
         {3, "mt_Pt30to35"}
    };
    sig_procs = {"EMB_Pt30to35_0.0"};
  }

  else if (categories == "Pt35to40"){
        cats["mt"] = {
         {4, "mt_Pt35to40"}
    };
    sig_procs = {"EMB_Pt35to40_0.0"};
  
  }

    else if (categories == "PtGt40"){
        cats["mt"] = {
         {5, "mt_PtGt40"}
    };
    sig_procs = {"EMB_PtGt40_0.0"};
   
  }

    else if (categories == "Inclusive"){
        cats["mt"] = {
         {6, "mt_Inclusive"}
    };
    sig_procs = {"EMB_Inclusive_0.0"};
    
  }

    else if (categories == "DM0"){
        cats["mt"] = {
         {7, "mt_DM0"}
    };
    sig_procs = {"EMB_DM0"};
    // for (int i = 0; i < es_shits_size; ++i) {
    //   sig_procs.push_back("EMB_"+es_shapes_cats_dm0[i]);
    // }
  }

    else if (categories == "DM1"){
        cats["mt"] = {
         {8, "mt_DM1"}
    };
    sig_procs = {"EMB_DM1"};
    // for (int i = 0; i < es_shits_size; ++i) {
    //   sig_procs.push_back("EMB_"+es_shapes_cats_dm1[i]);
    // }
  }

    else if (categories == "DM10_11"){
        cats["mt"] = {
         {9, "mt_DM10_11"}
    };
    sig_procs = {"EMB_DM10_11"};
    // for (int i = 0; i < es_shits_size; ++i) {
    //   sig_procs.push_back("EMB_"+es_shapes_cats_dm10_11[i]);
    // }
  }
  
  else{
    throw std::runtime_error("Given categorization is not known.");}
  // else if (categories == "pt_binned") {
  //   cats["mt"] = {
  //       {1, "mt_Pt20to25"}, {2, "mt_Pt25to30"}, {3, "mt_Pt30to35"},
  //       {4, "mt_Pt35to40"}, {5, "mt_PtGt40"},
  //   };
  // }
  // // pt and dm binned scale factors
  // else if (categories == "ptdm_binned") {
  //   cats["mt"] = {
  //       {1, "mt_Pt20to25_DM0"},   {2, "mt_Pt20to25_DM1"},
  //       {3, "mt_Pt20to25_DM10"},  {4, "mt_Pt25to30_DM0"},
  //       {5, "mt_Pt25to30_DM1"},   {6, "mt_Pt25to30_DM10"},
  //       {7, "mt_Pt30to35_DM0"},   {8, "mt_Pt30to35_DM1"},
  //       {9, "mt_Pt30to35_DM10"},  {10, "mt_Pt35to40_DM0"},
  //       {11, "mt_Pt35to40_DM1"},  {12, "mt_Pt35to40_DM10"},
  //       {13, "mt_Pt40to50_DM0"},  {14, "mt_Pt40to50_DM1"},
  //       {15, "mt_Pt40to50_DM10"}, {16, "mt_Pt50to70_DM0"},
  //       {17, "mt_Pt50to70_DM1"},  {18, "mt_Pt50to70_DM10"},
  //       {19, "mt_PtGt70_DM0"},    {20, "mt_PtGt70_DM1"},
  //       {21, "mt_PtGt70_DM10"},
  //   };
  // } else if (categories == "dm_binned") {
  //   cats["mt"] = {
  //       {6, "mt_DM0"},
  //       {7, "mt_DM1"},
  //       {8, "mt_DM10_11"},
  //   };
  // } else
  //   throw std::runtime_error("Given categorization is not known.");
  cats["mm"] = {
      {100, "mm_control_region"},
  };




  vector<string> masses = {"125"};
  std::cout <<":::::::::::THIS IS MY SIGNAL::::::" << sig_procs[0] << std:: endl;
  // Create combine harverster object
  ch::CombineHarvester cb;

  // Add observations and processes
  std::string era_tag;
  if (era == 2016)
    era_tag = "Run2016";
  else if (era == 2017)
    era_tag = "Run2017";
  else if (era == 2018)
    era_tag = "Run2018";
  else
    std::runtime_error("Given era is not implemented.");



  for (auto chn : chns) {
    cb.AddObservations({"*"}, {"htt"}, {era_tag}, {chn}, cats[chn]);
    cb.AddProcesses({"*"}, {"htt"}, {era_tag}, {chn}, bkg_procs[chn], cats[chn],
                    false);

    std::cout << "#################### I have to add this process: " << sig_procs[0] << std::endl;

    
    // cb.AddProcesses(masses, {"htt"}, {era_tag}, {chn}, sig_procs, cats[chn],
    //                 true);
    cb.AddProcesses(es_shifts, {"htt"}, {era_tag}, {chn}, sig_procs, cats[chn],
                    true);
    // adding ES shifts
    // cb.AddProcesses(es_shifts, {"htt"}, {era_tag}, {chn}, es_signals, cats[chn], true);
  }

  // Add systematics
  ch::AddTauIDRun2Systematics(cb, jetfakes, embedding, regional_jec, era);

  // Extract shapes from input ROOT files
  for (string chn : chns) {
    cb.cp().channel({chn}).backgrounds().ExtractShapes(
        input_dir[chn] + "htt_" + chn + ".inputs-sm-" + era_tag + postfix +
            ".root",
        "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC");
    if (chn == "mt") {

      for (unsigned int i = 1; i <= cats[chn].size(); i++) {
        cb.cp()
            .channel({chn})
            // .bin_id({static_cast<int>(i)})
            .process(sig_procs)
            .ExtractShapes(input_dir[chn] + "htt_" + chn + ".inputs-sm-" +
                               era_tag + postfix + ".root",
                           "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC");
      }
    }
    std::cout << "[INFO] Extracted shapes for channel " << chn << std::endl;
    std::cout << "[INFO] Used input file  "
              << input_dir[chn] + "htt_" + chn + ".inputs-sm-" + era_tag +
                     postfix + ".root"
              << std::endl;

    std::cout << "[INFO] I have to add this process: " << sig_procs[0] << std::endl;
    
    // // TODO: Uncomment for check of mm fit to obtain start value for full fit.
    // if (chn == "mm") {
    //     cb.cp().channel({chn}).bin_id({static_cast<int>(100)}).process(sig_procs).ExtractShapes(
    //         input_dir[chn] + "htt_" + chn + ".inputs-sm-" + era_tag + postfix
    //         + ".root",
    //         "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC");
    // }
  }

  // Delete processes with 0 yield
  cb.FilterProcs([&](ch::Process *p) {

    std::cout << ch::Process::PrintHeader << *p << std::endl;
    std::cout << "Process rate: " << p->rate() << std::endl;
    bool null_yield = !(p->rate() > 0.0);
    if (null_yield) {
      std::cout << "[WARNING] Removing process with null yield: \n ";
      std::cout << ch::Process::PrintHeader << *p << "\n";
      
      cb.FilterSysts([&](ch::Systematic *s) {
        bool remove_syst = (MatchingProcess(*p, *s));
        return remove_syst;
      });
    }
    return null_yield;
  });

  // Delete systematics with 0 yield since these result in a bogus norm error in
  // combine
  cb.FilterSysts([&](ch::Systematic *s) {
    if (s->type() == "shape") {
      if (s->shape_u()->Integral() == 0.0) {
        // std::cout
        //     << "[WARNING] Removing systematic with null yield in up shift:"
        //     << std::endl;
        // std::cout << ch::Systematic::PrintHeader << *s << "\n";
        return true;
      }
      if (s->shape_d()->Integral() == 0.0) {
        // std::cout
        //     << "[WARNING] Removing systematic with null yield in down shift:"
        //     << std::endl;
        // std::cout << ch::Systematic::PrintHeader << *s << "\n";
        return true;
      }
    }
    return false;
  });

  int count_lnN = 0;
  int count_all = 0;
  cb.cp().ForEachSyst([&count_lnN, &count_all](ch::Systematic *s) {
    if (TString(s->name()).Contains("scale") ||
        TString(s->name()).Contains("CMS_htt_boson_reso_met")) {
      count_all++;
      double err_u = 0.0;
      double err_d = 0.0;
      int nbins = s->shape_u()->GetNbinsX();
      double yield_u = s->shape_u()->IntegralAndError(1, nbins, err_u);
      double yield_d = s->shape_d()->IntegralAndError(1, nbins, err_d);
      double value_u = s->value_u();
      double value_d = s->value_d();
      if (std::abs(value_u - 1.0) + std::abs(value_d - 1.0) <
          err_u / yield_u + err_d / yield_d) {
        count_lnN++;
        // std::cout << "[WARNING] Replacing systematic by lnN:" << std::endl;
        // std::cout << ch::Systematic::PrintHeader << *s << "\n";
        s->set_type("lnN");
        bool up_is_larger = (value_u > value_d);
        if (value_u < 1.0)
          value_u = 1.0 / value_u;
        if (value_d < 1.0)
          value_d = 1.0 / value_d;
        if (up_is_larger) {
          value_u = std::sqrt(value_u * value_d);
          value_d = 1.0 / value_u;
        } else {
          value_d = std::sqrt(value_u * value_d);
          value_u = 1.0 / value_d;
        }
        std::cout << "Former relative integral up shift: " << s->value_u()
                  << "; New relative integral up shift: " << value_u
                  << std::endl;
        std::cout << "Former relative integral down shift: " << s->value_d()
                  << "; New relative integral down shift: " << value_d
                  << std::endl;
        s->set_value_u(value_u);
        s->set_value_d(value_d);
      }
    }
  });
  std::cout << "[WARNING] Turned " << count_lnN << " of " << count_all
            << " checked systematics into lnN:" << std::endl;

  // systematics disentanglement among different pT and DM bins            

  // for(int i=1; i<=5; i++){

  //   cb.cp().bin_id({i}).channel({"mt"}).RenameSystematic(cb,"CMS_ExtrapSSOS_mt_Run2018", "CMS_ExtrapSSOS_mt_mt_"+std::to_string(i)+"Run2018");

  // }

  // for(int i=7; i<=9; i++){

  //   cb.cp().bin_id({i}).channel({"mt"}).RenameSystematic(cb,"CMS_ExtrapSSOS_mt_Run2018", "CMS_ExtrapSSOS_mt_mt_"+std::to_string(i)+"Run2018");

  // }

  //   for(int i=1; i<=5; i++){

  //   cb.cp().bin_id({i}).channel({"mt"}).RenameSystematic(cb,"CMS_scale_t_emb_1prong1pizero_Run2018", "CMS_scale_t_emb_1prong1pizero_mt_"+std::to_string(i)+"Run2018");

  // }

  // Replacing observation with the sum of the backgrounds (Asimov data)
  // useful to be able to check this, so don't do the replacement
  // for these
  if (!real_data) {
    for (auto b : cb.cp().bin_set()) {
      std::cout << "[INFO] Replacing data with asimov in bin " << b << "\n";
      cb.cp().bin({b}).ForEachObs([&](ch::Observation *obs) {
        obs->set_shape(cb.cp().bin({b}).backgrounds().GetShape() +
                           cb.cp().bin({b}).signals().GetShape(),
                       true);
      });
    }
  }

  // At this point we can fix the negative bins
  std::cout << "[INFO] Fixing negative bins.\n";
  cb.ForEachProc([](ch::Process *p) {
    if (ch::HasNegativeBins(p->shape())) {
      auto newhist = p->ClonedShape();
      ch::ZeroNegativeBins(newhist.get());
      p->set_shape(std::move(newhist), false);
    }
  });

  cb.ForEachSyst([](ch::Systematic *s) {
    if (s->type().find("shape") == std::string::npos)
      return;
    if (ch::HasNegativeBins(s->shape_u()) ||
        ch::HasNegativeBins(s->shape_d())) {
      auto newhist_u = s->ClonedShapeU();
      auto newhist_d = s->ClonedShapeD();
      ch::ZeroNegativeBins(newhist_u.get());
      ch::ZeroNegativeBins(newhist_d.get());
      s->set_shapes(std::move(newhist_u), std::move(newhist_d), nullptr);
    }
  });

  // Perform auto-rebinning
  if (auto_rebin) {
    const auto threshold = 1.0;
    const auto tolerance = 1e-4;
    for (auto b : cb.cp().bin_set()) {
      std::cout << "[INFO] Rebin bin " << b << "\n";
      // Get shape of this category with sum of backgrounds
      auto shape = cb.cp().bin({b}).backgrounds().GetShape();
      // Push back last bin edge
      vector<double> binning;
      const auto num_bins = shape.GetNbinsX();
      binning.push_back(shape.GetBinLowEdge(num_bins + 1));
      // Now, go backwards through bins (from right to left) and merge a bin if
      // the background yield is below a given threshold.
      auto offset = shape.GetBinLowEdge(1);
      auto width = 1.0 - offset;
      auto c = 0.0;
      for (auto i = num_bins; i > 0; i--) {
        // Determine whether this is a boundary of an unrolled category
        // if it's a multiple of the width between minimum NN score and 1.0.
        auto low_edge = shape.GetBinLowEdge(i);
        auto is_boundary =
            fabs(fmod(low_edge - offset, width)) < tolerance ? true : false;
        c += shape.GetBinContent(i);
        if (is_boundary) { // If the lower edge is a boundary, set a bin edge.
          if (c <= threshold &&
              !(fabs(fmod(binning[0] - offset, width)) < tolerance ||
                fabs(fmod(binning[0] - offset, width)) - width <
                    tolerance)) { // Special case: If this bin is at a boundary
                                  // but it is below the threshold and the bin
                                  // above is not again a boundary, merge to the
                                  // right.
            binning.erase(binning.begin());
          }
          binning.insert(binning.begin(), low_edge);
          c = 0.0;
        } else { // If this is not a boundary, check whether the content is
                 // above the threshold.
          if (c > threshold) { // Set lower edge if the bin content is above the
                               // threshold.
            binning.insert(binning.begin(), low_edge);
            c = 0.0;
          }
        }
      }
      cb.cp().bin({b}).VariableRebin(binning);
    }
    // blind subcategories with to little events
    for (auto b : cb.cp().bin_set()) {
      // Get shape of this category with sum of backgrounds
      auto shape = cb.cp().bin({b}).backgrounds().GetShape();
      const auto num_bins = shape.GetNbinsX();
      for (auto i = num_bins; i > 0; i--) {
        if (shape.GetBinContent(i) < threshold) {
          std::cout << "[INFO] Blind bin " << i << " in " << b
                    << " due to insufficient population!"
                    << "\n";
          cb.cp().bin({b}).ForEachProc([i](ch::Process *p) {
            auto newhist = p->ClonedShape();
            newhist->SetBinContent(i, 0.0);
            newhist->SetBinError(i, 0.0);
            p->set_shape(std::move(newhist), false);
          });
          cb.cp().bin({b}).ForEachObs([i](ch::Observation *p) {
            auto newhist = p->ClonedShape();
            newhist->SetBinContent(i, 0.0);
            newhist->SetBinError(i, 0.0);
            p->set_shape(std::move(newhist), false);
          });
          cb.cp().bin({b}).ForEachSyst([i](ch::Systematic *s) {
            if (s->type().find("shape") == std::string::npos)
              return;
            auto newhist_u = s->ClonedShapeU();
            auto newhist_d = s->ClonedShapeD();
            newhist_u->SetBinContent(i, 0.0);
            newhist_u->SetBinError(i, 0.0);
            newhist_d->SetBinContent(i, 0.0);
            newhist_d->SetBinError(i, 0.0);
            s->set_shapes(std::move(newhist_u), std::move(newhist_d), nullptr);
          });
        }
      }
    }
  }

  if (manual_rebin_for_yields) {
    for (auto b : cb.cp().bin_set()) {
      std::cout << "Rebinning by hand for bin: " << b << std::endl;
      cb.cp().bin({b}).VariableRebin({0.0, 1.0});
    }
  }

  // Merge bins and set bin-by-bin uncertainties if no autoMCStats is used.
  if (classic_bbb) {
    auto bbb = ch::BinByBinFactory()
                   .SetAddThreshold(0.0)
                   .SetMergeThreshold(0.5)
                   .SetFixNorm(false);
    bbb.MergeBinErrors(cb.cp().backgrounds());
    bbb.AddBinByBin(cb.cp().backgrounds(), cb);
  }
  if (binomial_bbb) {
    auto bbb =
        ch::BinomialBinByBinFactory()
            .SetPattern(
                "CMS_$ANALYSIS_$CHANNEL_$BIN_$ERA_$PROCESS_binomial_bin_$#")
            .SetBinomialP(0.022)
            .SetBinomialN(1000.0)
            .SetFixNorm(false);
    bbb.AddBinomialBinByBin(cb.cp().channel({"em"}).process({"EMB"}), cb);
  }

  // This function modifies every entry to have a standardised bin name of
  // the form: {analysis}_{channel}_{bin_id}_{era}
  ch::SetStandardBinNames(cb, "$ANALYSIS_$CHANNEL_$BINID_$ERA");

  cb.SetAutoMCStats(cb, 0.);

  // Write out datacards. Naming convention important for rest of workflow. We
  // make one directory per chn-cat, one per chn and cmb. In this code we only
  // store the individual datacards for each directory to be combined later.
  string output_prefix = "output/";
  ch::CardWriter writer(output_prefix + output_folder + "/$TAG/$BIN.txt",
                        output_prefix + output_folder + "/$TAG/htt_input_" +
                            era_tag +categories + ".root");

  // We're not using mass as an identifier - which we need to tell the
  // CardWriter
  // otherwise it will see "*" as the mass value for every object and skip it
  writer.SetWildcardMasses({});

  // Set verbosity
  if (verbose)
    writer.SetVerbosity(3);

  // Write datacards combined and per channel
  writer.WriteCards("cmb", cb);

  for (auto chn : chns) {
    if (chn == std::string("mm")) {
      continue;
    }
    if (use_control_region) {
      writer.WriteCards(chn, cb.cp().channel({chn, "mm"}));
    } else {
      writer.WriteCards(chn, cb.cp().channel({chn}));
    }
    // per-category
    for (auto cat : cats[chn]) {
      if (use_control_region) {
        writer.WriteCards(
            "htt_" + cat.second,
            cb.cp()
                .channel({chn, "mm"})
                .bin_id({cat.first, 100})); //.attr({cat.second,"control"},
                                            //"cat"));
      } else {
        writer.WriteCards("htt_" + cat.second,
                          cb.cp().channel({chn}).bin_id({cat.first, 100}));
      }
    }
  }

  // if (verbose)
  cb.PrintAll();

  cout << "[INFO] Done producing datacards.\n";
}
