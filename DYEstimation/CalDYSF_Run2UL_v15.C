#include <iostream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

using namespace std;

struct EvtErr {
    double Evt;
    double Err;
};

// SF result for one step
struct StepSF {
    TString label;
    double SF_ee,  SF_ee_e;
    double SF_mm,  SF_mm_e;
    double SF_em,  SF_em_e;
    double NMC_ee, NMC_ee_e;
    double NMC_mm, NMC_mm_e;
    double NMC_em, NMC_em_e;
    double NDY_ee, NDY_ee_e;
    double NDY_mm, NDY_mm_e;
    double NDY_em, NDY_em_e;
};

// Per-period configuration
struct PeriodConfig {
    TString period;   // e.g. "UL2018"
    TString base;     // base directory
    vector<TString> DY_names;
    vector<TString> Data_MuMu;
    vector<TString> Data_ElEl;
    vector<TString> Data_MuEl;
};

double CalcError(TH1D* h, int min_, int max_) {
    double err = 0.0;
    for (int i = min_; i <= max_; ++i) err += pow(h->GetBinError(i), 2);
    return sqrt(err);
}

struct EvtErr GetZmassIn(TH1D* h) {
    int min_ = h->GetXaxis()->FindBin(76.0);
    int max_ = h->GetXaxis()->FindBin(105.9999);
    EvtErr r;
    r.Evt = h->Integral(min_, max_);
    r.Err = CalcError(h, min_, max_);
    return r;
}

struct EvtErr GetZmassOut(TH1D* h) {
    int min_ = h->GetXaxis()->FindBin(76.0);
    int max_ = h->GetXaxis()->FindBin(105.9999);
    EvtErr r;
    r.Evt = h->Integral() - h->Integral(min_, max_);
    double err = 0.0;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        if (i >= min_ && i <= max_) continue;
        err += pow(h->GetBinError(i), 2);
    }
    r.Err = sqrt(err);
    return r;
}

double UncR(double x1, double x1_unc, double x2, double x2_unc) {
    double y = x1 / x2;
    return sqrt(pow(y/x1,2)*pow(x1_unc,2) + pow(y/x2,2)*pow(x2_unc,2));
}

double UncKll(double x1, double x1_unc, double x2, double x2_unc) {
    return sqrt(0.25*(1.0/(x1*x2))*pow(x1_unc,2)
              + 0.25*(x1/pow(x2,3))*pow(x2_unc,2));
}

double UncNoutD(double R,      double R_unc,
                double Nin_ll, double Nin_ll_unc,
                double Nin_em, double Nin_em_unc,
                double k,      double k_unc) {
    return sqrt(pow((Nin_ll - 0.5*Nin_em*k), 2) * pow(R_unc,     2)
              + pow(R,                        2) * pow(Nin_ll_unc, 2)
              + pow(-0.5*R*k,                 2) * pow(Nin_em_unc, 2)
              + pow(-0.5*R*Nin_em,            2) * pow(k_unc,      2));
}

double UncSFem(double SF_ee, double SF_ee_unc, double SF_mm, double SF_mm_unc) {
    return sqrt(0.25*(SF_mm/SF_ee)*pow(SF_ee_unc,2)
              + 0.25*(SF_ee/SF_mm)*pow(SF_mm_unc,2));
}

double UncFinal(double N, double N_unc, double SF, double SF_unc) {
    return sqrt(pow(SF,2)*pow(N_unc,2) + pow(N,2)*pow(SF_unc,2));
}

TH1D* SumHistograms(const vector<TString>& filepaths, const TString& hname) {
    TH1D* h_sum = nullptr;
    for (const auto& path : filepaths) {
        TFile* f = TFile::Open(path);
        if (!f || f->IsZombie()) {
            cerr << "[WARNING] Cannot open: " << path << endl;
            continue;
        }
        TH1D* h = (TH1D*)f->Get(hname);
        if (!h) {
            cerr << "[WARNING] Histogram not found: " << hname << " in " << path << endl;
            f->Close();
            continue;
        }
        if (!h_sum) {
            h_sum = (TH1D*)h->Clone(hname + "_sum");
            h_sum->SetDirectory(0);
        } else {
            h_sum->Add(h);
        }
        f->Close();
    }
    return h_sum;
}

void PrintTable(const TString& tag,
                double N_inMC_ee,  double N_inMC_ee_e,
                double N_outMC_ee, double N_outMC_ee_e,
                double R_ee,       double R_ee_e,
                double k_ee,       double k_ee_e,
                double N_inD_ee,   double N_inD_ee_e,
                double N_out_ee,   double N_out_ee_e,
                double SF_ee,      double SF_ee_e,
                double N_inMC_mm,  double N_inMC_mm_e,
                double N_outMC_mm, double N_outMC_mm_e,
                double R_mm,       double R_mm_e,
                double k_mm,       double k_mm_e,
                double N_inD_mm,   double N_inD_mm_e,
                double N_out_mm,   double N_out_mm_e,
                double SF_mm,      double SF_mm_e,
                double N_inD_em,   double N_inD_em_e,
                double SF_em,      double SF_em_e,
                double N_MC_em,    double N_MC_em_e,
                double N_DY_em,    double N_DY_em_e) {
    cout << "\n";
    cout << "  DY estimation [" << tag << "] - >= 1 b-tag" << endl;
    cout << "  " << string(72, '-') << endl;
    cout << Form("  %-22s  %20s  %20s  %20s",
                 "Level: 1btag", "e+e-", "mu+mu-", "e+-mu-+") << endl;
    cout << "  " << string(72, '-') << endl;
    cout << Form("  %-22s  %20s  %20s",
                 "N_in (MC)",
                 Form("%.0f +- %.0f", N_inMC_ee, N_inMC_ee_e),
                 Form("%.0f +- %.0f", N_inMC_mm, N_inMC_mm_e)) << endl;
    cout << Form("  %-22s  %20s  %20s",
                 "N_out (MC)",
                 Form("%.0f +- %.0f", N_outMC_ee, N_outMC_ee_e),
                 Form("%.0f +- %.0f", N_outMC_mm, N_outMC_mm_e)) << endl;
    cout << Form("  %-22s  %20s  %20s",
                 "R_out/in (MC)",
                 Form("%.3f +- %.3f", R_ee, R_ee_e),
                 Form("%.3f +- %.3f", R_mm, R_mm_e)) << endl;
    cout << Form("  %-22s  %20s  %20s",
                 "k_ll",
                 Form("%.3f +- %.3f", k_ee, k_ee_e),
                 Form("%.3f +- %.3f", k_mm, k_mm_e)) << endl;
    cout << Form("  %-22s  %20s  %20s  %20s",
                 "N_in (Data)",
                 Form("%.0f +- %.0f", N_inD_ee, N_inD_ee_e),
                 Form("%.0f +- %.0f", N_inD_mm, N_inD_mm_e),
                 Form("%.0f +- %.0f", N_inD_em, N_inD_em_e)) << endl;
    cout << Form("  %-22s  %20s  %20s",
                 "N_out (est)",
                 Form("%.1f +- %.1f", N_out_ee, N_out_ee_e),
                 Form("%.1f +- %.1f", N_out_mm, N_out_mm_e)) << endl;
    cout << Form("  %-22s  %20s  %20s  %20s",
                 "SF (Data/MC)",
                 Form("%.3f +- %.3f", SF_ee, SF_ee_e),
                 Form("%.3f +- %.3f", SF_mm, SF_mm_e),
                 Form("%.3f +- %.3f", SF_em, SF_em_e)) << endl;
    cout << Form("  %-22s  %20s  %20s  %20s",
                 "Drell-Yan (MC)", "", "",
                 Form("%.1f +- %.1f", N_MC_em, N_MC_em_e)) << endl;
    cout << Form("  %-22s  %20s  %20s  %20s",
                 "Drell-Yan (Data)", "", "",
                 Form("%.1f +- %.1f", N_DY_em, N_DY_em_e)) << endl;
    cout << "  " << string(72, '-') << endl;
}

void DrawPlot(TH1D* h_data, TH1D* h_dy50, TH1D* h_dy10,
              const TString& xlabel, const TString& outname) {
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas(outname, outname, 800, 600);
    c->cd();

    h_data->GetXaxis()->SetRangeUser(0, 300);
    h_data->GetXaxis()->SetTitle(xlabel);
    h_data->GetYaxis()->SetTitle("Entries");
    h_data->SetLineColor(kBlack);
    h_data->SetMarkerStyle(20);
    h_data->Draw("E");

    h_dy50->GetXaxis()->SetRangeUser(0, 300);
    h_dy50->SetLineColor(kAzure - 4);
    h_dy50->SetFillColor(kAzure - 4);
    h_dy50->Draw("HIST same");

    if (h_dy10) {
        h_dy10->GetXaxis()->SetRangeUser(0, 300);
        h_dy10->SetLineColor(kRed - 4);
        h_dy10->SetFillColor(kRed - 4);
        h_dy10->Draw("HIST same");
    }

    h_data->Draw("E same");

    TLegend* leg = new TLegend(0.65, 0.65, 0.92, 0.90);
    leg->AddEntry(h_data, "Data",                  "lpe");
    leg->AddEntry(h_dy50, "DY M>50 (amcatnlo)",    "f");
    if (h_dy10)
    leg->AddEntry(h_dy10, "DY 10<M<50 (amcatnlo)", "f");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);
    leg->SetTextFont(42);
    leg->Draw();

    c->SetLogy();
    c->SaveAs(outname + ".pdf");
    cout << "[INFO] Saved: " << outname << ".pdf" << endl;
    delete c;
}

// Run DY estimation for one period, write results to txt file
void RunPeriod(const PeriodConfig& cfg, FILE* ftxt) {

    TString base = cfg.base;

    auto makePaths = [&](const vector<TString>& names, const TString& ch) {
        vector<TString> paths;
        for (const auto& n : names)
            paths.push_back(Form("%s/%s/%s/%s.root", base.Data(), ch.Data(), n.Data(), n.Data()));
        return paths;
    };

    vector<TString> DY50_only = { "DYJetsToLL_M_50" };
    vector<TString> DY10_only = { "DYJetsToLL_M_10To50" };

    auto dy_mumu   = makePaths(cfg.DY_names,    "MuMu");
    auto dy_elel   = makePaths(cfg.DY_names,    "ElEl");
    auto dy_muel   = makePaths(cfg.DY_names,    "MuEl");
    auto dat_mumu  = makePaths(cfg.Data_MuMu,   "MuMu");
    auto dat_elel  = makePaths(cfg.Data_ElEl,   "ElEl");
    auto dat_muel  = makePaths(cfg.Data_MuEl,   "MuEl");
    auto dy50_mumu = makePaths(DY50_only,       "MuMu");
    auto dy10_mumu = makePaths(DY10_only,       "MuMu");
    auto dy50_elel = makePaths(DY50_only,       "ElEl");
    auto dy10_elel = makePaths(DY10_only,       "ElEl");
    auto dy50_muel = makePaths(DY50_only,       "MuEl");
    auto dy10_muel = makePaths(DY10_only,       "MuEl");

    // Steps: 2=Dilepton, 3=>=2jets, 5=>=1btag
    vector<int>     steps      = { 2,          3,           5 };
    vector<TString> stepLabels = { "Dilepton", ">= 2 jets", ">= 1 b-tag" };

    // Load step 5 NoMETCut histograms for plotting (outside step loop)
    TH1D* plot_dat_mumu  = SumHistograms(dat_mumu,  "h_cf_NoMETCut_DiLepMass_5");
    TH1D* plot_dat_elel  = SumHistograms(dat_elel,  "h_cf_NoMETCut_DiLepMass_5");
    TH1D* plot_dat_muel  = SumHistograms(dat_muel,  "h_cf_NoMETCut_DiLepMass_5");
    TH1D* plot_dy50_mumu = SumHistograms(dy50_mumu, "h_cf_NoMETCut_DiLepMass_5");
    TH1D* plot_dy10_mumu = SumHistograms(dy10_mumu, "h_cf_NoMETCut_DiLepMass_5");
    TH1D* plot_dy50_elel = SumHistograms(dy50_elel, "h_cf_NoMETCut_DiLepMass_5");
    TH1D* plot_dy10_elel = SumHistograms(dy10_elel, "h_cf_NoMETCut_DiLepMass_5");
    TH1D* plot_dy50_muel = SumHistograms(dy50_muel, "h_cf_NoMETCut_DiLepMass_5");
    TH1D* plot_dy10_muel = SumHistograms(dy10_muel, "h_cf_NoMETCut_DiLepMass_5");

    // Store results for all steps
    vector<StepSF> sf_table;

    // Store step 5 results for final summary
    double SF_ee_final = 1.0, SF_ee_final_e = 0.0;
    double SF_mm_final = 1.0, SF_mm_final_e = 0.0;
    double SF_em_final = 1.0, SF_em_final_e = 0.0;
    double N_outMC_ee_final = 0.0, N_outMC_ee_final_e = 0.0;
    double N_outMC_mm_final = 0.0, N_outMC_mm_final_e = 0.0;
    double N_MC_em_final    = 0.0, N_MC_em_final_e    = 0.0;

    // Table 7 style header
    cout << "\n\n";
    cout << "========================================" << endl;
    cout << " " << cfg.period                          << endl;
    cout << "========================================" << endl;
    cout << "\n  DY SF stability across selection steps (" << cfg.period << ")" << endl;
    cout << "  " << string(66, '-') << endl;
    cout << Form("  %-14s  %18s  %18s  %18s",
                 "Selection", "e+e-", "mu+mu-", "e+-mu-+") << endl;
    cout << "  " << string(66, '-') << endl;

    for (int istep = 0; istep < (int)steps.size(); ++istep) {

        int step = steps[istep];
        TString h_nomet = Form("h_cf_NoMETCut_DiLepMass_%d", step);
        TString h_met   = Form("h_cf_DiLepMass_%d",          step);

        TH1D* h_nomet_dy_mumu  = SumHistograms(dy_mumu,  h_nomet);
        TH1D* h_nomet_dy_elel  = SumHistograms(dy_elel,  h_nomet);
        TH1D* h_nomet_dy_muel  = SumHistograms(dy_muel,  h_nomet);
        TH1D* h_nomet_dat_mumu = SumHistograms(dat_mumu, h_nomet);
        TH1D* h_nomet_dat_elel = SumHistograms(dat_elel, h_nomet);
        TH1D* h_nomet_dat_muel = SumHistograms(dat_muel, h_nomet);

        TH1D* h_met_dy_mumu    = SumHistograms(dy_mumu,  h_met);
        TH1D* h_met_dy_elel    = SumHistograms(dy_elel,  h_met);
        TH1D* h_met_dy_muel    = SumHistograms(dy_muel,  h_met);
        TH1D* h_met_dat_mumu   = SumHistograms(dat_mumu, h_met);
        TH1D* h_met_dat_elel   = SumHistograms(dat_elel, h_met);
        TH1D* h_met_dat_muel   = SumHistograms(dat_muel, h_met);

        // k factors from NoMETCut MC
        double N_inMC_ee_nom   = GetZmassIn(h_nomet_dy_elel).Evt;
        double N_inMC_ee_nom_e = GetZmassIn(h_nomet_dy_elel).Err;
        double N_inMC_mm_nom   = GetZmassIn(h_nomet_dy_mumu).Evt;
        double N_inMC_mm_nom_e = GetZmassIn(h_nomet_dy_mumu).Err;

        double k_ee   = sqrt(N_inMC_ee_nom / N_inMC_mm_nom);
        double k_ee_e = UncKll(N_inMC_ee_nom, N_inMC_ee_nom_e, N_inMC_mm_nom, N_inMC_mm_nom_e);
        double k_mm   = sqrt(N_inMC_mm_nom / N_inMC_ee_nom);
        double k_mm_e = UncKll(N_inMC_mm_nom, N_inMC_mm_nom_e, N_inMC_ee_nom, N_inMC_ee_nom_e);

        double SF_ee_nomet = 1.0, SF_ee_nomet_e = 0.0;
        double SF_mm_nomet = 1.0, SF_mm_nomet_e = 0.0;
        double SF_em_nomet = 1.0, SF_em_nomet_e = 0.0;
        double N_MC_em_nomet = 0.0, N_MC_em_nomet_e = 0.0;
        double SF_ee_met = 1.0, SF_ee_met_e = 0.0;
        double SF_mm_met = 1.0, SF_mm_met_e = 0.0;
        double N_outMC_ee_met = 0.0, N_outMC_ee_met_e = 0.0;
        double N_outMC_mm_met = 0.0, N_outMC_mm_met_e = 0.0;

        for (int imet = 0; imet < 2; ++imet) {

            TH1D* dy_ee_h  = (imet == 0) ? h_nomet_dy_elel  : h_met_dy_elel;
            TH1D* dy_mm_h  = (imet == 0) ? h_nomet_dy_mumu  : h_met_dy_mumu;
            TH1D* dy_em_h  = (imet == 0) ? h_nomet_dy_muel  : h_met_dy_muel;
            TH1D* dat_ee_h = (imet == 0) ? h_nomet_dat_elel : h_met_dat_elel;
            TH1D* dat_mm_h = (imet == 0) ? h_nomet_dat_mumu : h_met_dat_mumu;
            TH1D* dat_em_h = (imet == 0) ? h_nomet_dat_muel : h_met_dat_muel;

            double N_inMC_ee    = GetZmassIn (dy_ee_h).Evt;
            double N_inMC_ee_e  = GetZmassIn (dy_ee_h).Err;
            double N_outMC_ee   = GetZmassOut(dy_ee_h).Evt;
            double N_outMC_ee_e = GetZmassOut(dy_ee_h).Err;
            double R_ee         = N_outMC_ee / N_inMC_ee;
            double R_ee_e       = UncR(N_outMC_ee, N_outMC_ee_e, N_inMC_ee, N_inMC_ee_e);

            double N_inMC_mm    = GetZmassIn (dy_mm_h).Evt;
            double N_inMC_mm_e  = GetZmassIn (dy_mm_h).Err;
            double N_outMC_mm   = GetZmassOut(dy_mm_h).Evt;
            double N_outMC_mm_e = GetZmassOut(dy_mm_h).Err;
            double R_mm         = N_outMC_mm / N_inMC_mm;
            double R_mm_e       = UncR(N_outMC_mm, N_outMC_mm_e, N_inMC_mm, N_inMC_mm_e);

            double N_inDat_ee   = GetZmassIn(dat_ee_h).Evt;
            double N_inDat_ee_e = GetZmassIn(dat_ee_h).Err;
            double N_inDat_mm   = GetZmassIn(dat_mm_h).Evt;
            double N_inDat_mm_e = GetZmassIn(dat_mm_h).Err;
            double N_inDat_em   = GetZmassIn(dat_em_h).Evt;
            double N_inDat_em_e = GetZmassIn(dat_em_h).Err;

            double N_out_ee   = R_ee * (N_inDat_ee - 0.5 * N_inDat_em * k_ee);
            double N_out_ee_e = UncNoutD(R_ee, R_ee_e, N_inDat_ee, N_inDat_ee_e,
                                         N_inDat_em, N_inDat_em_e, k_ee, k_ee_e);
            double N_out_mm   = R_mm * (N_inDat_mm - 0.5 * N_inDat_em * k_mm);
            double N_out_mm_e = UncNoutD(R_mm, R_mm_e, N_inDat_mm, N_inDat_mm_e,
                                         N_inDat_em, N_inDat_em_e, k_mm, k_mm_e);

            double SF_ee   = N_out_ee / N_outMC_ee;
            double SF_ee_e = UncR(N_out_ee, N_out_ee_e, N_outMC_ee, N_outMC_ee_e);
            double SF_mm   = N_out_mm / N_outMC_mm;
            double SF_mm_e = UncR(N_out_mm, N_out_mm_e, N_outMC_mm, N_outMC_mm_e);
            double SF_em   = sqrt(SF_ee * SF_mm);
            double SF_em_e = UncSFem(SF_ee, SF_ee_e, SF_mm, SF_mm_e);

            double N_MC_em   = dy_em_h->Integral();
            double N_MC_em_e = 0.0;
            for (int i = 1; i <= dy_em_h->GetNbinsX(); ++i)
                N_MC_em_e += pow(dy_em_h->GetBinError(i), 2);
            N_MC_em_e = sqrt(N_MC_em_e);

            double N_DY_em   = N_MC_em * SF_em;
            double N_DY_em_e = UncFinal(N_MC_em, N_MC_em_e, SF_em, SF_em_e);

            // Print detailed table for step 5 only
            if (step == 5) {
                TString tag = (imet == 0) ? "NoMETCut" : "METCut";
                PrintTable(tag,
                           N_inMC_ee,  N_inMC_ee_e,  N_outMC_ee, N_outMC_ee_e,
                           R_ee,       R_ee_e,        k_ee,       k_ee_e,
                           N_inDat_ee, N_inDat_ee_e,
                           N_out_ee,   N_out_ee_e,    SF_ee,      SF_ee_e,
                           N_inMC_mm,  N_inMC_mm_e,  N_outMC_mm, N_outMC_mm_e,
                           R_mm,       R_mm_e,        k_mm,       k_mm_e,
                           N_inDat_mm, N_inDat_mm_e,
                           N_out_mm,   N_out_mm_e,    SF_mm,      SF_mm_e,
                           N_inDat_em, N_inDat_em_e,
                           SF_em,      SF_em_e,
                           N_MC_em,    N_MC_em_e,
                           N_DY_em,    N_DY_em_e);
            }

            if (imet == 0) {
                SF_ee_nomet    = SF_ee;      SF_ee_nomet_e    = SF_ee_e;
                SF_mm_nomet    = SF_mm;      SF_mm_nomet_e    = SF_mm_e;
                SF_em_nomet    = SF_em;      SF_em_nomet_e    = SF_em_e;
                N_MC_em_nomet  = N_MC_em;    N_MC_em_nomet_e  = N_MC_em_e;
            } else {
                SF_ee_met      = SF_ee;      SF_ee_met_e      = SF_ee_e;
                SF_mm_met      = SF_mm;      SF_mm_met_e      = SF_mm_e;
                N_outMC_ee_met = N_outMC_ee; N_outMC_ee_met_e = N_outMC_ee_e;
                N_outMC_mm_met = N_outMC_mm; N_outMC_mm_met_e = N_outMC_mm_e;
            }
        } // end imet loop

        // Table 7 style row
        cout << Form("  %-14s  %18s  %18s  %18s",
                     stepLabels[istep].Data(),
                     Form("%.3f +- %.3f", SF_ee_met,   SF_ee_met_e),
                     Form("%.3f +- %.3f", SF_mm_met,   SF_mm_met_e),
                     Form("%.3f +- %.3f", SF_em_nomet, SF_em_nomet_e)) << endl;

        // Compute final yields for this step
        double NDY_ee_step   = N_outMC_ee_met * SF_ee_met;
        double NDY_ee_step_e = UncFinal(N_outMC_ee_met, N_outMC_ee_met_e, SF_ee_met, SF_ee_met_e);
        double NDY_mm_step   = N_outMC_mm_met * SF_mm_met;
        double NDY_mm_step_e = UncFinal(N_outMC_mm_met, N_outMC_mm_met_e, SF_mm_met, SF_mm_met_e);
        double NDY_em_step   = N_MC_em_nomet  * SF_em_nomet;
        double NDY_em_step_e = UncFinal(N_MC_em_nomet, N_MC_em_nomet_e, SF_em_nomet, SF_em_nomet_e);

        // Push this step into sf_table
        StepSF row;
        row.label   = stepLabels[istep];
        row.SF_ee   = SF_ee_met;      row.SF_ee_e  = SF_ee_met_e;
        row.SF_mm   = SF_mm_met;      row.SF_mm_e  = SF_mm_met_e;
        row.SF_em   = SF_em_nomet;    row.SF_em_e  = SF_em_nomet_e;
        row.NMC_ee  = N_outMC_ee_met; row.NMC_ee_e = N_outMC_ee_met_e;
        row.NMC_mm  = N_outMC_mm_met; row.NMC_mm_e = N_outMC_mm_met_e;
        row.NMC_em  = N_MC_em_nomet;  row.NMC_em_e = N_MC_em_nomet_e;
        row.NDY_ee  = NDY_ee_step;    row.NDY_ee_e = NDY_ee_step_e;
        row.NDY_mm  = NDY_mm_step;    row.NDY_mm_e = NDY_mm_step_e;
        row.NDY_em  = NDY_em_step;    row.NDY_em_e = NDY_em_step_e;
        sf_table.push_back(row);

        // Save step 5 results for final summary
        if (step == 5) {
            SF_ee_final      = SF_ee_met;      SF_ee_final_e      = SF_ee_met_e;
            SF_mm_final      = SF_mm_met;      SF_mm_final_e      = SF_mm_met_e;
            SF_em_final      = SF_em_nomet;    SF_em_final_e      = SF_em_nomet_e;
            N_outMC_ee_final = N_outMC_ee_met; N_outMC_ee_final_e = N_outMC_ee_met_e;
            N_outMC_mm_final = N_outMC_mm_met; N_outMC_mm_final_e = N_outMC_mm_met_e;
            N_MC_em_final    = N_MC_em_nomet;  N_MC_em_final_e    = N_MC_em_nomet_e;
        }
    } // end istep loop
    cout << "  " << string(66, '-') << endl;

    // Final DY estimates
    double N_DY_ee   = N_outMC_ee_final * SF_ee_final;
    double N_DY_ee_e = UncFinal(N_outMC_ee_final, N_outMC_ee_final_e, SF_ee_final, SF_ee_final_e);
    double N_DY_mm   = N_outMC_mm_final * SF_mm_final;
    double N_DY_mm_e = UncFinal(N_outMC_mm_final, N_outMC_mm_final_e, SF_mm_final, SF_mm_final_e);
    double N_DY_muel   = N_MC_em_final * SF_em_final;
    double N_DY_muel_e = UncFinal(N_MC_em_final, N_MC_em_final_e, SF_em_final, SF_em_final_e);

    cout << "\n  Final DY Estimates (>= 1 b-tag)" << endl;
    cout << "  " << string(66, '-') << endl;
    cout << Form("  %-22s  %18s  %18s  %18s",
                 "Channel", "e+e-", "mu+mu-", "e+-mu-+") << endl;
    cout << "  " << string(66, '-') << endl;
    cout << Form("  %-22s  %18s  %18s  %18s",
                 "SF",
                 Form("%.4f +- %.4f", SF_ee_final, SF_ee_final_e),
                 Form("%.4f +- %.4f", SF_mm_final, SF_mm_final_e),
                 Form("%.4f +- %.4f", SF_em_final, SF_em_final_e)) << endl;
    cout << Form("  %-22s  %18s  %18s  %18s",
                 "N_MC",
                 Form("%.1f +- %.1f", N_outMC_ee_final, N_outMC_ee_final_e),
                 Form("%.1f +- %.1f", N_outMC_mm_final, N_outMC_mm_final_e),
                 Form("%.1f +- %.1f", N_MC_em_final,    N_MC_em_final_e)) << endl;
    cout << Form("  %-22s  %18s  %18s  %18s",
                 "N_DY (Data-driven)",
                 Form("%.1f +- %.1f", N_DY_ee,   N_DY_ee_e),
                 Form("%.1f +- %.1f", N_DY_mm,   N_DY_mm_e),
                 Form("%.1f +- %.1f", N_DY_muel, N_DY_muel_e)) << endl;
    cout << "  " << string(66, '-') << endl;
    cout << "  Note: MuEl SF = NoMETCut SF, MuMu/ElEl SF = METCut SF" << endl;

    // Write to txt file (all steps)
    fprintf(ftxt, "# %s\n", cfg.period.Data());
    fprintf(ftxt, "# key  SF  SF_err  N_MC  N_MC_err  N_DY  N_DY_err\n");
    for (int is = 0; is < (int)sf_table.size(); ++is) {
        fprintf(ftxt, "%s_%s_ElEl  %.6f  %.6f  %.4f  %.4f  %.4f  %.4f\n",
                cfg.period.Data(), sf_table[is].label.Data(),
                sf_table[is].SF_ee,   sf_table[is].SF_ee_e,
                sf_table[is].NMC_ee,  sf_table[is].NMC_ee_e,
                sf_table[is].NDY_ee,  sf_table[is].NDY_ee_e);
        fprintf(ftxt, "%s_%s_MuMu  %.6f  %.6f  %.4f  %.4f  %.4f  %.4f\n",
                cfg.period.Data(), sf_table[is].label.Data(),
                sf_table[is].SF_mm,   sf_table[is].SF_mm_e,
                sf_table[is].NMC_mm,  sf_table[is].NMC_mm_e,
                sf_table[is].NDY_mm,  sf_table[is].NDY_mm_e);
        fprintf(ftxt, "%s_%s_MuEl  %.6f  %.6f  %.4f  %.4f  %.4f  %.4f\n",
                cfg.period.Data(), sf_table[is].label.Data(),
                sf_table[is].SF_em,   sf_table[is].SF_em_e,
                sf_table[is].NMC_em,  sf_table[is].NMC_em_e,
                sf_table[is].NDY_em,  sf_table[is].NDY_em_e);
    }
    fprintf(ftxt, "#\n");

    // Draw plots
    cout << "\n[INFO] Drawing plots for " << cfg.period << "..." << endl;
    DrawPlot(plot_dat_mumu, plot_dy50_mumu, plot_dy10_mumu,
             "M_{#mu^{+}#mu^{-}} (GeV)", "NanoAODv15_" + cfg.period + "_mumu_nomet");
    DrawPlot(plot_dat_elel, plot_dy50_elel, plot_dy10_elel,
             "M_{e^{+}e^{-}} (GeV)", "NanoAODv15_" + cfg.period + "_elel_nomet");
    DrawPlot(plot_dat_muel, plot_dy50_muel, plot_dy10_muel,
             "M_{e#mu} (GeV)", "NanoAODv15_" + cfg.period + "_muel_nomet");
}

void CalDYSF_Run2UL_v15() {

    TString basedir = "/Users/seungkyuha/SK_Analysis_v2/CPV/Run2Full_v1/AnalysisCodes/BGEstimation_v1/NanoAODv15_DYEst_v1";

    // -------------------------------------------------------
    // UL2017 configuration
    // -------------------------------------------------------
    PeriodConfig cfg2017;
    cfg2017.period = "UL2017";
    cfg2017.base   = basedir + "/UL2017";
    cfg2017.DY_names = { "DYJetsToLL_M_10To50", "DYJetsToLL_M_50" };
    cfg2017.Data_MuMu = {
        "Data_DoubleMuon_Run2017B", "Data_DoubleMuon_Run2017C",
        "Data_DoubleMuon_Run2017D", "Data_DoubleMuon_Run2017E",
        "Data_DoubleMuon_Run2017F",
        "Data_SingleMuon_Run2017B", "Data_SingleMuon_Run2017C",
        "Data_SingleMuon_Run2017D", "Data_SingleMuon_Run2017E",
        "Data_SingleMuon_Run2017F"
    };
    cfg2017.Data_ElEl = {
        "Data_DoubleEG_Run2017B",        "Data_DoubleEG_Run2017C",
        "Data_DoubleEG_Run2017D",        "Data_DoubleEG_Run2017E",
        "Data_DoubleEG_Run2017F",
        "Data_SingleElectron_Run2017B",  "Data_SingleElectron_Run2017C",
        "Data_SingleElectron_Run2017D",  "Data_SingleElectron_Run2017E",
        "Data_SingleElectron_Run2017F"
    };
    cfg2017.Data_MuEl = {
        "Data_MuonEG_Run2017B",       "Data_MuonEG_Run2017C",
        "Data_MuonEG_Run2017D",       "Data_MuonEG_Run2017E",
        "Data_MuonEG_Run2017F",
        "Data_SingleElectron_Run2017B", "Data_SingleElectron_Run2017C",
        "Data_SingleElectron_Run2017D", "Data_SingleElectron_Run2017E",
        "Data_SingleElectron_Run2017F",
        "Data_SingleMuon_Run2017B",   "Data_SingleMuon_Run2017C",
        "Data_SingleMuon_Run2017D",   "Data_SingleMuon_Run2017E",
        "Data_SingleMuon_Run2017F"
    };

    // -------------------------------------------------------
    // UL2018 configuration
    // -------------------------------------------------------
    PeriodConfig cfg2018;
    cfg2018.period = "UL2018";
    cfg2018.base   = basedir + "/UL2018";
    cfg2018.DY_names = { "DYJetsToLL_M_10To50", "DYJetsToLL_M_50" };
    cfg2018.Data_MuMu = {
        "Data_DoubleMuon_Run2018A", "Data_DoubleMuon_Run2018B",
        "Data_DoubleMuon_Run2018C", "Data_DoubleMuon_Run2018D",
        "Data_SingleMuon_Run2018A", "Data_SingleMuon_Run2018B",
        "Data_SingleMuon_Run2018C", "Data_SingleMuon_Run2018D"
    };
    cfg2018.Data_ElEl = {
        "Data_EGamma_Run2018A", "Data_EGamma_Run2018B",
        "Data_EGamma_Run2018C", "Data_EGamma_Run2018D"
    };
    cfg2018.Data_MuEl = {
        "Data_MuonEG_Run2018A",     "Data_MuonEG_Run2018B",
        "Data_MuonEG_Run2018C",     "Data_MuonEG_Run2018D",
        "Data_SingleMuon_Run2018A", "Data_SingleMuon_Run2018B",
        "Data_SingleMuon_Run2018C", "Data_SingleMuon_Run2018D",
        "Data_EGamma_Run2018A",     "Data_EGamma_Run2018B",
        "Data_EGamma_Run2018C",     "Data_EGamma_Run2018D"
    };

    // Open output txt file
    TString outTxt = "DYSF_Run2UL_NanoAODv15.txt";
    FILE* ftxt = fopen(outTxt.Data(), "w");
    if (!ftxt) {
        cerr << "[ERROR] Cannot open " << outTxt << " for writing" << endl;
        return;
    }
    fprintf(ftxt, "# DY Scale Factors Run2 UL (>= 1 b-tag)\n");
    fprintf(ftxt, "# Note: MuMu/ElEl SF from METCut, MuEl SF from NoMETCut\n");
    fprintf(ftxt, "# Columns: key  SF  SF_err  N_MC  N_MC_err  N_DY  N_DY_err\n");
    fprintf(ftxt, "#\n");

    // -------------------------------------------------------
    // UL2016PreVFP configuration
    // -------------------------------------------------------
    PeriodConfig cfg2016Pre;
    cfg2016Pre.period = "UL2016PreVFP";
    cfg2016Pre.base   = basedir + "/UL2016PreVFP";
    cfg2016Pre.DY_names = { "DYJetsToLL_M_10To50", "DYJetsToLL_M_50" };
    // NanoAODv15: PreVFP eras are named *-HIPM; Run2016B-HIPM (ver1) has no
    // golden-JSON lumis (all jobs EMPTY_INPUT, no merged file), so only
    // Run2016B-HIPM-v2 is used for era B.
    cfg2016Pre.Data_MuMu = {
        "Data_DoubleMuon_Run2016B-HIPM-v2", "Data_DoubleMuon_Run2016C-HIPM",
        "Data_DoubleMuon_Run2016D-HIPM",    "Data_DoubleMuon_Run2016E-HIPM",
        "Data_DoubleMuon_Run2016F-HIPM",
        "Data_SingleMuon_Run2016B-HIPM-v2", "Data_SingleMuon_Run2016C-HIPM",
        "Data_SingleMuon_Run2016D-HIPM",    "Data_SingleMuon_Run2016E-HIPM",
        "Data_SingleMuon_Run2016F-HIPM"
    };
    cfg2016Pre.Data_ElEl = {
        "Data_DoubleEG_Run2016B-HIPM-v2",       "Data_DoubleEG_Run2016C-HIPM",
        "Data_DoubleEG_Run2016D-HIPM",          "Data_DoubleEG_Run2016E-HIPM",
        "Data_DoubleEG_Run2016F-HIPM",
        "Data_SingleElectron_Run2016B-HIPM-v2", "Data_SingleElectron_Run2016C-HIPM",
        "Data_SingleElectron_Run2016D-HIPM",    "Data_SingleElectron_Run2016E-HIPM",
        "Data_SingleElectron_Run2016F-HIPM"
    };
    cfg2016Pre.Data_MuEl = {
        "Data_MuonEG_Run2016B-HIPM-v2",         "Data_MuonEG_Run2016C-HIPM",
        "Data_MuonEG_Run2016D-HIPM",            "Data_MuonEG_Run2016E-HIPM",
        "Data_MuonEG_Run2016F-HIPM",
        "Data_SingleElectron_Run2016B-HIPM-v2", "Data_SingleElectron_Run2016C-HIPM",
        "Data_SingleElectron_Run2016D-HIPM",    "Data_SingleElectron_Run2016E-HIPM",
        "Data_SingleElectron_Run2016F-HIPM",
        "Data_SingleMuon_Run2016B-HIPM-v2",     "Data_SingleMuon_Run2016C-HIPM",
        "Data_SingleMuon_Run2016D-HIPM",        "Data_SingleMuon_Run2016E-HIPM",
        "Data_SingleMuon_Run2016F-HIPM"
    };

    // -------------------------------------------------------
    // UL2016PostVFP configuration
    // -------------------------------------------------------
    PeriodConfig cfg2016Post;
    cfg2016Post.period = "UL2016PostVFP";
    cfg2016Post.base   = basedir + "/UL2016PostVFP";
    cfg2016Post.DY_names = { "DYJetsToLL_M_10To50", "DYJetsToLL_M_50" };
    cfg2016Post.Data_MuMu = {
        "Data_DoubleMuon_Run2016F",      "Data_DoubleMuon_Run2016G",
        "Data_DoubleMuon_Run2016H",
        "Data_SingleMuon_Run2016F",      "Data_SingleMuon_Run2016G",
        "Data_SingleMuon_Run2016H"
    };
    cfg2016Post.Data_ElEl = {
        "Data_DoubleEG_Run2016F",          "Data_DoubleEG_Run2016G",
        "Data_DoubleEG_Run2016H",
        "Data_SingleElectron_Run2016F",    "Data_SingleElectron_Run2016G",
        "Data_SingleElectron_Run2016H"
    };
    cfg2016Post.Data_MuEl = {
        "Data_MuonEG_Run2016F",            "Data_MuonEG_Run2016G",
        "Data_MuonEG_Run2016H",
        "Data_SingleElectron_Run2016F",    "Data_SingleElectron_Run2016G",
        "Data_SingleElectron_Run2016H",
        "Data_SingleMuon_Run2016F",        "Data_SingleMuon_Run2016G",
        "Data_SingleMuon_Run2016H"
    };

    // Run each period
    RunPeriod(cfg2016Pre,  ftxt);
    RunPeriod(cfg2016Post, ftxt);
    RunPeriod(cfg2017, ftxt);
    RunPeriod(cfg2018, ftxt);

    fclose(ftxt);
    cout << "\n[INFO] Saved: " << outTxt << endl;
}
