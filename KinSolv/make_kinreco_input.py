#!/usr/bin/env python3
"""Build a KinReco input file from ssb_analysis outputs made with MakeKinRecoInput : "True".

Sums the KinReco_* histograms over all given files (e.g. the TTbar_Signal outputs of
the dielec / dimuon / muelec jobs of one era), clips negative bins to 0 (TH1::GetRandom
needs non-negative contents) and writes them with the names KinematicReconstruction::
loadData() expects. Optionally draws a shape comparison against a reference file.

  python3 KinSolv/make_kinreco_input.py --out KinSolv/KinReco_input_UL2018.root \
      --ref KinSolv/Kin_input.root --plots plots_kinreco_UL2018 \
      output/<study>/UL2018/*/UL2018_TTbar_Signal*.root
"""
import argparse
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

SOLVER_HISTS = [
    "KinReco_W_mass_step0",
    "KinReco_mbl_true_step0",
    "KinReco_d_angle_jet_step7",
    "KinReco_d_angle_lep_step7",
    "KinReco_fE_jet_step7",
    "KinReco_fE_lep_step7",
]
DIAG_HISTS = [
    "KinRecoDiag_fE_jet_vs_pt_step7",
    "KinRecoDiag_fE_lep_vs_pt_step7",
    "KinRecoDiag_dR_lep_step7",
    "KinRecoDiag_dR_jet_step7",
    "KinRecoDiag_count",
]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="ssb_analysis output files (TTbar_Signal, MakeKinRecoInput=True)")
    ap.add_argument("--out", required=True, help="output KinReco input file")
    ap.add_argument("--ref", default=None, help="reference input file to compare with (e.g. KinSolv/Kin_input.root)")
    ap.add_argument("--plots", default=None, help="directory for comparison plots (needs --ref)")
    ap.add_argument("--keep-diag", action="store_true", help="also store the KinRecoDiag_* histograms")
    args = ap.parse_args()

    names = SOLVER_HISTS + (DIAG_HISTS if args.keep_diag else [])
    summed = {}
    nfiles = 0
    for path in args.inputs:
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            print(f"[WARN] cannot open {path} - skipped")
            continue
        if not f.Get(SOLVER_HISTS[0]):
            print(f"[WARN] {path} has no KinReco_* histograms (MakeKinRecoInput off?) - skipped")
            f.Close()
            continue
        nfiles += 1
        for n in names:
            h = f.Get(n)
            if not h:
                continue
            if n in summed:
                summed[n].Add(h)
            else:
                summed[n] = h.Clone(n)
        f.Close()
    if nfiles == 0:
        sys.exit("[ERROR] no usable input file")
    missing = [n for n in SOLVER_HISTS if n not in summed]
    if missing:
        sys.exit(f"[ERROR] missing histograms: {missing}")

    print(f"[INFO] summed {nfiles} files")
    for n in SOLVER_HISTS:
        h = summed[n]
        nneg = 0
        for b in range(0, h.GetNbinsX() + 2):
            if h.GetBinContent(b) < 0:
                h.SetBinContent(b, 0.)
                h.SetBinError(b, 0.)
                nneg += 1
        print(f"  {n:28s} integral={h.Integral():14.1f} mean={h.GetMean():8.4f} rms={h.GetRMS():8.4f}"
              + (f"  (clipped {nneg} negative bins)" if nneg else ""))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    fout = ROOT.TFile(args.out, "RECREATE")
    for n in names:
        if n in summed:
            summed[n].Write(n)
    fout.Close()
    print(f"[INFO] wrote {args.out}")

    if args.ref and args.plots:
        compare(summed, args.ref, args.plots)


def compare(new, ref_path, outdir):
    os.makedirs(outdir, exist_ok=True)
    fref = ROOT.TFile.Open(ref_path)
    c = ROOT.TCanvas("c", "c", 700, 600)
    for n in SOLVER_HISTS:
        href = fref.Get(n)
        if not href:
            continue
        a = new[n].Clone(n + "_new")
        b = href.Clone(n + "_ref")
        for h in (a, b):
            if h.Integral() > 0:
                h.Scale(1. / h.Integral())
        a.SetLineColor(ROOT.kRed + 1)
        a.SetLineWidth(2)
        b.SetLineColor(ROOT.kBlack)
        b.SetLineWidth(2)
        b.SetTitle(n)
        b.SetMaximum(1.25 * max(a.GetMaximum(), b.GetMaximum()))
        b.Draw("hist")
        a.Draw("hist same")
        leg = ROOT.TLegend(0.55, 0.75, 0.88, 0.88)
        leg.AddEntry(b, f"ref (mean {href.GetMean():.3f})", "l")
        leg.AddEntry(a, f"new (mean {new[n].GetMean():.3f})", "l")
        leg.Draw()
        c.SaveAs(os.path.join(outdir, n + ".png"))
    fref.Close()
    print(f"[INFO] comparison plots in {outdir}")


if __name__ == "__main__":
    main()
