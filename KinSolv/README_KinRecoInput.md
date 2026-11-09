# KinReco input histograms

`KinematicReconstruction::loadData()` reads six histograms. The file is set with the config key
`KinRecoInputFile` (default `./KinSolv/Kin_input.root`, the old 2016 80X inputs).

| histogram | definition |
|---|---|
| `KinReco_W_mass_step0` | gen W mass, W+ from t and W- from tbar (last copy) |
| `KinReco_mbl_true_step0` | gen m(l+ b) and m(l- bbar). Lepton: final state (status 1), prompt, from hard process, not from tau, comes from the W. b: the b quark that is a direct top daughter (parton). |
| `KinReco_fE_lep_step7` | E_gen / E_reco, reco lepton matched to the same-charge gen lepton (dR < 0.1) |
| `KinReco_d_angle_lep_step7` | 3D angle between the reco and gen lepton |
| `KinReco_fE_jet_step7` | E_b / E_jet, selected jet matched to the top b quark (closest jet, dR < 0.4, the two b's must not share a jet) |
| `KinReco_d_angle_jet_step7` | 3D angle between the jet and the b quark |

step0 means every gen event of the channel, before any reco cut. step7 means after the full dilepton
selection with >= 1 b-tag, using the same `Lep`/`AnLep`/`jets` that go into the solver.
The histograms are filled with sign(Generator_weight) only. Binning is the same as in the old file.

## Producing them

1. In the era's `dielec.config`, `dimuon.config` and `muelec.config`, set `MakeKinRecoInput : "True"`.
2. Run `TTbar_Signal` for all three channels as usual. The `KinReco_*` and `KinRecoDiag_*` histograms
   are written into the normal output file. The flag is ignored for every other sample and for data.
3. Merge the outputs and compare with the old file:

       python3 KinSolv/make_kinreco_input.py --out KinSolv/KinReco_input_UL2018.root \
           --ref KinSolv/Kin_input.root --plots plots_kinreco_UL2018 <TTbar_Signal outputs of the 3 channels>

4. Set `MakeKinRecoInput : "False"` again and enable
   `KinRecoInputFile : "./KinSolv/KinReco_input_UL2018.root"` in all configs of that era (MC and data).

Check `KinRecoDiag_count` (matching efficiencies) and `KinRecoDiag_fE_*_vs_pt_step7` before using a new file.
After switching, the KinReco data/MC scale factors need to be re-derived.
