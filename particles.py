### particles aka ###
### qpdb aka ###
### quality particle data-base ###

tau_decay_modes_simple = {
  # mu leptonic
  "mu- nu(mu)~ nu(tau)": [{"mu-": 1, "nu(mu)~": 1, "nu(tau)": 1}],                                # 17.4%
  # e leptonic
  "e- nu(e)~ nu(tau)": [{"e-": 1, "nu(e)~": 1, "nu(tau)": 1}],                                    # 17.8%
  # 1-prong 0pi0
  "h- nu(tau)": [{"pi-": 1, "nu(tau)": 1},
                 {"K-": 1, "nu(tau)": 1}],                                                        # 11.5%
  # 1-prong 1pi0
  "h- pi0 nu(tau)": [{"pi-": 1, "gamma": 2, "nu(tau)": 1},
                     {"pi-": 1, "e+": 1, "e-": 1, "gamma": 1, "nu(tau)": 1},
                     {"K-": 1, "gamma": 2, "nu(tau)": 1}],                                        # 25.9%
  # 1-prong 2pi0
  "h- 2pi0 nu(tau)": [{"pi-": 1, "gamma": 4, "nu(tau)": 1},
                      {"pi-": 1, "e+": 1, "e-": 1, "gamma": 3, "nu(tau)": 1},
                      {"K-": 1, "gamma": 4, "nu(tau)": 1}],                                       #  9.7%
  # 1-prong 3pi0
  "h- 3pi0 nu(tau)": [{"pi-": 1, "gamma": 6, "nu(tau)": 1},
                      {"pi-": 1, "e+": 1, "e-": 1, "gamma": 5, "nu(tau)": 1},
                      {"K-": 1, "gamma": 6, "nu(tau)": 1}],                                       #  1.1%
  # 3-prong 0pi0
  "pi- h- h+ nu(tau)": [{"pi-": 2, "pi+": 1, "nu(tau)": 1},
                        {"pi-": 1, "K-": 1, "K+": 1, "nu(tau)": 1}],                              #  9.1%
  # 3(1)-prong 1pi0(omega)
  "pi- omega/(pi- p+ pi0) nu(tau)": [{"pi-": 2, "pi+": 1, "gamma": 2, "nu(tau)": 1},
                                     {"pi-": 2, "pi+": 1, "e+": 1, "e-": 1, "gamma": 1, "nu(tau)": 1},
                                     {"pi-": 1, "gamma": 3, "nu(tau)": 1}],                       #  4.5%
  # 3(1)-prong 2pi0(1pi0 omega)
  #"h- omega/(h- h+ pi0) pi0 nu(tau)": [{"pi-": 2, "pi+": 1, "gamma": 4, "nu(tau)": 1},
  #{"pi-": 1, "K-": 1, "K+": 1, "gamma": 4, "nu(tau)": 1},
  #{"pi-": 2, "pi+": 1, "e+": 1, "e-": 1, "gamma": 3, "nu(tau)": 1},
  #]          #  0.5%
  # other                                                                                            3.0%
}

tau_decay_topologies = {
  "mu- nu(mu)~ nu(tau)": {"mu-": 1, "nu(mu)~": 1, "nu(tau)": 1},                                  # 17.4%
  "e- nu(e)~ nu(tau)": {"e-": 1, "nu(e)~": 1, "nu(tau)": 1},                                      # 17.8%
  "pi- nu(tau)": {"pi-": 1, "nu(tau)": 1},                                                        # 10.8%
  "K- nu(tau)": {"K-": 1, "nu(tau)": 1},                                                          #  0.7%
  "pi- pi0(a) nu(tau)": {"pi-": 1, "gamma": 2, "nu(tau)": 1},                                     # 25.2%
  "pi- pi0(b) nu(tau)": {"pi-": 1, "e+": 1, "e-": 1, "gamma": 1, "nu(tau)": 1},                   #  0.3%
  "K- pi0(a) nu(tau)": {"K-": 1, "gamma": 2, "nu(tau)": 1},                                       #  0.4%
  "pi- 2pi0(aa) nu(tau)": {"pi-": 1, "gamma": 4, "nu(tau)": 1},                                   #  9.0%
  "pi- 2pi0(ab) nu(tau)": {"pi-": 1, "e+": 1, "e-": 1, "gamma": 3, "nu(tau)": 1},                 #  0.1%
  "K- 2pi0(aa) nu(tau)": {"K-": 1, "gamma": 4, "nu(tau)": 1},                                     #  0.6%
  "pi- 3pi0(aaa) nu(tau)": {"pi-": 1, "gamma": 6, "nu(tau)": 1},                                  #  1.0%
  "pi- 3pi0(aab) nu(tau)": {"pi-": 1, "e+": 1, "e-": 1, "gamma": 5, "nu(tau)": 1},                #  0.01%
  "K- 3pi0(aaa) nu(tau)": {"K-": 1, "gamma": 6, "nu(tau)": 1},                                    #  0.05%
  #"pi- K~0 nu(tau)": {"pi-": 1, "K~0": 1, "nu(tau)": 1}, #  0.8%
  #"K- K~0 nu(tau)": {"K-": 1, "K~0": 1, "nu(tau)": 1}, #  0.2%
  "pi- pi- pi+ nu(tau)": {"pi-": 2, "pi+": 1, "nu(tau)": 1},                                      #  9.0%
  "pi- [pi- pi+ pi0(a)]/[omega(aa)]  nu(tau)": {"pi-": 2, "pi+": 1, "gamma": 2, "nu(tau)": 1},    #  4.37%
  "pi- [pi- pi+ pi0(b)]/[omega(ab)] nu(tau)": {"pi-": 2, "pi+": 1, "e+": 1, "e-": 1, "gamma": 1, "nu(tau)": 1}, #  0.06%
  "pi- omega(ba) nu(tau)": {"pi-": 1, "gamma": 3, "nu(tau)": 1},                                  #  0.16%
  # other                                                                                            3.1%
}

pi0_decay_modes = [{"gamma": 2},                    # 98.8%
                   {"e+": 1, "e-": 1, "gamma": 1,}  #  1.2%
                   ] # https://pdglive.lbl.gov/Particle.action?node=S009&home=sumtabM

omega_decay_modes = [{"pi+": 1, "pi-": 1, "pi0": 1}, # 89.2%
                     {"pi0": 1, "gamma": 1},        #  8.4%
                     {"pi+": 1, "pi-": 1},          #  1.5%
                     ] # https://pdglive.lbl.gov/Particle.action?node=M001&home=sumtabM

tau_decay_modes = { # 96.2% of modes categorized
  "mu- nu(mu)~ nu(tau)": {"mu-": 1, "nu(mu)~": 1, "nu(tau)": 1},          # 17.4%
  "e- nu(e)~ nu(tau)": {"e-": 1, "nu(e)~": 1, "nu(tau)": 1},              # 17.8%
  # ---
  "pi- nu(tau)": {"pi-": 1, "nu(tau)": 1},                                # 10.8%
  "K- nu(tau)": {"K-": 1, "nu(tau)": 1},                                  #  0.7%
  # ---
  "pi- pi0 nu(tau)": {"pi-": 1, "pi0": 1, "nu(tau)": 1},                  # 25.5%
  "K- pi0 nu(tau)": {"K-": 1, "pi0": 1, "nu(tau)": 1},                    #  0.4%
  # ---
  "pi- 2pi0 nu(tau)": {"pi-": 1, "pi0": 2, "nu(tau)": 1},                 #  9.3%
  "K- 2pi0 nu(tau)": {"K-": 1, "pi0": 2, "nu(tau)": 1},                   #  0.6%
  # ---
  "pi- 3pi0 nu(tau)": {"pi-": 1, "pi0": 3, "nu(tau)": 1},                 #  1.0%
  "K- 3pi0 nu(tau)": {"K-": 1, "pi0": 3, "nu(tau)": 1},                   #  0.04%
  # ...
  "pi- K~0 nu(tau)": {"pi-": 1, "K~0": 1, "nu(tau)": 1},                  #  0.8%
  "K- K~0 nu(tau)": {"K-": 1, "K~0": 1, "nu(tau)": 1},                    #  0.2%
  # ...
  "pi- pi- pi+ nu(tau)": {"pi-": 2, "pi+": 1, "nu(tau)": 1},              #  9.0%
  "pi- pi- pi+ pi0 nu(tau)": {"pi-": 2, "pi+": 1, "pi0": 1, "nu(tau)": 1},#  2.7%
  # ...
  "pi- omega nu(tau)": {"pi-": 1, "omega(782)": 1, "nu(tau)": 1},         #  1.9%
  # other                                                                    1.9%
} # https://pdglive.lbl.gov/Particle.action?node=S035&init=0

pdg_code = {
  "electron": 11,
  "electron neutrino": 15,
  "muon": 13,
  "muon neutrino": 14,
  "tau": 15,
  "tau neutrino": 16,
  "argon" : 1000180400,
  "pion0" : 111,
  "pion+" : 211,
  "pion-" : -211,
}

mass = {
  "electron": 0.510998,
  "muon": 105.658,
  "tau": 1776.86,
}

neut_code_to_name = {
  1 : "CCQE",
  -1 : "CCQE",
  #2 : "2p2h",
  2 : "MEC", #MEC2p2h
  #-2 : "2p2h",
  -2 : "MEC",
  #11 : "CCRes1pi+",
  11 : "ResCCNuNeutronPiPlus",
  #13 : "CCRes1pi+",
  13 : "ResCCNuProtonPiPlus",
  12 : "CCRes1pi0",
  -12 : "CCRes1pi0",
  -11 : "CCRes1pi-",
  -13 : "CCRes1pi-",
  15 : "CCDif1pi+",
  -15 : "CCDif1pi-",
  16 : "CCCoh1pi+",
  17 : "CCRes1gamma",
  -17 : "CCRes1gamma",
  #21 : "CCNpi",
  21 : "CCNpi (Res)",
  #-21 : "CCNpi",
  -21 : "CCNpi (Res)",
  22 : "CCRes1eta0",
  -22 : "CCRes1eta0",
  23 : "CCRes1K0",
  -23 : "CCRes1K+",
  26 : "CCDis",
  -26 : "CCDis",
  31 : "NCRes1pi0",
  -31 : "NCRes1pi0",
  32 : "NCRes1pi0",
  -32 : "NCRes1pi0",
  33 : "NCRes1pi0",
  -33 : "NCRes1pi0",
  34 : "NCRes1pi+",
  -34 : "NCRes1pi+",
  35 : "NCDif1pi0",
  -35 : "NCDif1pi0",
  36 : "NCCoh1pi0",
  -36 : "NCCoh1pi0",
  38 : "NCRes1gamma",
  -38 : "NCRes1gamma",
  39 : "NCRes1gamma",
  -39 : "NCRes1gamma",
  41 : "NCNpi",
  -41 : "NCNpi",
  42 : "NCRes1eta0",
  -42 : "NCRes1eta0",
  43 : "NCRes1eta0",
  -43 : "NCRes1eta0",
  44 : "NCRes1K0",
  -44 : "NCRes1K0",
  45 : "NCRes1K+",
  -45 : "NCRes1K+",
  46 : "NCDIS",
  -46 : "NCDIS",
  51: "NCEL", #1p1h
  -51: "NCEL",
  52: "NCEL",
  -52: "NCEL",
}
scattering_code_to_name = {
  #0 : "QE",
  0 : "QE (other)",
  #1 : "Res",
  1 : "Res (other)",
  #2 : "DIS",
  2 : "DIS (other)",
  3 : "Coh",
  4 : "CohElastic",
  5 : "ElectronScattering",
  6 : "IMDAnnihilation",
  7 : "InverseBetaDecay",
  8 : "GlashowResonance",
  9 : "AMNuGamma",
  10 : "MEC", #MEC2p2h
  11 : "Diffractive",
  12 : "EM",
  13 : "WeakMix",
}
nuance_code_to_name = {
  1 : "CCQE",
  2 : "NCQE",
  3 : "ResCCNuProtonPiPlus",
  4 : "ResCCNuNeutronPiPlus",
  5 : "ResCCNuNeutronPiPlus",
  6 : "ResNCNuProtonPi0",
  7 : "ResNCNuProtonPiPlus",
  8 : "ResNCNuNeutronPi0",
  9 : "ResNCNuNeutronPiMinus",
  10 : "ResCCNuBarNeutronPiMinus",
  11 : "ResCCNuBarProtonPi0",
  12 : "ResCCNuBarProtonPiMinus",
  13 : "ResNCNuBarProtonPi0",
  14 : "ResNCNuBarProtonPiPlus",
  15 : "ResNCNuBarNeutronPi0",
  16 : "ResNCNuBarNeutronPiMinus",
  17 : "ResCCNuDeltaPlusPiPlus",
  21 : "ResCCNuDelta2PlusPiMinus",
  28 : "ResCCNuBarDelta0PiMinus",
  32 : "ResCCNuBarDeltaMinusPiPlus",
  39 : "ResCCNuProtonRhoPlus",
  41 : "ResCCNuNeutronRhoPlus",
  46 : "ResCCNuBarNeutronRhoMinus",
  48 : "ResCCNuBarNeutronRho0",
  53 : "ResCCNuSigmaPlusKaonPlus",
  55 : "ResCCNuSigmaPlusKaon0",
  60 : "ResCCNuBarSigmaMinusKaon0",
  62 : "ResCCNuBarSigma0Kaon0",
  67 : "ResCCNuProtonEta",
  70 : "ResCCNuBarNeutronEta",
  73 : "ResCCNuKaonPlusLambda0",
  76 : "ResCCNuBarKaon0Lambda0",
  79 : "ResCCNuProtonPiPlusPiMinus",
  80 : "ResCCNuProtonPi0Pi0",
  85 : "ResCCNuBarNeutronPiPlusPiMinus",
  86 : "ResCCNuBarNeutronPi0Pi0",
  90 : "ResCCNuBarProtonPi0Pi0",
  91 : "CCDIS",
  92 : "NCDIS",
  93 : "UnUsed1",
  94 : "UnUsed2",
  95 : "CCQEHyperon",
  96 : "NCCOH",
  97 : "CCCOH",
  98 : "NuElectronElastic",
  99 : "InverseMuDecay",
  100: "MEC", #MEC2p2h
}

