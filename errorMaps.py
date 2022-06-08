from collections import OrderedDict

error_bands = OrderedDict()

error_bands["Detector"]                = ["AngleXZ",
                                          "AngleYZ",
                                          "LYAtt",
                                          "LYRay",
                                          "LY",
                                          "Recom2",
                                          "SCE",
                                          "WireX",
                                          "WireYZ"]

error_bands["Flux"]                    = ["expskin_FluxUnisim",
                                          "horncurrent_FluxUnisim",
                                          "kminus_PrimaryHadronNormalization",
                                          "kplus_PrimaryHadronFeynmanScaling",
                                          "kzero_PrimaryHadronSanfordWang",
                                          "nucleoninexsec_FluxUnisim",
                                          "nucleonqexsec_FluxUnisim",
                                          "nucleontotxsec_FluxUnisim",
                                          "piminus_PrimaryHadronSWCentralSplineVariation",
                                          "pioninexsec_FluxUnisim",
                                          "pionqexsec_FluxUnisim",
                                          "piontotxsec_FluxUnisim",
                                          "piplus_PrimaryHadronSWCentralSplineVariation"]

error_bands["GEANT4"]                  = ["piminus",
                                          "piplus",
                                          "proton"]

error_bands["GENIE"]                   = ["All_UBGenie",
                                          "AxFFCCQEshape_UBGenie",
                                          "DecayAngMEC_UBGenie",
                                          "NormCCCOH_UBGenie",
                                          "NormNCCOH_UBGenie",
                                          "RPA_CCQE_UBGenie",
                                          "Theta_Delta2Npi_UBGenie",
                                          "VecFFCCQEshape_UBGenie"]

error_bands["POT"]                     = ["POT_variation"]

error_bands["target"]                  = ["target_variation"]

