import os
user = os.environ['USER']

#TOP_DIRECTORY = "/tmp/" + user + "/HmumuTuplesVer_2_4/"
TOP_DIRECTORY = "/Users/lukasadamek/cernbox/EOPTuples/haddedTogether/EOPTuples_v04/"

INPUT = \
{ "PythiaJetJet": ("Pythia8 LowMu and Dijet", [\
TOP_DIRECTORY + "user.luadamek.mc16_13TeV.361020.jetjet.Lukas_NewSpookyTreeLabelsOct30_hist",\
TOP_DIRECTORY + "user.luadamek.mc16_13TeV.361021.jetjet.Lukas_NewSpookyTreeLabelsOct30_hist",\
TOP_DIRECTORY + "user.luadamek.mc16_13TeV.361022.jetjet.Lukas_NewSpookyTreeLabelsOct30_hist",\
]),\
"LowMuData": ("2017 Data, Low-<#mu> Run 341294",[\
TOP_DIRECTORY + "user.luadamek.luadamek.data17_13TeV.00341294.lukas_NewSpookyTreeLabelsOct30_hist",\
])}
