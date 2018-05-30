#' PAM 50 Gene Names Lookup Table
#' 
#' PAM50 gene names from `pam50_annotation.txt` at
#' <https://genome.unc.edu/pubsup/breastGEO/PAM50.zip>
data_pam50_annotation <- function() {
  pam50_annotation <- tibble::tribble(
    ~pcrID,    ~UniGene, ~EntrezGene, ~GeneName, ~InProliferation,
    "ACTR3B", "Hs.647117",      57180L,  "ACTR3B",              "N",
    "ANLN",  "Hs.62180",      54443L,    "ANLN",              "N",
    "BAG1", "Hs.377484",        573L,    "BAG1",              "N",
    "BCL2", "Hs.150749",        596L,    "BCL2",              "N",
    "BIRC5", "Hs.514527",        332L,   "BIRC5",              "Y",
    "BLVRA", "Hs.488143",        644L,   "BLVRA",              "N",
    "CCNB1",  "Hs.23960",        891L,   "CCNB1",              "Y",
    "CCNE1", "Hs.244723",        898L,   "CCNE1",              "N",
    "CDC20", "Hs.524947",        991L,   "CDC20",              "Y",
    "CDC6", "Hs.405958",        990L,    "CDC6",              "N",
    "CDCA1", "Hs.651950",      83540L,    "NUF2",              "Y", #<<
    "CDH3", "Hs.461074",       1001L,    "CDH3",              "N",
    "CENPF", "Hs.497741",       1063L,   "CENPF",              "N",
    "CEP55",  "Hs.14559",      55165L,   "CEP55",              "Y",
    "CXXC5", "Hs.189119",      51523L,   "CXXC5",              "N",
    "EGFR", "Hs.488293",       1956L,    "EGFR",              "N",
    "ERBB2", "Hs.446352",       2064L,   "ERBB2",              "N",
    "ESR1", "Hs.208124",       2099L,    "ESR1",              "N",
    "EXO1", "Hs.498248",       9156L,    "EXO1",              "N",
    "FGFR4", "Hs.165950",       2264L,   "FGFR4",              "N",
    "FOXA1", "Hs.163484",       3169L,   "FOXA1",              "N",
    "FOXC1", "Hs.348883",       2296L,   "FOXC1",              "N",
    "GPR160", "Hs.231320",      26996L,  "GPR160",              "N",
    "GRB7",  "Hs.86859",       2886L,    "GRB7",              "N",
    "KIF2C",  "Hs.69360",      11004L,   "KIF2C",              "N",
    "KNTC2", "Hs.414407",      10403L,   "NDC80",              "Y", #<<
    "KRT14", "Hs.654380",       3861L,   "KRT14",              "N",
    "KRT17",   "Hs.2785",       3872L,   "KRT17",              "N",
    "KRT5", "Hs.694210",       3852L,    "KRT5",              "N",
    "MAPT", "Hs.101174",       4137L,    "MAPT",              "N",
    "MDM2", "Hs.567303",       4193L,    "MDM2",              "N",
    "MELK", "Hs.184339",       9833L,    "MELK",              "N",
    "MIA", "Hs.646364",       8190L,     "MIA",              "N",
    "MKI67",  "Hs.80976",       4288L,   "MKI67",              "Y",
    "MLPH", "Hs.102406",      79083L,    "MLPH",              "N",
    "MMP11", "Hs.143751",       4320L,   "MMP11",              "N",
    "MYBL2", "Hs.179718",       4605L,   "MYBL2",              "N",
    "MYC", "Hs.202453",       4609L,     "MYC",              "N",
    "NAT1", "Hs.591847",          9L,    "NAT1",              "N",
    "ORC6L",  "Hs.49760",      23594L,   "ORC6L",              "N",
    "PGR", "Hs.368072",       5241L,     "PGR",              "N",
    "PHGDH", "Hs.487296",      26227L,   "PHGDH",              "N",
    "PTTG1", "Hs.350966",       9232L,   "PTTG1",              "Y",
    "RRM2", "Hs.226390",       6241L,    "RRM2",              "Y",
    "SFRP1", "Hs.695991",       6422L,   "SFRP1",              "N",
    "SLC39A6",  "Hs.79136",      25800L, "SLC39A6",              "N",
    "TMEM45B", "Hs.504301",     120224L, "TMEM45B",              "N",
    "TYMS", "Hs.592338",       7298L,    "TYMS",              "Y",
    "UBE2C",  "Hs.93002",      11065L,   "UBE2C",              "Y",
    "UBE2T",   "Hs.5199",      29089L,   "UBE2T",              "N"
  )
}