* mayo.tcx.RData
  variable name:  mtx.tcx
    15160 genes:  "A1BG" "A1BG-AS1" "A2M" "A2M-AS1"  "A2ML1"    "A4GALT"
    264 samples:  "11492_TCX" "6810_TCX"  "1046_TCX"  "1924_TCX"  "1926_TCX"  "6913_TCX" ...
        fivenum:  -5.00730562 -0.66780069  0.00698823  0.68114114  4.31540280
          stats:  31M  (20 jul 2017)
  from /local/Cory/Alzheimers/synapse.windsorized/mayo.tcx.RData


* rosmap.fcx.RData
  variable name: mtx.fcx
          stats: 38M (08 oct 2017)
    15538 genes: "A1BG" "A1BG-AS1" "A2M" "A2M-AS1" "A2ML1" "A4GALT" ...
    638 samples: "s01_120405" "s02_120405" "s03_120405" "s04_120405" "s05_120405" "s07_120410"

* SCH_11923_B01_GRM_WGS_2017-04-27_5.recalibrated_variants.vcf.gz, .tbi
  variable name: NA
          stats: 3.5G, (27 mar 2017)
           from: https://www.synapse.org/#!Synapse:syn10995976
           rows: 5:88010133_A/T, 5:88010224_G/A, 5:88010323_G/G, 5:88010829_C/T, ...
        samples: S1005_TCX S1010_TCX S1015_TCX S1019_TCX S1029_TCX ...
          notes: hg19

* tbl.snp.hg38.score-ref-alt.RData
  variable name: tbl.snp
          stats: 7k (06 dec 2017)
           from: 020516_TableForCorySeth_AD_eQTL_Loci.xlsx
           rows: 155 snps,  "rs7721099"  "rs7703782"  "rs10514301" "rs10067451" "rs383883"   "rs13162708"
        columns: "rsid" "chrom", "pos", "score", "iupac", "ref", "alt", "A1", "CER_Beta",
	         "CER_P", "TX_Beta", "TX_P", "IGAP_A1", "IGAP_OR", "IGAP_Pvalue", "RegulomeDB",
		 "Rsquared_rs254776", "Dprime_rs254776"
          notes: lifted over to hg38 from spreadsheet's hg19

* mef2c.tf.5kb.RData
  4733 Dec 14 11:46

