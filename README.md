# Poptools

- **scanning_signalsnp**

  This program is used to scan the significant signals of SNPs in gemma analysis result. It will extract the SNPs with p-value less than the threshold and their annotation from snpEff annotation file. The output will be a table with SNP, p-value, Pve, -Log, and their annotation.

  Run `scanning_signalsnp --help` to show the program's usage guide.
  ```bash
   Usage:./scanning_signalsnp -g <gemma_output_file> -s <snpEff_annotation_file> -n <sample_number> [-t <threshold>] [-pre <prefix>] [-o <output_path>] [-h]
   Required options:
      -g, --gemma  Intput the result of gemma analysis
      -s, --snpAnn    Intput the result of snpEff annotation
      -n, --number  The sample number in gemma model
   Optional options:
      -t, --threshold  The threshold for p-value. default: 0.05/total_snps
      -pre, --prefix  Prefix of the output
      -o, --output The output path
      -h, --help      Display this help message
  ```


  **Example:**
  ```bash
  ./scanning_signalsnp -g ./gemma_output.accoc.txt -s ./snpEff_annotation.txt -n 1000 -t 7 -pre test -o ./output_path
  ```
  **The log will show the following information:**
  ```bash
   [2024-06-22 20:01:02] INFO: Gemma file: gemma_output.txt
   [2024-06-22 20:01:02] INFO: SnpEff annotation file: ./snpEff_annotation.txt
   [2024-06-22 20:01:02] INFO: Prefix: test
   [2024-06-22 20:01:02] INFO: Output file: /output_path
   [2024-06-22 20:01:02] INFO: Sample number: 1000
   [2024-06-22 20:01:02] INFO: Total snps: 1000000
   [2024-06-22 20:01:02] INFO: The threshold for p-value (-log10) is: 7.000000
   [2024-06-22 20:01:08] INFO: The number of signal snps found is 9
   -------------  ------------  ------------  ---------
   SNP            P-value       Pve           -Log
   -------------  ------------  ------------  ---------
   chr6:4814333 | 1.841745e-08 | 0.037303 | 7.734770
   chr6:4814334 | 1.841745e-08 | 0.037303 | 7.734770
   chr6:4819462 | 1.841745e-08 | 0.037303 | 7.734770
   chr6:4826854 | 7.037685e-10 | 0.045325 | 9.152570
   chr6:4827971 | 3.024143e-08 | 0.036078 | 7.519398
   chr6:4827976 | 1.730176e-08 | 0.037428 | 7.761909
   chr6:5086022 | 7.953511e-08 | 0.033730 | 7.099441
   chr7:8325427 | 8.009874e-08 | 0.033715 | 7.096375
   chr12:1591044 | 5.179233e-08 | 0.034843 | 7.285735
   -------------  ------------  ------------  ---------
   [2024-06-22 20:01:08] INFO: getting the annotation of signal snps from SnpEff annotation file ...
   [2024-06-22 20:01:08] INFO: The result file is ./output_path/test.scanning_signalsnp.txt
   [2024-06-22 20:01:08] INFO: INFO: Done
  ```

  **The format of gemma_output.accoc.txt is:**  
  note: Currently the script only applies to the results of gemma's LMM. The first line must be the same as the example below.
  ```
   chr     rs      ps      n_miss  allele1 allele0 af      beta    se      logl_H1 l_remle p_wald
   1       chr1:12363      12363   0       T       G       0.081   2.045538e+00    1.545180e+00    -9.536663e+02   2.022674e+00    1.867041e-01
  ```