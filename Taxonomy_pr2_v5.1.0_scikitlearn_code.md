# Taxonomic classification for Functional Traits of Thalassiosira

Functional Traits of Thalassiosira samples from OME Run4 Sequencing Data taxonomic classification with pr2 v5.1.0 and scikit learn naive bayesian classifier.

## 18S v4: (target length 270 bp, but ranges from 200-416 bp)

```
# start qiime conda environment
conda activate qiime2-amplicon-2024.5
```

### Scikit learn pr2 ref db

Import data

```
# unzip pr2 files:
gunzip pr2_version_5.0.0_SSU_mothur.*

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path ref_db/pr2/pr2_version_5.1.0_SSU_mothur.fasta  \
  --output-path ref_db/pr2/pr2_v5.1.0_SSU_seqs.qza

# import taxonomy
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path ref_db/pr2/pr2_version_5.1.0_SSU_mothur.tax \
  --output-path ref_db/pr2/pr2_v5.1.0_SSU_tax.qza
```
Extract region of interest

```
# extact-reads to classify
qiime feature-classifier extract-reads \
  --i-sequences ref_db/pr2/pr2_v5.1.0_SSU_seqs.qza \
  --p-f-primer CCAGCASCYGCGGTAATTCC \
  --p-r-primer ACTTTCGTTCTTGATYR \
  --p-min-length 200 \
  --p-max-length 420 \
  --o-reads ref_db/pr2/pr2_v5.1.0_18Sv4_seqs.qza

qiime feature-table tabulate-seqs \
  --i-data ref_db/pr2/pr2_v5.1.0_18Sv4_seqs.qza \
  --o-visualization ref_db/pr2/pr2_v5.1.0_18Sv4_seqs.qzv
```
After sequences with primers are extracted:
135,169 seqs, 201-420 length, 374 mean length

Evaluate reference: (pr2 18S v4 only)
* Seqs in PR2: 309,663
* Seqs 18S v4: 144,121 seqs, 201-420 seq length, 373 mean length

Now extract any reads that were maybe missing primers

```
# First steps 
 qiime rescript extract-seq-segments \
    --i-input-sequences ref_db/pr2/pr2_v5.1.0_SSU_seqs.qza \
    --i-reference-segment-sequences ref_db/pr2/pr2_v5.1.0_18Sv4_seqs.qza \
    --p-perc-identity 0.8 \
    --p-min-seq-len 100 \
    --p-threads 5 \
    --o-extracted-sequence-segments ref_db/pr2/pr2_v5.1.0_18Sv4_seqs-01.qza \
    --o-unmatched-sequences ref_db/pr2/pr2_v5.1.0_18Sv4_unmatched_seqs-01.qza \
    --verbose 

qiime feature-table tabulate-seqs \
  --i-data ref_db/pr2/pr2_v5.1.0_18Sv4_seqs-01.qza \
  --o-visualization ref_db/pr2/pr2_v5.1.0_18Sv4_seqs-01.qzv

```
After first round of extract-seq-segments:
200,695 seqs, 186-507 length, 363 mean length

Now extract any reads that were maybe missing primers, Round 2.

```
 qiime rescript extract-seq-segments \
    --i-input-sequences ref_db/pr2/pr2_v5.1.0_SSU_seqs.qza \
    --i-reference-segment-sequences ref_db/pr2/pr2_v5.1.0_18Sv4_seqs-01.qza \
    --p-perc-identity 0.9 \
    --p-min-seq-len 100 \
    --p-threads 5 \
    --o-extracted-sequence-segments ref_db/pr2/pr2_v5.1.0_18Sv4_seqs-02.qza \
    --o-unmatched-sequences ref_db/pr2/pr2_v5.1.0_18Sv4_unmatched_seqs-02.qza \
    --verbose 

qiime feature-table tabulate-seqs \
  --i-data ref_db/pr2/pr2_v5.1.0_18Sv4_seqs-02.qza \
  --o-visualization ref_db/pr2/pr2_v5.1.0_18Sv4_seqs-02.qzv

```
After second round of extract-seq-segments:
212,382 seqs, 185-520 length, 364 mean length

Now extract any reads that were maybe missing primers, Round 3.

```
 qiime rescript extract-seq-segments \
    --i-input-sequences ref_db/pr2/pr2_v5.1.0_SSU_seqs.qza \
    --i-reference-segment-sequences ref_db/pr2/pr2_v5.1.0_18Sv4_seqs-02.qza \
    --p-perc-identity 0.9 \
    --p-min-seq-len 100 \
    --p-threads 6 \
    --o-extracted-sequence-segments ref_db/pr2/pr2_v5.1.0_18Sv4_seqs-03.qza \
    --o-unmatched-sequences ref_db/pr2/pr2_v5.1.0_18Sv4_unmatched_seqs-03.qza \
    --verbose 

# Compare counts of seqs
qiime feature-table tabulate-seqs \
  --i-data ref_db/pr2/pr2_v5.1.0_18Sv4_seqs-03.qza \
  --o-visualization ref_db/pr2/pr2_v5.1.0_18Sv4_seqs-03.qzv
```
After third round of extract-seq-segments:
214,116 seqs, 183-522 length, 365 mean length

* Note: if you're able to extract sequence segments from ~50-70% of your original query sequences, you likely have enough to stop the segment extraction iterations and move on to other downstream curation tasks. (From online rescript guide)
*After 2nd extract-seq-segments iteration, we have extracted ~70% of original sequences so stop here!*


Using second extract-seq-segments iteration:
Get taxa that remain:

```
qiime rescript filter-taxa \
    --i-taxonomy ref_db/pr2/pr2_v5.1.0_SSU_tax.qza \
    --m-ids-to-keep-file ref_db/pr2/pr2_v5.1.0_18Sv4_seqs-02.qza \
    --o-filtered-taxonomy ref_db/pr2/pr2_v5.1.0_18Sv4_tax-02.qza
```

### Dereplicate V4 reads with -uniq mode

Dereplicate the extracted reads, retaining identical sequence records that have differing taxonomies. 

```
qiime rescript dereplicate \
    --i-sequences ref_db/pr2/pr2_v5.1.0_18Sv4_seqs-02.qza  \
    --i-taxa ref_db/pr2/pr2_v5.1.0_18Sv4_tax-02.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences ref_db/pr2/pr2_v5.1.0_18Sv4_uniq-seqs.qza \
    --o-dereplicated-taxa ref_db/pr2/pr2_v5.1.0_18Sv4_uniq-tax.qza
```

### Evaluate Reference
```
qiime metadata tabulate \
    --m-input-file pr2_v5.1.0_18Sv4_seqs.qza \
    --o-visualization pr2_v5.1.0_18Sv4_seqs.qzv
    
qiime metadata tabulate \
    --m-input-file pr2_v5.1.0_18Sv4_uniq-seqs.qza \
    --o-visualization pr2_v5.1.0_18Sv4_uniq-seqs.qzv
    
qiime rescript evaluate-taxonomy \
    --i-taxonomies pr2_v5.1.0_18Sv4_uniq-tax.qza \
    --o-taxonomy-stats pr2_v5.1.0_18Sv4_uniq-tax.qzv

# code below did not work in my qiime environment
# qiime rescript evaluate-seqs \
#    --i-sequences pr2_v5.1.0_18Sv4_uniq-seqs.qza \
#    --p-kmer-lengths 8 4 2 \
#    --o-visualization pr2_v5.1.0_18Sv4_uniq-seqs-eval.qzv
```

212,382 seqs (not dereplicated)

129,520 seqs (Derep)

309,663 originally

### Train classifier

Use qiime2-2024.10.1

```
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads pr2_v5.1.0_18Sv4_uniq-seqs.qza \
  --i-reference-taxonomy pr2_v5.1.0_18Sv4_uniq-tax.qza \
  --o-classifier pr2_v5.1.0_18Sv4_uniq-classifier.qza

```

### Use classifier on Run4 sequencing data:

#### Use qiime2-amplicon-2024.10

```
# import fasta sequences
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path OME_Run4/Run4_18Sv4_ASVs.fasta \
  --output-path OME_Run4/Run4_18Sv4_ASVs.qza

# run classifier
qiime feature-classifier classify-sklearn \
   --i-read OME_Run4/Run4_18Sv4_ASVs.qza \
   --i-classifier pr2_v5.1.0_18Sv4_uniq-classifier.qza \
   --p-n-jobs 5 \
   --o-classification OME_Run4/Run4_18Sv4_tax_pr2.5.1.0_sklearn.qza

# export taxonomy
qiime tools export \
  --input-path OME_Run4/Run4_18Sv4_tax_pr2.5.1.0_sklearn.qza \
  --output-path OME_Run4/Run4_18Sv4_tax_pr2.5.1.0_sklearn.tsv

``` 