-- The bidirectional role of DTA mutational burden in MDS outcomes: a retrospective study and machine learning analysis
-------------------------------------- DATASET CREATION AND PROCESSING ------------------------------------------------

----> NOTE: IMPORT DATA FROM CSV FILES USING PGADMIN4 GUI (RIGHT AFTER EACH TABLE CREATION)


---------- From cBioPortal MDS dataset
-- Create table samples
DROP TABLE IF EXISTS final_project.table_samples;
CREATE TABLE final_project.table_samples
 (patient_id         CHARACTER VARYING,
 sample_id           CHARACTER VARYING,
 mds_type            CHARACTER VARYING,
 who_2016            CHARACTER VARYING,
 bm_blast            TEXT,
 pb_blast            TEXT,
 hbg                 TEXT,
 plt                 TEXT,
 wbc                 TEXT,
 anc                 TEXT,
 monocytes           TEXT,
 ringed_sideroblasts TEXT,
 cytogenetics        TEXT,
 cyto_ipssr          TEXT,
 complex_karyotype   TEXT,
 chr_tp53_locus      TEXT,
 cnacs_chrarm_loss   TEXT,
 cnacs_gene_loss     TEXT,
 cnacs_chrarm_gain   TEXT,
 cnacs_gene_gain     TEXT,
 cnacs_chrarm_upd    TEXT,
 cnacs_gene_upd      TEXT,
 flt3_itd            TEXT,
 mll_ptd             TEXT,
 ipssr               CHARACTER VARYING,
 ipssr_score         TEXT,
 ipmms               CHARACTER VARYING,
 ipmms_score         TEXT,
 tmb_nonsynonymous   TEXT
);

---> IMPORT DATA FROM CSV FILE USING PGADMIN4 GUI!

-- Create patients table
DROP TABLE IF EXISTS final_project.patients;
CREATE TABLE final_project.patients
 (patient_id CHARACTER VARYING,
 sex         CHARACTER VARYING,
 age         TEXT,
 os_months   TEXT,
 os_status   CHARACTER VARYING,
 lfs_months  TEXT,
 lfs_status  CHARACTER VARYING
);

----> IMPORT DATA FROM CSV FILE USING PGADMIN4 GUI!


-- Create mutations table
--> to do: convert NA, X and other symbols to blank
DROP TABLE IF EXISTS final_project.mutations;
CREATE TABLE final_project.mutations    
(hugo_symbol                    CHARACTER VARYING,
entrez_gene_id                  INTEGER,
center                          TEXT,
ncbi_build                      CHARACTER VARYING,
chromosome                      TEXT,
start_position                  INTEGER,
end_position                    INTEGER,
strand                          CHARACTER VARYING,
consequence                     CHARACTER VARYING,
variant_classification          CHARACTER VARYING,
variant_type                    CHARACTER VARYING,
reference_allele                CHARACTER VARYING,
tumor_seq_allele1               CHARACTER VARYING,
tumor_seq_allele2               CHARACTER VARYING,
dbsnp_rs                        CHARACTER VARYING,
db_snp_val_status               TEXT,
sample_id                       CHARACTER VARYING,
matched_norm_sample             TEXT,
matched_norm_seq_allele1        TEXT,
matched_norm_seq_allele2        TEXT,
tumor_validation_allele1        TEXT,
tumor_validation_allele2        TEXT,
match_norm_validation_allele1   TEXT,
match_norm_validation_allele2   TEXT,
verification_status             TEXT,
validation_status               TEXT,
mutation_status                 TEXT,
sequencing_phase                TEXT,
sequence_source                 TEXT,
validation_methods              TEXT,
score                           TEXT,
bam_file                        TEXT,
sequencer                       TEXT,
t_ref_count                     INTEGER,
t_alt_count                     INTEGER,
n_ref_count                     INTEGER,
n_alt_count                     INTEGER,
hgvsc                           CHARACTER VARYING,
hgvsp                           CHARACTER VARYING,
hgvsp_short                     CHARACTER VARYING,
transcript_id                   CHARACTER VARYING,
ref_seq                         CHARACTER VARYING,
protein_position                INTEGER,
codons                          CHARACTER VARYING,
exon_number                     TEXT,
cbp_driver                      CHARACTER VARYING,
cbp_driver_annotation           CHARACTER VARYING,
annotation_status               CHARACTER VARYING
);

-- IMPORT DATA FROM CSV FILE USING PGADMIN4 GUI

--> Pre-processing HGVS coding (Ensembl) column (split column by ':')
ALTER TABLE final_project.mutations -- add columns
  ADD COLUMN IF NOT EXISTS ensembl_short text,
  ADD COLUMN IF NOT EXISTS hgvs_coding  text;

UPDATE final_project.mutations
SET
  ensembl_short = btrim(split_part(hgvsc, ':', 1)), -- split left side)
  hgvs_coding   = NULLIF(btrim(split_part(hgvsc, ':', 2)), '') -- split (right side)
WHERE hgvsc IS NOT NULL;


--> adding VAF (variant allele frequency)
ALTER TABLE final_project.mutations
ADD COLUMN IF NOT EXISTS vaf NUMERIC;

UPDATE final_project.mutations
SET vaf = ROUND((t_alt_count::NUMERIC / (t_alt_count + t_ref_count)) * 100, 2);



---------- From ClinVar

----- ASXL1
DROP TABLE IF EXISTS final_project.clinvar_asxl1;
CREATE TABLE final_project.clinvar_asxl1
(ref_seq                    CHARACTER VARYING,
hugo_symbol                 CHARACTER VARYING,
hgvsp_short                 CHARACTER VARYING,
disease                     CHARACTER VARYING,
accession                   CHARACTER VARYING,
chromosome                  INTEGER,
start_position              TEXT,
grch38_chromosome           INTEGER,
grch38_position             TEXT,
clinvar_id                  INTEGER,
allele_id                   INTEGER,
dbsnp_id                    CHARACTER VARYING,
canonical_spdi              CHARACTER VARYING,
variant_type                CHARACTER VARYING,
consequence                 CHARACTER VARYING,
germline_classification     CHARACTER VARYING,
germline_date_eval          TEXT,
germline_review_status      CHARACTER VARYING,
somatic_impact              CHARACTER VARYING,
somatic_impact_date_eval    TEXT,
somatic_review_status       CHARACTER VARYING,
oncogenicity_classif        CHARACTER VARYING,
oncogenicity_date_eval      TEXT,
oncogenicity_review_status  CHARACTER VARYING
);

----> IMPORT DATA FROM CSV FILE USING PGADMIN4 GUI!


--> Pre-processing HGVS coding (Ensembl) column (split column by ':', followed by ' ') 
-- First chunck
ALTER TABLE final_project.clinvar_asxl1 -- add columns
  ADD COLUMN IF NOT EXISTS ref_seq1 text,
  ADD COLUMN IF NOT EXISTS hgvs_coding_tmp text;

UPDATE final_project.clinvar_asxl1
SET
  ref_seq1 = btrim(split_part(ref_seq, ':', 1)), -- split left side)
  hgvs_coding_tmp = NULLIF(btrim(split_part(ref_seq, ':', 2)), '') -- split (right side)
WHERE ref_seq IS NOT NULL;

-- Second chunck
ALTER TABLE final_project.clinvar_asxl1 -- add columns
  ADD COLUMN IF NOT EXISTS hgvs_coding_asxl1    text,
  ADD COLUMN IF NOT EXISTS protein text;

UPDATE final_project.clinvar_asxl1
SET
  hgvs_coding_asxl1 = btrim(split_part(hgvs_coding_tmp, ' ', 1)), -- split left side)
  protein = NULLIF(btrim(split_part(hgvs_coding_tmp, ' ', 2)), '') -- split (right side)
WHERE hgvs_coding_tmp IS NOT NULL;


----- DNMT3A
DROP TABLE IF EXISTS final_project.clinvar_dnmt3a;
CREATE TABLE final_project.clinvar_dnmt3a
(ref_seq                    CHARACTER VARYING,
hugo_symbol                 CHARACTER VARYING,
hgvsp_short                 CHARACTER VARYING,
disease                     CHARACTER VARYING,
accession                   CHARACTER VARYING,
chromosome                  INTEGER,
start_position              TEXT,
grch38_chromosome           INTEGER,
grch38_position             TEXT,
clinvar_id                  INTEGER,
allele_id                   INTEGER,
dbsnp_id                    CHARACTER VARYING,
canonical_spdi              CHARACTER VARYING,
variant_type                CHARACTER VARYING,
consequence                 CHARACTER VARYING,
germline_classification     CHARACTER VARYING,
germline_date_eval          TEXT,
germline_review_status      CHARACTER VARYING,
somatic_impact              CHARACTER VARYING,
somatic_impact_date_eval    TEXT,
somatic_review_status       CHARACTER VARYING,
oncogenicity_classif        CHARACTER VARYING,
oncogenicity_date_eval      TEXT,
oncogenicity_review_status  CHARACTER VARYING
);

----> IMPORT DATA FROM CSV FILE USING PGADMIN4 GUI!

--> Pre-processing HGVS coding (Ensembl) column (split column by ':', followed by ' ') 
-- First chunck
ALTER TABLE final_project.clinvar_dnmt3a -- add columns
  ADD COLUMN IF NOT EXISTS ref_seq1 text,
  ADD COLUMN IF NOT EXISTS hgvs_coding_tmp text;

UPDATE final_project.clinvar_dnmt3a
SET
  ref_seq1 = btrim(split_part(ref_seq, ':', 1)), -- split left side)
  hgvs_coding_tmp = NULLIF(btrim(split_part(ref_seq, ':', 2)), '') -- split (right side)
WHERE ref_seq IS NOT NULL;

-- Second chunck
ALTER TABLE final_project.clinvar_dnmt3a -- add columns
  ADD COLUMN IF NOT EXISTS hgvs_coding_dnmt3a    text,
  ADD COLUMN IF NOT EXISTS protein text;

UPDATE final_project.clinvar_dnmt3a
SET
  hgvs_coding_dnmt3a = btrim(split_part(hgvs_coding_tmp, ' ', 1)), -- split left side)
  protein = NULLIF(btrim(split_part(hgvs_coding_tmp, ' ', 2)), '') -- split (right side)
WHERE hgvs_coding_tmp IS NOT NULL;


----- TET2
DROP TABLE IF EXISTS final_project.clinvar_tet2;
CREATE TABLE final_project.clinvar_tet2
(ref_seq                    CHARACTER VARYING,
hugo_symbol                 CHARACTER VARYING,
hgvsp_short                 CHARACTER VARYING,
disease                     CHARACTER VARYING,
accession                   CHARACTER VARYING,
chromosome                  INTEGER,
start_position              TEXT,
grch38_chromosome           INTEGER,
grch38_position             TEXT,
clinvar_id                  INTEGER,
allele_id                   INTEGER,
dbsnp_id                    CHARACTER VARYING,
canonical_spdi              CHARACTER VARYING,
variant_type                CHARACTER VARYING,
consequence                 CHARACTER VARYING,
germline_classification     CHARACTER VARYING,
germline_date_eval          TEXT,
germline_review_status      CHARACTER VARYING,
somatic_impact              CHARACTER VARYING,
somatic_impact_date_eval    TEXT,
somatic_review_status       CHARACTER VARYING,
oncogenicity_classif        CHARACTER VARYING,
oncogenicity_date_eval      TEXT,
oncogenicity_review_status  CHARACTER VARYING
);

----> IMPORT DATA FROM CSV FILE USING PGADMIN4 GUI!

--> Pre-processing HGVS coding (Ensembl) column (split column by ':', followed by ' ') 
-- First chunck
ALTER TABLE final_project.clinvar_tet2 -- add columns
  ADD COLUMN IF NOT EXISTS ref_seq1 text,
  ADD COLUMN IF NOT EXISTS hgvs_coding_tmp text;

UPDATE final_project.clinvar_tet2
SET
  ref_seq1 = btrim(split_part(ref_seq, ':', 1)), -- split left side)
  hgvs_coding_tmp = NULLIF(btrim(split_part(ref_seq, ':', 2)), '') -- split (right side)
WHERE ref_seq IS NOT NULL;

-- Second chunck
ALTER TABLE final_project.clinvar_tet2 -- add columns
  ADD COLUMN IF NOT EXISTS hgvs_coding_tet2    text,
  ADD COLUMN IF NOT EXISTS protein text;

UPDATE final_project.clinvar_tet2
SET
  hgvs_coding_tet2 = btrim(split_part(hgvs_coding_tmp, ' ', 1)), -- split left side)
  protein = NULLIF(btrim(split_part(hgvs_coding_tmp, ' ', 2)), '') -- split (right side)
WHERE hgvs_coding_tmp IS NOT NULL;


---------- From gnomAD dataset

------ ASXL1
DROP TABLE IF EXISTS final_project.gnomad_asxl1;
CREATE TABLE final_project.gnomad_asxl1
(gnomad_id                  CHARACTER VARYING,
chromosome                  CHARACTER VARYING,
start_position              INTEGER,
dbsnp_id                    CHARACTER VARYING,
reference_allele            CHARACTER VARYING,
alternate_allele            CHARACTER VARYING,
sequence_source             CHARACTER VARYING,
exome_filters               CHARACTER VARYING,
genome_filters              CHARACTER VARYING,
transcript                  CHARACTER VARYING,
hvgs_consequence            CHARACTER VARYING, 
hgvsp                       CHARACTER VARYING,
hgvs_coding                 CHARACTER VARYING,
consequence                 CHARACTER VARYING,
germline_classification     CHARACTER VARYING,
clinvar_id                  CHARACTER VARYING,
flags                       CHARACTER VARYING,
allele_count                INTEGER,
allele_numer                INTEGER,
allele_frequency            NUMERIC,
homozygote                  INTEGER,
heterozygote                INTEGER,
filters_joint               CHARACTER VARYING,
group_max_faf_group         CHARACTER VARYING,
group_max_faf_freq          NUMERIC,
cadd                        NUMERIC,
revel_max                   NUMERIC,
spliceai_ds_max             NUMERIC,
pangolin_largest_ds         NUMERIC,
phylop                      NUMERIC,
sift_max                    NUMERIC,
polyphen_max                NUMERIC
);

----> IMPORT DATA FROM CSV FILE USING PGADMIN4 GUI!

----- DNMT3A
DROP TABLE IF EXISTS final_project.gnomad_dnmt3a;
CREATE TABLE final_project.gnomad_dnmt3a
(gnomad_id                  CHARACTER VARYING,
chromosome                  CHARACTER VARYING,
start_position              INTEGER,
dbsnp_id                    CHARACTER VARYING,
reference_allele            CHARACTER VARYING,
alternate_allele            CHARACTER VARYING,
sequence_source             CHARACTER VARYING,
exome_filters               CHARACTER VARYING,
genome_filters              CHARACTER VARYING,
transcript                  CHARACTER VARYING,
hvgs_consequence            CHARACTER VARYING, 
hgvsp                       CHARACTER VARYING,
hgvs_coding                 CHARACTER VARYING,
consequence                 CHARACTER VARYING,
germline_classification     CHARACTER VARYING,
clinvar_id                  CHARACTER VARYING,
flags                       CHARACTER VARYING,
allele_count                INTEGER,
allele_numer                INTEGER,
allele_frequency            NUMERIC,
homozygote                  INTEGER,
heterozygote                INTEGER,
filters_joint               CHARACTER VARYING,
group_max_faf_group         CHARACTER VARYING,
group_max_faf_freq          NUMERIC,
cadd                        NUMERIC,
revel_max                   NUMERIC,
spliceai_ds_max             NUMERIC,
pangolin_largest_ds         NUMERIC,
phylop                      NUMERIC,
sift_max                    NUMERIC,
polyphen_max                NUMERIC
);

----> IMPORT DATA FROM CSV FILE USING PGADMIN4 GUI!


------ TET2
DROP TABLE IF EXISTS final_project.gnomad_tet2;
CREATE TABLE final_project.gnomad_tet2
(gnomad_id                  CHARACTER VARYING,
chromosome                  CHARACTER VARYING,
start_position              INTEGER,
dbsnp_id                    CHARACTER VARYING,
reference_allele            CHARACTER VARYING,
alternate_allele            CHARACTER VARYING,
sequence_source             CHARACTER VARYING,
exome_filters               CHARACTER VARYING,
genome_filters              CHARACTER VARYING,
transcript                  CHARACTER VARYING,
hvgs_consequence            CHARACTER VARYING, 
hgvsp                       CHARACTER VARYING,
hgvs_coding                 CHARACTER VARYING,
consequence                 CHARACTER VARYING,
germline_classification     CHARACTER VARYING,
clinvar_id                  CHARACTER VARYING,
flags                       CHARACTER VARYING,
allele_count                INTEGER,
allele_numer                INTEGER,
allele_frequency            NUMERIC,
homozygote                  INTEGER,
heterozygote                INTEGER,
filters_joint               CHARACTER VARYING,
group_max_faf_group         CHARACTER VARYING,
group_max_faf_freq          NUMERIC,
cadd                        NUMERIC,
revel_max                   NUMERIC,
spliceai_ds_max             NUMERIC,
pangolin_largest_ds         NUMERIC,
phylop                      NUMERIC,
sift_max                    NUMERIC,
polyphen_max                NUMERIC
);

----> IMPORT DATA FROM CSV FILE USING PGADMIN4 GUI!



--------------- Defining cohort -------------------
-- Number of samples per patient
SELECT patient_id, COUNT(sample_id)
FROM final_project.table_samples
GROUP BY patient_id;


---> Conclusion: although in real-world one patient can have multiple samples, for this study, each patient contributes with only one sample

---- Cohort: patient with MDS data extracted from cBioPortal
-- Exclusion criteria: no benign/ likely benign DTA variants per ClinVar
-- Primary outcomes: mostly in patients table - os_months, os_status (correlation with DTA mutation burden)
-- Secondary outcomes: mostly in patients table - lfs_months, lfs_status (correlation with DTA mutation burden)



---------------------- EXCLUSION CRITERIA ------------------------------------
----------------------------------------------------------------------------
-- Create cleaned samples table excluding benign/ likely benign variants

DROP TABLE IF EXISTS temp_cleaned_samples;
CREATE TEMP TABLE temp_cleaned_samples AS

WITH benign_only AS (
    SELECT DISTINCT s.sample_id
    FROM final_project.table_samples s
    JOIN final_project.mutations m
      ON s.sample_id = m.sample_id
    LEFT JOIN final_project.clinvar_asxl1 c1
      ON m.hugo_symbol = 'ASXL1'
     AND m.hgvs_coding = c1.hgvs_coding_asxl1
    LEFT JOIN final_project.clinvar_dnmt3a c2
      ON m.hugo_symbol = 'DNMT3A'
     AND m.hgvs_coding = c2.hgvs_coding_dnmt3a
    LEFT JOIN final_project.clinvar_tet2 c3
      ON m.hugo_symbol = 'TET2'
     AND m.hgvs_coding = c3.hgvs_coding_tet2
    WHERE m.hugo_symbol IN ('ASXL1','DNMT3A','TET2')
      AND (c1.germline_classification ILIKE '%benign%'
         OR c2.germline_classification ILIKE '%benign%'
         OR c3.germline_classification ILIKE '%benign%'))

SELECT s.*
FROM final_project.table_samples s
WHERE NOT EXISTS (
    SELECT 1
    FROM benign_only b
    WHERE b.sample_id = s.sample_id);

SELECT * FROM temp_cleaned_samples; -- check as needed




---------------------- JOININGS AND PROCESSING ----------------------
---------------------------------------------------------------------

---------- Step 1: defining type of mutation for each gene (truncating) from table_mutations

-- ASXL1
DROP TABLE IF EXISTS tmp_asxl1;
CREATE TEMP TABLE tmp_asxl1 AS
SELECT DISTINCT m.sample_id, m.exon_number, m.cbp_driver_annotation, m.hgvs_coding, m.hugo_symbol,m.vaf,
CASE WHEN m.consequence ILIKE ANY (ARRAY['%splice%','%stop%','%start%','%frameshift%','%lost%','%UTR%'])
  THEN 1
  ELSE 0
  END AS type_of_mutation
FROM final_project.mutations m
WHERE m.hugo_symbol = 'ASXL1';

-- DNMT3A
DROP TABLE IF EXISTS tmp_dnmt3a;
CREATE TEMP TABLE tmp_dnmt3a AS
SELECT DISTINCT m.sample_id, m.exon_number, m.cbp_driver_annotation, m.hgvs_coding, m.hugo_symbol, m.vaf, 
CASE WHEN m.consequence ILIKE ANY (ARRAY['%splice%','%stop%','%start%','%frame%','%lost%','%UTR%'])
  THEN 1
  ELSE 0
  END AS type_of_mutation
FROM final_project.mutations m
WHERE m.hugo_symbol = 'DNMT3A';


-- TET2
DROP TABLE IF EXISTS tmp_tet2;
CREATE TEMP TABLE tmp_tet2 AS
SELECT DISTINCT m.sample_id, m.exon_number, m.cbp_driver_annotation, m.hgvs_coding, m.hugo_symbol, m.vaf,
CASE WHEN m.consequence ILIKE ANY (ARRAY['%splice%','%stop%','%start%','%frame%','%lost%','%UTR%'])
  THEN 1
  ELSE 0
  END AS type_of_mutation
FROM final_project.mutations m
WHERE m.hugo_symbol = 'TET2';


----------- Step 2: Union function applied to mutation tables with DNMT3A, TET2 and ASXL1 data

DROP TABLE IF EXISTS tmp_union;
CREATE TEMP TABLE tmp_union AS
SELECT sample_id, exon_number, cbp_driver_annotation, hgvs_coding, type_of_mutation, hugo_symbol, vaf, 'ASXL1'::text as mutation
FROM tmp_asxl1
UNION ALL
SELECT sample_id, exon_number, cbp_driver_annotation, hgvs_coding, type_of_mutation, hugo_symbol, vaf, 'DNMT3A'::text as mutation
FROM tmp_dnmt3a
UNION ALL
SELECT sample_id, exon_number, cbp_driver_annotation, hgvs_coding, type_of_mutation, hugo_symbol, vaf, 'TET2'::text as mutation
FROM tmp_tet2;

-- SELECT * FROM tmp_union; (check as needed)


----------- Step 3: Joining ClinVar and GnomAD datasets per each gene (excluding benign/ likely benign variants)

----------- ASXL1
DROP TABLE IF EXISTS tmp_clinvar_gnomad_asxl1;
CREATE TEMP TABLE tmp_clinvar_gnomad_asxl1 AS
SELECT
  btrim(regexp_replace(c.hugo_symbol,   '[|-].*$', '')) AS hugo_symbol, 
  g.hgvs_coding,
  CASE WHEN c.germline_classification ILIKE '%pathogenic%' 
        THEN 1 ELSE 0  
        END AS pathogenic_asxl1,
  g.revel_max AS revel_asxl1,
  g.spliceai_ds_max AS spliceai_asxl1,
  g.cadd AS cadd_asxl1,
  g.phylop AS phylop_asxl1
FROM final_project.gnomad_asxl1 g
JOIN final_project.clinvar_asxl1 c ON g.hgvs_coding = c.hgvs_coding_asxl1
WHERE c.germline_classification NOT ILIKE '%benign%';


----------- DNMT3A
DROP TABLE IF EXISTS tmp_clinvar_gnomad_dnmt3a;
CREATE TEMP TABLE tmp_clinvar_gnomad_dnmt3a AS
SELECT
  btrim(regexp_replace(c.hugo_symbol,   '[|-].*$', '')) AS hugo_symbol, 
  g.hgvs_coding,
  CASE WHEN c.germline_classification ILIKE '%pathogenic%' 
    THEN 1 ELSE 0  
    END AS pathogenic_dnmt3a,
  g.revel_max AS revel_dnmt3a,
  g.spliceai_ds_max AS spliceai_dnmt3a,
  g.cadd AS cadd_dnmt3a,
  g.phylop AS phylop_dnmt3a
FROM final_project.gnomad_dnmt3a g
JOIN final_project.clinvar_dnmt3a c ON g.hgvs_coding = c.hgvs_coding_dnmt3a
WHERE c.germline_classification NOT ILIKE '%benign%';



----------- TET2
DROP TABLE IF EXISTS tmp_clinvar_gnomad_tet2;
CREATE TEMP TABLE tmp_clinvar_gnomad_tet2 AS
SELECT
  btrim(regexp_replace(c.hugo_symbol,   '[|-].*$', '')) AS hugo_symbol, 
  g.hgvs_coding,
  CASE WHEN c.germline_classification ILIKE '%pathogenic%' 
    THEN 1 ELSE 0  
    END AS pathogenic_tet2,
  g.revel_max AS revel_tet2,
  g.spliceai_ds_max AS spliceai_tet2,
  g.cadd AS cadd_tet2,
  g.phylop AS phylop_tet2
FROM final_project.gnomad_tet2 g
JOIN final_project.clinvar_tet2 c ON g.hgvs_coding = c.hgvs_coding_tet2
WHERE c.germline_classification NOT ILIKE '%benign%';


----------- Step 4: Apply exclusion criteria to tmp_union table (remove benign/ likely benign variants)


----------- Step 5: Joining Union Mutation table with GnomAD/ClinVar

DROP TABLE IF EXISTS tmp_mutation_clinvar_gnomad;
CREATE TEMP TABLE tmp_mutation_clinvar_gnomad AS
SELECT
  m.sample_id,
  -- Presence flags (0/1) -> "does this gene appear at all?""
MAX(CASE WHEN mutation = 'ASXL1'  THEN 1 ELSE 0 END) AS asxl1,
MAX(CASE WHEN mutation = 'DNMT3A' THEN 1 ELSE 0 END) AS dnmt3a,
MAX(CASE WHEN mutation = 'TET2'   THEN 1 ELSE 0 END) AS tet2,

-- Per-gene mutation counts -> "how many mutaitons in each gene (per-gene counts)"
COUNT(*) FILTER (WHERE mutation='ASXL1')  AS asxl1_count,
COUNT(*) FILTER (WHERE mutation='DNMT3A') AS dnmt3a_count,
COUNT(*) FILTER (WHERE mutation='TET2')   AS tet2_count,

-- Total DTA mutation count -> "how many total DTA mutations per sample_id, even if it is in the same gene"
COUNT(*) FILTER (WHERE mutation IN ('ASXL1','DNMT3A','TET2')) AS n_dta,

-- Number of different DTA genes mutated (0–3)
( MAX(CASE WHEN mutation='ASXL1'  THEN 1 ELSE 0 END)
+ MAX(CASE WHEN mutation='DNMT3A' THEN 1 ELSE 0 END)
+ MAX(CASE WHEN mutation='TET2'   THEN 1 ELSE 0 END) ) AS n_different_dta,

-- Presence flags (0/1) -> "is there any truncating variant?""
MAX(CASE WHEN type_of_mutation = 1 THEN 1 ELSE 0 END) AS truncating_variant,

-- Per-gene mutation counts -> "how many truncating variant in each gene (per-gene counts)"
COUNT(*) FILTER (WHERE mutation='ASXL1' AND type_of_mutation = 1)  AS asxl1_truncating_variant, 
COUNT(*) FILTER (WHERE mutation='DNMT3A' AND type_of_mutation = 1) AS dnmt3a_truncating_variant,
COUNT(*) FILTER (WHERE mutation='TET2' AND type_of_mutation = 1)  AS tet2_truncating_variant,

-- Total DTA mutation count -> "how many total DTA truncating variant per sample_id, even if it is in the same gene"
COUNT(*) FILTER (WHERE type_of_mutation = 1) AS n_truncating_variant,

-- Pathogenic sums based on ClinVar
COUNT(*) FILTER (WHERE pathogenic_asxl1=1)  AS pathogenic_asxl1,
COUNT(*) FILTER (WHERE pathogenic_dnmt3a=1)  AS pathogenic_dnmt3a,
COUNT(*) FILTER (WHERE pathogenic_tet2=1)   AS pathogenic_tet2,


 -- conditional averages: if >1 non-NULL, SUM/COUNT; else the single value
CASE
    WHEN COUNT(m.vaf) > 1
      THEN ROUND(SUM(m.vaf)::numeric / COUNT(m.vaf), 2)
    ELSE MAX(m.vaf)
  END AS vaf_avg,

  -- Cross-source totals (no outer COALESCE)
-- REVEL
ROUND(
  SUM(COALESCE(a.revel_asxl1, d.revel_dnmt3a, t.revel_tet2))::numeric
  / NULLIF(COUNT(COALESCE(a.revel_asxl1, d.revel_dnmt3a, t.revel_tet2)), 0)
, 2) AS revel_avg,

-- CADD
ROUND(
  SUM(COALESCE(a.cadd_asxl1, d.cadd_dnmt3a, t.cadd_tet2))::numeric
  / NULLIF(COUNT(COALESCE(a.cadd_asxl1, d.cadd_dnmt3a, t.cadd_tet2)), 0)
, 2) AS cadd_avg,

-- PhyloP 
ROUND(
  SUM(COALESCE(a.phylop_asxl1, d.phylop_dnmt3a, t.phylop_tet2))::numeric
  / NULLIF(COUNT(COALESCE(a.phylop_asxl1, d.phylop_dnmt3a, t.phylop_tet2)), 0)
, 2) AS phylop_avg

FROM tmp_union m
  LEFT JOIN tmp_clinvar_gnomad_asxl1 a
    ON m.hugo_symbol = a.hugo_symbol AND m.hgvs_coding = a.hgvs_coding
  LEFT JOIN tmp_clinvar_gnomad_dnmt3a d
    ON m.hugo_symbol = d.hugo_symbol AND m.hgvs_coding = d.hgvs_coding
  LEFT JOIN tmp_clinvar_gnomad_tet2 t
    ON m.hugo_symbol = t.hugo_symbol AND m.hgvs_coding = t.hgvs_coding
GROUP BY m.sample_id;

SELECT * from tmp_mutation_clinvar_gnomad;



----------- Step 5: Final cohort table joining patients, samples and mutation data


DROP TABLE IF EXISTS final_project.cohort_ml_v4;
CREATE TABLE final_project.cohort_ml_v4 AS
SELECT
  p.patient_id,
  m.sample_id,
  p.os_months,
  p.os_status,
  p.lfs_months,
  p.lfs_status,
  p.sex,
  p.age,
  s.mds_type,
  s.who_2016,
  s.bm_blast,
  s.pb_blast,
  s.hbg,
  s.plt,
  s.wbc,
  s.anc,
  s.monocytes,
  s.cyto_ipssr,
  s.complex_karyotype,
CASE -- for TP53 mutation alteration status as needed
   WHEN EXISTS (
        SELECT 1
        FROM final_project.mutations mut
        WHERE mut.sample_id = s.sample_id
          AND mut.hugo_symbol = 'TP53'
   )
   OR LOWER(s.chr_tp53_locus) IN ('cnloh','del','iso17q')
   THEN 1
   ELSE 0
END AS tp53_altered,
  s.flt3_itd,
  s.mll_ptd,
  s.ipssr,
  s.ipssr_score,
  s.ipmms,
  s.ipmms_score,  
  s.tmb_nonsynonymous,
  m.asxl1,
  m.dnmt3a,
  m.tet2,
  m.asxl1_count,
  m.dnmt3a_count,
  m.tet2_count,
  m.n_dta,
  m.n_different_dta,
  m.truncating_variant,
  m.asxl1_truncating_variant,
  m.dnmt3a_truncating_variant,
  m.tet2_truncating_variant,
  m.n_truncating_variant,
  m.pathogenic_asxl1,
  m.pathogenic_dnmt3a,
  m.pathogenic_tet2,
  m.vaf_avg,
  m.revel_avg,
  m.cadd_avg,
  m.phylop_avg
  FROM temp_cleaned_samples s
  JOIN final_project.patients p ON s.patient_id = p.patient_id -- cleaned samples table after exclusion criteria
  LEFT JOIN tmp_mutation_clinvar_gnomad m ON s.sample_id = m.sample_id; -- change to inner join if the goal is only DTA mutation population


SELECT * FROM final_project.cohort_ml_v4; -- export as csv as "mds_dta_cohort_all_mds.csv"


-------------------------------------- DEMOGRAPHICS  ------------------------------------------------
-- Number of patients in cBioPortal MDS dataset
SELECT DISTINCT COUNT(patient_id) FROM final_project.table_samples;
-- N = 3,323

-- Number of patients in the final cohort
SELECT DISTINCT COUNT(patient_id) FROM final_project.cohort_ml_v4;
-- N = 3,300 patients after exclusion criteria (benign/ likely benign DTA variants removed)

-- Number of individuals with mutation in at least one DTA gene
SELECT COUNT(DISTINCT sample_id) -- since one sample per patient
FROM final_project.mutations
WHERE hugo_symbol IN ('ASXL1', 'TET2', 'DNMT3A');
-- N = 2,000 samples with mutation in at least one DTA gene (therefore, 1,323 patients without DTA mutation)

-- Number of individuals with at least one DTA mutation in the studied cohort after applying exclusion criteria
SELECT COUNT(DISTINCT c.patient_id)
FROM final_project.cohort_ml_v4 c
WHERE c.n_dta >= 1;
-- N = 1,977 patients with at least one DTA mutation after exclusion criteria 

-- Number of individuals without mutation in any DTA gene
SELECT COUNT(DISTINCT s.patient_id)
FROM final_project.table_samples s
WHERE NOT EXISTS (
    SELECT 1
    FROM final_project.mutations m
    WHERE m.sample_id = s.sample_id
      AND m.hugo_symbol IN ('ASXL1', 'TET2', 'DNMT3A')
);
-- N = 1,323 patients without mutation in any DTA gene as expected


-- Number of individuals with benign/ likely benign DTA mutations
-- COUNT EXCLUSION CRITERIA - BENIGN/ LIKELY BENIGN

SELECT COUNT(DISTINCT s.patient_id)
FROM final_project.table_samples s
JOIN final_project.mutations m
  ON s.sample_id = m.sample_id
LEFT JOIN final_project.clinvar_asxl1 c1
  ON m.hugo_symbol = 'ASXL1'
 AND m.hgvs_coding = c1.hgvs_coding_asxl1
LEFT JOIN final_project.clinvar_dnmt3a c2
  ON m.hugo_symbol = 'DNMT3A'
 AND m.hgvs_coding = c2.hgvs_coding_dnmt3a
LEFT JOIN final_project.clinvar_tet2 c3
  ON m.hugo_symbol = 'TET2'
 AND m.hgvs_coding = c3.hgvs_coding_tet2
WHERE m.hugo_symbol IN ('ASXL1','DNMT3A','TET2')
  AND (
        (c1.germline_classification ILIKE '%benign%')
     OR (c2.germline_classification ILIKE '%benign%')
     OR (c3.germline_classification ILIKE '%benign%')
  );

