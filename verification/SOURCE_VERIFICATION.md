# Model inputs and their sources

Published values used by the model (R/params.R and R/inputs.R), each checked against the primary publication or
prescribing information, with its source and page and its use in the model. Each section below becomes a block of the
Inputs_and_sources sheet of Supplementary Data 1 (scripts/12_article_tables.R).

## Fixed inputs

| Input | Value used | Source |
|---|---|---|
| Attrition: no endocrine-based second line after first-line progression at 6 months or later | 57% if progression at 6 to 12 months, 29% if after 12 months | Varnier 2025 (ESME-MBC), Suppl. Table A.2: second line by first-line PFS (75/131 and 75/257 not endocrine-based) |
| Ribociclib 600 mg: CYP3A substrate AUC ratio (sets the everolimus interaction) | 5.2 | ribociclib (Kisqali) prescribing information p30: midazolam AUC 3.8-fold at 400 mg (measured), 5.2-fold at 600 mg (predicted); Bardia 2021 (TRINITI-1): everolimus troughs 2-3-fold (ribociclib 250 mg) and 4-5-fold (300 mg) above expected |
| Palbociclib: CYP3A substrate AUC ratio | 1.61 | palbociclib (Ibrance) prescribing information: midazolam AUC +61% |
| Everolimus fraction of clearance by CYP3A4/P-gp (fm) | 0.933 | everolimus (Afinitor) prescribing information: ketoconazole raises everolimus AUC 15-fold (fm at least 1 - 1/15) |
| Palbociclib stopped for adverse events (PALOMA-2) | 7.4% | Finn 2016 (PALOMA-2) p7: 33/444 permanently stopped palbociclib |
| Median CDK4/6 inhibitor treatment duration in the toxicity calibration | palbociclib 603 d (19.8 months); ribociclib 396 d (13.0 months); abemaciclib 460 d (15.1 months) | palbociclib prescribing information p8; Hortobagyi 2016 (MONALEESA-2) p4 (median exposure in the report that gives the toxicity rates used); abemaciclib (Verzenio) prescribing information p13 |
| MONARCH 3 grade 2 or higher diarrhoea | 36.7% (calibration target, with grade 3 = 9.5%) | Goetz 2017 p4: grade 2 27.2% + grade 3 9.5% |
| MONARCH 3 grade 3-4 neutropenia | 21.1% | Goetz 2017 p1 (primary analysis, adverse-event based, as for PALOMA-2 and MONALEESA-2) |
| PALOMA-2 / MONALEESA-2 grade 3-4 neutropenia | 66.4% / 59.3% | Finn 2016 Table 2; Hortobagyi 2016 Table 3 |
| Tumour assessment schedule in the trial emulations | as published for each trial (section Tumour assessment schedules); model default where the report gives none | methods sections of the trial reports |
| TRINITI-1 ribociclib schedule | continuous (28 of 28 days) | Bardia 2021 pp. 1, 3: once daily, days 1-28 |
| postMONARCH and MAINTAIN eligibility | no minimum duration of prior CDK4/6 inhibitor therapy | Kalinsky 2025 p2; Kalinsky 2023 p2 |
| MAINTAIN comparator | own fulvestrant arm (MAINTAIN schedule and eligibility) | Kalinsky 2023 |
| evERA control-arm treatment duration (toxicity calibration) | 140 d (4.6 months) | Mayer 2026 (evERA) p7 |
| Mouthwash scenarios | no mouthwash x1.63 / x1.86 (BOLERO-2 vs evERA control); strict prophylaxis x0.47 / x0.62 (SWISH vs evERA control) | BOLERO-2 (Rugo 2014; everolimus prescribing information), SWISH (Rugo 2017), evERA (Mayer 2026) |

## Calibration targets: PALOMA-2, BOLERO-2, evERA control arm, acelERA and MONALEESA-2

| Target | Published value | Source | Value used in the model |
|---|---|---|---|
| PALOMA-2 median PFS, palbociclib + letrozole (investigator; extended follow-up, median about 38 months) | 27.6 months | Rugo 2019 (Breast Cancer Res Treat 174:719-729, abstract; doi 10.1007/s10549-018-05125-4) | as published |
| PALOMA-2 median PFS, placebo + letrozole | 14.5 (12.9 to 17.1) months | Finn 2016 p8; Rugo 2019 | as published |
| PALOMA-2 HR (extended follow-up) | 0.563 | Rugo 2019 | as published |
| PALOMA-2 objective response, measurable disease | 55.3% (n = 338) / 44.4% (n = 171) | Finn 2016 pp. 8-9, Table 3 | as published |
| BOLERO-2 median PFS, everolimus + exemestane (final analysis, local review) | 7.8 (6.9 to 8.5) months | Yardley 2013 pp. 2, 5; everolimus prescribing information p43 | as published |
| BOLERO-2 median PFS, exemestane | 3.2 (2.8 to 4.1) months | same | as published |
| BOLERO-2 HR | 0.45 (0.38 to 0.54) | same | as published |
| evERA control median PFS, overall / ESR1-mutated | 5.5 (4.0 to 5.6) / 5.5 (3.8 to 5.6) months | Mayer 2026 p7 | as published |
| acelERA median PFS, giredestrant / physician's choice | 5.6 / 5.4 months | Martín 2024 (acelERA) p6 | as published |
| acelERA HR, overall / ESR1-mutated | 0.81 (0.60 to 1.10) / 0.60 (0.35 to 1.03) | Martín 2024 p6 | as published |
| MONALEESA-2 median overall survival, ribociclib + letrozole | 63.9 (52.4 to 71.0) months | Hortobagyi 2022 pp. 1, 4 | as published |
| evERA control 18-month overall survival (survival calibration) | 62.0% (54.3 to 69.7) | Mayer 2026 p8 | as published |
| persevERA, PALOMA-3 and PADA-1 step 1 | see the sections on these trials | - | calibration targets |

## Toxicity calibration targets

| Target | Published value | Source | Value used in the model |
|---|---|---|---|
| Grade 3-4 neutropenia, palbociclib (PALOMA-2) | 66.4% (grade 3 56.1% + grade 4 10.4%) | Finn 2016 p1, Table 2 | 66.4% |
| Grade 3-4 neutropenia, ribociclib (MONALEESA-2) | 59.3% (Hortobagyi 2016); 63.8% (Hortobagyi 2022) | Hortobagyi 2016 p1, Table 3 | 59.3% |
| Grade 3-4 neutropenia, abemaciclib (MONARCH 3) | 21.1% (primary analysis); 23.9% (final analysis); 22.1% (laboratory-based, primary analysis) | Goetz 2017 p1, Table 3; Johnston 2019 (MONARCH 3 final analysis, npj Breast Cancer 5:5) p1 | 21.1% (primary analysis, adverse-event based like the other two) |
| Palbociclib stopped for adverse events | 7.4% palbociclib; 9.7% any study drug | Finn 2016 p7 | 7.4% |
| CDK4/6 inhibitor median treatment duration | palbociclib 19.8 months; ribociclib 13.0 months; abemaciclib 15.1 months | palbociclib prescribing information p8; Hortobagyi 2016 p4; abemaciclib prescribing information p13 | drug-specific values |
| Ribociclib grade 3-4 ALT increase | 9.3% (7.5% + 1.8%) | Hortobagyi 2016 pp. 6, 9 | as published |
| Ribociclib QTcF above 480 ms | 3.3% (11 patients) | Hortobagyi 2016 pp. 6, 9 | as published |
| Abemaciclib grade 3 diarrhoea | 9.5% (no grade 4) | Goetz 2017 pp. 1, 4 | as published |
| Abemaciclib grade 2 or higher diarrhoea | 36.7% (27.2% + 9.5%, primary analysis); 39.8% (final analysis) | Goetz 2017 p4; Johnston 2019 p5 | 36.7% |
| evERA control: stomatitis any grade / grade 3 | 48.9% / 3.2% | Mayer 2026 Table 2, p11 | as published |
| evERA control: pneumonitis any grade | 12.4% | same | as published |
| evERA control: grade 3-4 hyperglycaemia, grade 3 diarrhoea, grade 3-4 anaemia | 2.2%, 1.1%, 8.6% | same | as published |
| evERA control: everolimus / endocrine therapy stopped for adverse events | 11.8% / 6.5% | Mayer 2026 p10 | as published |
| evERA experimental arm: bradycardia; giredestrant stopped for adverse events | 3.8%; 8.2% | Mayer 2026 p10 | as published |
| evERA treatment durations | standard endocrine therapy 4.6 months (everolimus 4.1); giredestrant 7.1 months (everolimus 6.0) | Mayer 2026 p7; appendix Table S4 | 140 d and 216 d (endocrine partner; everolimus can be stopped earlier in the model) |

## Held-out validation targets

| Item | Published value | Source | Value used in the model |
|---|---|---|---|
| evERA median PFS, giredestrant + everolimus: ESR1-mutated / overall / no ESR1 mutation | 10.0 (8.1 to 12.9) / 8.8 (6.6 to 9.6) / 5.7 months | Mayer 2026 pp. 7, 10 | as published |
| evERA HR: ESR1-mutated / overall / no ESR1 mutation | 0.38 (0.27 to 0.54) / 0.56 (0.44 to 0.71) / 0.84 (0.59 to 1.18) | Mayer 2026 pp. 7, 10 | as published |
| evERA 6- and 12-month PFS, giredestrant + everolimus vs control | ESR1-mutated 66.1% vs 38.1% and 40.5% vs 15.2%; overall 57.4% vs 39.6% and 34.1% vs 18.1% | Mayer 2026 Fig 2 | as published |
| evERA objective response, giredestrant + everolimus vs control | ESR1-mutated 26.6% vs 13.8%; overall 23.8% vs 11.7% | Mayer 2026 pp. 8-9 | as published |
| evERA 18-month overall survival, giredestrant + everolimus vs control | overall 75.1% vs 62.0%; ESR1-mutated 71.5% vs 50.9% | Mayer 2026 pp. 8-9 | as published |
| evERA ESR1-mutated share | 55.5% (207/373) | Mayer 2026 p5 | as published |
| evERA experimental-arm toxicity | stomatitis 47.3% any grade, 2.7% grade 3; pneumonitis 17.0% any grade, 2.7% grade 3 or higher; everolimus stopped for adverse events 17.0% | Mayer 2026 Table 2 | as published |
| BOLERO-4 median PFS, everolimus + letrozole | 22.0 (18.1 to 25.1) months | Royce 2018 pp. 1, 4 | as published |
| postMONARCH HR | 0.73 (0.57 to 0.95) | Kalinsky 2025 pp. 1, 4 | as published |
| postMONARCH median PFS | 6.0 / 5.3 months | Kalinsky 2025 | as published |
| MAINTAIN HR | 0.57 (0.39 to 0.85) | Kalinsky 2023 pp. 1, 5 | as published |
| TRINITI-1 median PFS | 5.7 (3.6 to 9.1) months | Bardia 2021 p4 | as published |
| SERENA-6, PADA-1 step 2, MIRACLE, PARSIFAL and FALCON | see the sections on these trials | - | validation |

## Tumour assessment schedules

| Trial | Schedule | Source |
|---|---|---|
| PALOMA-2 | every 12 weeks (11 to 13) | Finn 2016 p3 |
| MONALEESA-2 | every 8 weeks for 18 months, then every 12 weeks | Hortobagyi 2016 p3 |
| BOLERO-2 | every 6 weeks | Baselga 2012 pp. 4-5 |
| BOLERO-4 | every 8 weeks | Royce 2018 p2 |
| evERA, acelERA, SERENA-6 | every 8 weeks for 18 months, then every 12 weeks | Mayer 2026 p3; Martín 2024 p3; Bidard 2025 p3 |
| postMONARCH | every 8 weeks for 12 months, then every 12 weeks | Kalinsky 2025 p3 |
| MAINTAIN | every 12 weeks | Kalinsky 2023 pp. 2-3 |
| TRINITI-1 | every 8 weeks for 12 months, then every 12 weeks | Bardia 2021 p2 |
| persevERA, PARSIFAL, FALCON | not in the abstracts; model default (every 8 weeks for 18 months, then every 12 weeks) | - |
| In silico strategy trial (arms A to D) | every 8 weeks for 18 months, then every 12 weeks | current trial practice (evERA, SERENA-6) |

## Pharmacology

| Item | Published value | Source | Value used in the model |
|---|---|---|---|
| Palbociclib half-life | 29 (+/-5) h | palbociclib prescribing information p24 | as published |
| Ribociclib half-life | effective 32.0 h | ribociclib prescribing information p29 | as published |
| Abemaciclib half-life | 18.3 h | abemaciclib prescribing information p26 | as published |
| Giredestrant half-life | 43.0 h at 30 mg (range 25.8 to 43.0 h) | Malhi 2023 pp. 3-4 | as published (approximate) |
| Everolimus half-life | about 30 h | everolimus prescribing information | as published |
| Palbociclib: midazolam AUC | +61% | palbociclib prescribing information p25 | 1.61 |
| Ribociclib: midazolam AUC | 3.8-fold at 400 mg (measured); 5.2-fold at 600 mg (physiologically based prediction) | ribociclib prescribing information p30 (section 12.3) | 5.2 at 600 mg; the model then gives 4.8 at 400 mg and 4.4 at 300 mg |
| Ribociclib + everolimus 2.5 mg (TRINITI-1) | everolimus trough 2-3-fold (ribociclib 250 mg) and 4-5-fold (300 mg) above expected | Bardia 2021 p4 | consistent with the value used |
| Abemaciclib CYP3A4 inhibition | no clinically meaningful effect on CYP3A4 substrates | abemaciclib prescribing information p28 | as published |
| Everolimus with ketoconazole / erythromycin / verapamil | AUC 15-fold / 4.4-fold / 3.5-fold | everolimus prescribing information p41 | fm = 0.933 |
| Giredestrant interaction with everolimus | no data; low CYP3A induction potential; no interaction with palbociclib | Malhi 2023 pp. 1, 5 | not reported; varied in the uncertainty analysis |
| Dose levels and modification rules | palbociclib 125/100/75, ribociclib 600/400/200, abemaciclib 150/100/50 mg; everolimus 10/5/2.5 mg; rules as in the prescribing information | prescribing information | as published |

## Real-world attrition

| Source | Population | Not endocrine-based second line | Use in the model |
|---|---|---|---|
| Varnier 2025 (ESME-MBC, France) | 538 patients who received a second line after a first-line CDK4/6 inhibitor | first-line PFS below 6 months: 77.9%; 6 to 12 months: 57.3% (75/131); above 12 months: 29.2% (75/257); 6 months or more pooled 38.7% | base case: 57% (6 to 12 months), 29% (above 12 months) |
| Mahtani 2026 (6 countries, physician survey) | 767 alive on the second line | first line below 12 months: 39% chemotherapy; 12 months or more: 25% | supports duration dependence |
| Raghavendra 2026 (MD Anderson) | 1290 with a first-line CDK4/6 inhibitor | 38% chemotherapy + 7% investigational (all durations) | consistent |
| Read 2024 (French national claims) | 1036 with a second line after a CDK4/6 inhibitor + aromatase inhibitor | 50.6% chemotherapy-based (all durations, short follow-up) | consistent |
| MONALEESA-2 (trial) | 267 treated after ribociclib | 31.8% chemotherapy-containing first subsequent therapy | consistent |

## persevERA: first-line giredestrant + palbociclib vs letrozole + palbociclib (abstract, ASCO 2026, LBA1006)

| Item | Published value | Source |
|---|---|---|
| Population | n = 992 (495 / 497); no prior systemic therapy for advanced disease; recurrent with treatment-free interval above 12 months 69.3%; de novo 19.1%; visceral 60.4%; median age 63 | Turner 2026 (abstract LBA1006) p1 |
| Investigator PFS, HR | 0.89 (0.76 to 1.05), p = 0.155; 623 events; median follow-up 52.2 months | Turner 2026 p1 |
| Investigator PFS, medians | 33.1 (30.2 to 38.3) vs 28.2 (25.0 to 33.1) months | Turner 2026 p1 |
| Objective response rate | 60.2% vs 58.8% | Turner 2026 p1 |
| Overall survival | HR 1.03 (0.83 to 1.28); medians not reached | Turner 2026 p1 |
| Duration of response | 38.5 vs 30.4 months, HR 0.80 (0.63 to 1.01) | Turner 2026 p1 |
| Endocrine therapy stopped for side effects | 6.5% (giredestrant) vs 5.5% (letrozole) | Turner 2026 p1 |

## persevERA: data of the conference presentation

| Item | Giredestrant + palbociclib | Letrozole + palbociclib | Use in the model |
|---|---|---|---|
| PFS at 12 / 24 / 36 months | 77.8 / 59.7 / 45.8% | 77.3 / 57.3 / 41.9% | calibration targets |
| Overall survival at 24 / 36 months | 84.9 / 74.1% | 84.2 / 74.1% | descriptive |
| Best response: complete / partial / stable / progressive disease | 4.1 / 56.1 / 31.4 / 4.5% (n = 465) | 4.5 / 54.4 / 30.9 / 6.2% (n = 469) | not used |
| Median duration of treatment | giredestrant 29.8; palbociclib 28.3 months | letrozole 25.2; palbociclib 24.2 months | exposure in the toxicity-only emulations of discontinuation (giredestrant 29.8, letrozole 25.2 months) |
| Mean dose intensity, palbociclib | 81.3% | 81.8% | not used |
| Bradycardia | 11.7% | 1.8% | not used |
| Design | HR 0.77 targeted with 89% power, about 620 PFS events; efficacy boundary HR 0.85 |  | - |
| Source | Seymour 2026 (OncLive) and Targeted Oncology Staff 2026 (Targeted Oncology), reports of the presentation with identical values; not in the abstract |  |  |

## SONIA: CDK4/6 inhibitor in first line vs second line (Sonke, Nature 2024)

| Item | Published value | Source |
|---|---|---|
| Design | 1,050 patients: 524 CDK4/6 inhibitor first (aromatase inhibitor + CDK4/6 inhibitor, then fulvestrant), 526 CDK4/6 inhibitor second (aromatase inhibitor, then fulvestrant + CDK4/6 inhibitor); imaging every 12 weeks; intended CDK4/6 inhibitor palbociclib 91%, ribociclib 8%, abemaciclib 1% | Sonke 2024 pp. 2, 3 (Table 1), 8 |
| Population | median age 64; visceral 56%; bone-only 17%; de novo 35%; previous (neo)adjuvant endocrine therapy 49% | Sonke 2024 p2 |
| PFS1 | 24.7 (22.1 to 28.3) vs 16.1 (14.4 to 17.1) months; HR 0.59 (0.51 to 0.69) | Sonke 2024 p3 |
| PFS2 (primary end point) | 31.0 (26.4 to 34.5) vs 26.8 (24.1 to 30.4) months; HR 0.87 (0.74 to 1.03), P = 0.10; median follow-up 37.3 months | Sonke 2024 p2 |
| Overall survival | 45.9 vs 53.7 months; HR 0.98 (0.80 to 1.20); 372 deaths | Sonke 2024 p2 |
| Started the protocol second line after stopping line 1 | 246 / 315 (78%) vs 349 / 404 (86%); chemotherapy instead 31 vs 21; off-protocol endocrine therapy 9 vs 4; no second-line therapy at cut-off 29 vs 30 | Sonke 2024 p2 |
| CDK4/6 inhibitor treatment duration | median 24.6 vs 8.1 months; restricted mean 28.5 vs 12.4 months | Sonke 2024 p2 |
| Safety | grade 3 or higher events 2,763 vs 1,591; patients with at least one grade 3 or higher event 83% vs 64%; CDK4/6 inhibitor dose reduction 43% vs 39%; CDK4/6 inhibitor stopped while continuing endocrine therapy 9% vs 6% | Sonke 2024 p3 |

## SONIA circulating tumour DNA analysis (Jongbloed, Nat Med 2025)

| Item | Published value | Source |
|---|---|---|
| Cohort | 409 patients with baseline cell-free DNA (220 first-line / 189 second-line strategy); high aneuploidy score (5 or more, about 10% tumour DNA or more) in 141 (34%) | Jongbloed 2025 p3 |
| High score | PFS1 21.1 vs 8.9 months; PFS2 24.7 vs 18.0 months; PFS2 HR 0.58 (0.38 to 0.88) | Jongbloed 2025 pp. 3, 4 |
| Low score | PFS1 23.7 vs 23.1 months; PFS2 32.5 vs 39.0 months; PFS2 HR 1.36 (0.95 to 1.96) | Jongbloed 2025 pp. 3, 4 |
| Strategy x score interaction (PFS2) | P = 0.004 | Jongbloed 2025 p3 |

## First-line fulvestrant (validation; values from the abstracts)

| Trial | Comparison | Published value | Source |
|---|---|---|---|
| PARSIFAL (phase 2, endocrine-sensitive, n = 486) | fulvestrant + palbociclib vs letrozole + palbociclib | median PFS 27.9 vs 32.8 months; HR 1.13 (0.89 to 1.45); response 46.5% vs 50.2% | Llombart-Cussac 2021, JAMA Oncol 7:1791-1799, doi 10.1001/jamaoncol.2021.4301 |
| FALCON (phase 3, endocrine-therapy-naive, n = 462) | fulvestrant 500 mg vs anastrozole | median PFS 16.6 vs 13.8 months; HR 0.797 (0.637 to 0.999) | Robertson 2016, Lancet 388:2997-3005, doi 10.1016/S0140-6736(16)32389-3 |

## PALOMA-3: fulvestrant with palbociclib or placebo after endocrine therapy (Cristofanilli, Lancet Oncol 2016)

| Item | Published value | Source (page) | Use in the model |
|---|---|---|---|
| Design | double-blind phase 3, 2:1 (347 / 174); progression on or within 1 month of endocrine therapy for advanced disease, or on or within 12 months of adjuvant therapy; one prior chemotherapy line allowed; pre/perimenopausal women with goserelin | pp. 1, 3 | emulation: line 1 letrozole alone, then fulvestrant +/- palbociclib for every patient who progresses |
| Population | postmenopausal 79%; visceral 59-60%; documented sensitivity to prior endocrine therapy 79%; prior aromatase inhibitor 85-87%; at least one prior line for advanced disease 79%; chemotherapy for advanced disease 33-37% | pp. 5-6 (Table 1), p7 (Fig 2B-C) | - |
| Tumour assessment | every 8 weeks (+/- 7 days) for the first year, then every 12 weeks | p3 | emulation schedule |
| Investigator PFS (final analysis, median follow-up 8.9 months) | 9.5 (9.2 to 11.0) vs 4.6 (3.5 to 5.6) months; HR 0.46 (0.36 to 0.59); 145 / 114 events | pp. 1, 9 | calibration targets |
| Central review (random 40% sample) | HR 0.37 (0.23 to 0.59) | supplement p19 | not used |
| Endocrine-sensitivity subgroups | sensitive HR 0.42 (0.32 to 0.56), median 10.2 vs 4.2; not sensitive HR 0.64 (0.39 to 1.07), 7.4 vs 5.4; interaction P = 0.13 | p9 (Fig 3) | descriptive |
| Objective response, measurable disease | 24.6% (66 / 268) vs 10.9% (15 / 138) | p10, p8 (Table 2) | calibration targets |
| Discontinuation for adverse events (any study treatment) | 14 (4%) vs 3 (2%) | p11 | descriptive |
| Grade 3-4 neutropenia | 223 / 345 (65%) vs 1 (1%) | p1, p12 (Table 3) | descriptive |
| ESR1, overall survival, treatment duration | not reported in this paper | - | - |

## MIRACLE: first-line letrozole with or without everolimus after tamoxifen (Fan, JAMA Oncol 2021)

| Item | Published value | Source | Use in the model |
|---|---|---|---|
| Design | open-label phase 2, China, 199 premenopausal patients (101 / 98) with recurrence within 24 months of adjuvant tamoxifen or progression on tamoxifen for advanced disease; letrozole 2.5 mg + goserelin 3.6 mg every 28 days, with or without everolimus 10 mg; crossover allowed | pp. 1-3 | emulation: first-line letrozole +/- everolimus, no mouthwash prophylaxis |
| Tumour assessment | every 8 weeks, then every 12 weeks after the first year if stable or responding | p3 | emulation schedule |
| Investigator PFS | 19.4 (16.3 to 22.0) vs 12.9 (7.6 to 15.7) months; HR 0.64 (0.46 to 0.89); P = .008 | pp. 1, 4, 5 (Fig 2A) | validation (HR within the 95% CI; medians descriptive) |
| PFS, supplement eFigure 2 | 16.5 vs 10.0 months, HR 0.638 (0.458 to 0.890); 67 / 75 events | Supplement 2 p3 | medians differ from the main text with the same HR, so the HR is the criterion |
| Endocrine-resistance subgroups | secondary resistance HR 0.50 (0.32 to 0.79); primary resistance HR 0.78 (0.45 to 1.37) | p4; Supplement 2 p3 | descriptive |
| Objective response, measurable disease | 50.0% (33 / 66) vs 39.3% (24 / 61) | pp. 4-5 (Table 2) | validation (descriptive) |
| Adverse events leading to discontinuation | 9 (8.9%) vs 0 / 92 | p5 | validation (within 5 points) |
| Stomatitis, any grade / grade 3 or higher | 79.2% / 9.9% | Supplement 2 p4 (eTable) | validation (no-mouthwash scenario; within 10 points) |
| Treatment duration, follow-up | not reported | - | - |

## PADA-1 (Bidard, Lancet Oncol 2022) and its ESR1 kinetics analysis (Cabel, Ann Oncol 2026)

| Item | Published value | Source | Use in the model |
|---|---|---|---|
| Design, step 1 | 1017 patients, first-line aromatase inhibitor + palbociclib; ctDNA (multiplex ddPCR, ESR1 codons 380, 536, 537, 538) at inclusion, after one cycle, then every two cycles; imaging after three cycles, then every two cycles | Bidard pp. 1, 3; Cabel p7 | emulation schedule |
| Rising mutation (definition) | new detection after a negative baseline, or a significant increase of a baseline mutation; a baseline mutation alone is not rising | Bidard supplement Table 1; Cabel p3 | baseline-positive patients are not counted as rising (the 3-fold rule is not emulated) |
| Detection limit | not reported; lowest allele fraction at detection 0.1%, median 0.8% (range 0.1-25.8); not re-detected at randomisation (about one month later) in 75/161 (46.6%) | Cabel p8 | effective limit calibrated (tum.LOD_pd); allele fraction at detection reported, not fitted |
| Baseline ESR1 mutation | 33/1017 = 3.2% (2.2 to 4.5) | Cabel p3 | calibration target |
| Rising mutation during step 1 | 279/1017 (27%) at the primary lock (median follow-up 34.5 months); update 283; 283 of the 678 patients with progression or a rising mutation (41.7%) | Bidard pp. 1, 6; Cabel pp. 3, 7 | calibration target: 41.7% (monitoring stopped 30-52 months after inclusion) |
| Synchronous progression (before detection or within 30 days after) | 71/283 (25%); primary analysis 69/279 | Cabel p3 (Fig 1); Bidard p4 (Fig 1) | calibration target: 25% |
| Time on aromatase inhibitor + palbociclib before randomisation | median 15.6 months (IQR 10.1 to 22.7) | Bidard p7 | calibration target (model: detection + 30 days) |
| Step-1 PFS | median 26.8 months (24.8 to 28.6); update 27.4 | Bidard p7; Cabel p3 | calibration target: 26.8 |
| Net cumulative incidence of rising mutation (Fig 2D) | about 12.9 / 28.2 / 38.6% at 12 / 24 / 36 months (read from the curve) | Cabel p5 | descriptive |
| Step 2 (randomised, n = 172) | fulvestrant + palbociclib 11.9 (9.1 to 13.6) vs aromatase inhibitor + palbociclib 5.7 (3.9 to 7.5) months; stratified HR 0.61 (0.43 to 0.86); update 12.8 vs 5.9 months, HR 0.54 (0.38 to 0.75); PFS2 29.4 vs 14 months, HR 0.35 | Bidard pp. 1, 7; Cabel pp. 7, 9 | validation (held out) |
| Step 2 eligibility | rising mutation without RECIST progression on imaging within 4 weeks before randomisation; detection sample about one month before randomisation (Fig S6, read from the figure) | Bidard p3; Cabel supplement | emulation: randomisation 30 days after detection, imaging check at detection and at randomisation |

## SERENA-6 extended analysis (Turner 2026)

| Item | Published value | Source | Use in the model |
|---|---|---|---|
| Design | first-line aromatase inhibitor + CDK4/6 inhibitor for 6 months or more; ctDNA (Guardant360 CDx) every 2-3 months; ESR1 mutation without progression (imaging within 28 days before randomisation): camizestrant + the same CDK4/6 inhibitor vs continuing the aromatase inhibitor + CDK4/6 inhibitor, double-blind; imaging every 8 weeks for 18 months, then every 12 weeks | Turner 2026 pp. 1, 3-4; protocol pp. 31, 38 | emulation (giredestrant as the SERD proxy) |
| Detection limit | not reported ("highly sensitive assay"); median allele fraction at a subsequent-test detection about 0.26%, first-test 0.63% (read from the figure) | Turner 2026 p6 (Fig 1C) | effective limit = PADA-1 limit x 0.26 / 0.8 |
| Detection to randomisation | not reported (results within 21 days, screening within 28 days) | protocol pp. 132-133 | assumption: 35 days |
| Screening | 3256 tested; 548 ESR1 mutations (16.8%); concurrent progression 117/548 (21%); 315 randomised (157 / 158); start of first-line therapy to detection, median 22 months | Turner 2026 pp. 1, 7-8 | descriptive |
| CDK4/6 inhibitor | palbociclib 75.6%, ribociclib 14.9%, abemaciclib 9.5% | Turner 2026 p7 (Table 1) | model: palbociclib |
| PFS (extended analysis) | 16.8 (14.7 to 19.4) vs 9.2 (7.2 to 9.7) months; HR 0.45 (0.34 to 0.59); 24-month PFS 34.9% vs 14.2% | Turner 2026 pp. 8-9 | validation |
| PFS, interim analysis (Bidard 2025) | 16.0 vs 9.2 months; HR 0.44 (0.31 to 0.60) | Turner 2026 p3 | not used (the extended analysis is the validation target) |
| PFS2 | 25.7 vs 19.1 months; HR 0.63 (0.46 to 0.86) | Turner 2026 p10 | not emulated (subsequent therapy not modelled) |
| Overall survival (immature) | 41.2 vs 40.2 months; HR 0.87 (0.57 to 1.30); 95 deaths, 58% information fraction | Turner 2026 p11 | validation (HR within the 95% CI) |
| Adverse events | treatment stopped for adverse events: camizestrant / aromatase inhibitor 1.3% vs 3.2%; CDK4/6 inhibitor 2 vs 4 patients; photopsia 21.3% vs 7.7%; bradycardia 10.3% vs 0% | Turner 2026 pp. 11-12 | not used |
