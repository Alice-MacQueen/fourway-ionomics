# PCE; new phytologist

Ionomic analysis was run on _380_ plants at three locations for 14 elemental compositions. 

Forage analysis has also been conducted by collaborators at 10 locations, producing estimates of cellulose content, lignin content, nitrogen content, and other variables that still need to be defined. We were interested in determining if we could use forage analysis values to predict ionomic quantities at these 10 sites, after training a model on both ionomic and forage analysis data at the three sites with ionomic data. Given the relative cost of each method (ionomic data costs $5-$6 per sample, and forage analysis is much cheaper), using forage data as a proxy for ionomic data could enable (larger) follow-up studies.

Ca44, Se78, Sr88, and Cd111 had bimodal distributions of abundance within plant samples across all three sites. Many had to have outlier samples removed `(did these fall into any sort of pattern of site or individual?)`

For each ion, heritability and V_a (additive genetic variance) varied extensively between the three sites. Na23, Mg26, Al27, P31, K39, Ca44, Mn55, Cu63, Rb85, Sr88, and Cd111 had heritabilities of 0.2 or greater for at least one site, while B11, Fe54, Co59, Zn66, As75, Se78, and Mo98 did not. PKLE had heritabilities higher than 0.2 for 7 elements (none unique to PKLE); CLMB for 9 elements (Al27 uniquely high at CLMB), and KBSM for 10 elements (Cd111 and K39 uniquely high heritability here). Mg26, P31, Ca44, Mn55, and Rb85 had heritabilities that were not significantly different between all sites; Mg26, Mn55, and Rb85 in particular were the most similar.

In terms of additive genetic variance (V_a) and error variance (V_e), PKLE had the most variance for Na23, Mg26, Al27, K39, Ca44, Zn66, As75, Rb85, Sr88, and Mo98, 10 elements, and 7 with significant V_a, but only 3 with H^2 > 0.2 at PKLE - thus, the variance in the ion quantities at this site were likely due to soil characteristics in most cases, and not under explicit genetic control. CLMB had the most variance for P31, Fe54, Mn55, Co59, Cu63, and Se78, 6 elements, and 3 with significant V_a, and 3 with H^2 > 0.2. KBSM had the most variance for B11 and Cd111, 2 elements, and 1 with significant V_a, Cd111, which also had the highest heritability (0.5) at this site - thus, the variance in Cd111 quantities at this site were affected by soil characteristics but were also under genetic control. 

## Elemental abundance correlations (phenotypic)

Of 18 total ions, four elements were not analyzed further due to machine errors in detecting elemental concentrations. 2 didn't run due to negative abundances, gave too few values to converge. 2 of the 16 remaining traits had no QTL - As75 & B11 did not have QTL. Co59 & Se78 didn't work due to negative values. 

At CLMB, the majority of elements were not strongly correlated with one another. Ca44 and Sr88 were highly correlated (0.77), as were Co59 and Cu63 (0.57), Fe54 and Al27 (0.51), Mn55 and Sr88 (0.5), Ca and Zn (0.43), Co59 and Fe54 (0.43), and Cu63 and Fe54 (0.44), Al27 and Cu63 (0.44). I believe these are all 2+ elements.

At PKLE, B11 and Sr88 were highly correlated (0.74), Al27 and B11 (0.65), Al27 and Sr88 (0.6), Ca44 with Sr88 and Zn66 (0.55, 0.49), and P31 with Rb85 (0.42).

Phenotypic correlations were typically larger at KBSM than at PKLE or CLMB.
As75 with Se78 (0.64), B11 with Mo98 & Na23 (0.69, 0.58), Ca44 with Sr88 and Zn66 (0.69 and 0.63), Cu63 with P31 (0.43), Fe54 with P31 (0.43), K39 with Rb85 (0.79), Mg26 with Sr88 (0.49), Mn55 with Sr88 (0.45), and P31 with Rb85 (0.41). 

Within sites, all strong (>0.3) correlations were positive. However, when correlations were made across all sites, there were three strong negative correlations: Ca44 and Cd111 (-0.49), Ca44 and Mn55 (-0.47), and Cd111 and Na23 (-0.38). 

`It would be nice to have correlation plots that are rearranged, for ease of interpretation...`

`For the QTL summary plot, the colors are too similar and the key goes off the page - I can't interpret the overlaps. Also, there are some small black dots which don't seem to correlate with any QTL.`

QTL regions on most chromosomes showed QxE. 16 of 18 chromosomes had QTL, and 14 of 18 had QTL with QxE.

## LOD scores

QTL analysis of the 14 remaining elemental compositions found between two and 14 QTL regions per element that had a LOD threshold above `3.5`. P31 had 14 QTL regions, while Cd111, Fe54, Mo98, and Na23 had 2.

B11, Fe54, Co59, Zn66, As75, Se78, and Mo98 had low heritabilities (<0.2) at all three sites. As75 & B11 did not have QTL. Co59 & Se78 were dropped from the analysis due to negative values for concentrations. Fe54, Zn66, and Mo98 had 2, 3, and 2 QTL, respectively.

The largest LOD score was for Cd111 on Chr02K@32.06 `(Mb? what?)`. In the C x D cross, this region had large effects on Cd111 levels at KBSM and CLMB. This QTL overlapped a QTL for Sr88, but this QTL fell into a much wider interval, and barring a mechanism affecting these elements the same way, I doubt these have pleiotropy.

## Overlapping QTL

`The QTL on Chr01N seems to be rendered incorrectly.`

Overall, there were 21 sets of QTL that had overlapping intervals, and 20 QTL that did not overlap another ionomics QTL. Rb85 Chr01K, Mg26 Chr01N, Mg26 on Chr02K, K39 on Chr02N, two P31 QTL on Chr03K, Mg26 on Chr03N, two P31 and a Cu63 on Chr04K, Mg26 on Chr04N, Rb85 on Chr05K, P31, Na23, and Mg26 on Chr05N, Mg26 on Chr07N, Sr88 on Chr08K, Rb85 on Chr08N, Sr88 on Chr09K, and K39 on Chr09N. 

77 QTL total
5 P31 single, 9 overlap
6 Mg26 single, 3 overlap
2 Sr88 single, 7 overlap
3 Rb85 single, 5 overlap
2 K39 single, 5 overlap
0 Ca44 single, 5 overlap
1 Cu63 single, 4 overlap
0 Mn55 single, 5 overlap
0 Al27 single, 4 overlap
0 Zn66 single, 3 overlap
0 Cd111 single, 2 overlap
0 Fe54 single, 2 overlap
0 Mo98 single, 2 overlap
1 Na23 single, 1 overlap

More Mg26 non-overlapping QTL than expected (6, and 2.3 expected). I bet some of these could be related to Mg's essential function in chlorophyll molecules. Has more QTL 
Fewer Ca44 non-overlapping QTL than expected (0, and 1.3 expected), and fewer Al27 non-overlapping QTL than expected (0, 1 expected). I bet Ca might overlap a lot with other 2A elements, or 2+ elements, like Mg, Sr, Zn, even Fe or Cu. Not sure about Al...

On Chr01K, there were two sets of overlapping QTL - P31, Cu63, and K39 near 14.42, and Cu63 and P31 near 71-75.
On Chr02K, there was one set of overlapping QTL - Cd111 and Sr88 @ ~32.
Chr02N had three sets of overlapping QTL - Mg26 and P31 from 0-12, Ca44, Sr88, and Zn66 from 20-33, Rb85, P31, and Cd111 from 76-79.
Chr03K had one set of overlapping QTL
Chr03N had one set of overlapping QTL
Chr04N had one set of overlapping QTL
Chr05K had two sets of overlapping QTL
Chr07K had two sets of overlapping QTL
Chr07N had one set of overlapping QTL
Chr08K had one set of overlapping QTL
Chr09K had three sets of overlapping QTL
Chr09N had three sets of overlapping QTL


## Ion functions and/or toxicity

What are major functions for each of these elements in plant cells?
Al Aluminum toxicity is a major constraint for crop production in acidic soil worldwide (Panda, Baluska, Matsumoto 2009). The target of Al toxicity is the root tip, in which Al exposure causes inhibition of cell elongation and cell division, leading to root stunting accompanied by reduced water and nutrient uptake. Al3+ is the most toxic form. Approximately 50% of arable land is negatively impacted by Al toxicity due to acidic soil.

Ca Calcium uptake by the plant is passive and does not require energy input. Calcium uptake is directly related to the plant transpiration rate. Conditions of high humidity, cold, and low transpiration rates may result in calcium deficiency. Calcium promotes proper plant cell elongation, strengthens cell wall structure by forming calcium pectate compounds, helps protect plants against heat stress by improving stomata function and participating in induction of HSPs. Usually soils with higher pH levels contain more available calcium. Calcium competes with other positively charged ions, such as Na+, K+, and Mg2+. 

Cd Cadmium is a environmental pollutant with toxic effects on plants. It can displace essential metals from a wealth of metalloproteins and disturb normal physiological processes and cause severe developmental aberrance. ROS overproduction, chloroplast structure change, cytoskeleton effects, vesicular trafficking effects, cell wall formation effects. Plants have cadmium pumps and transport cadmium into leaf vacuoles. (Wan and Zhang 2012)

Cu In plants, copper (Cu) acts as essential cofactor of numerous proteins. Cu is an essential player in electron transport, we also review the recent insights into the molecular mechanisms controlling chloroplastic and mitochondrial Cu transport and homeostasis. We finally highlight the involvement of numerous Cu-proteins and Cu-dependent activities in the properties of one of the major Cu-accumulation sites in plants: the cell wall. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4859090/ The rise of photosynthetic organisms on Earth has driven a progressive accumulation of oxygen in the environment. This oxidative atmosphere led to a decreased solubility of iron by the formation of iron oxides and to the progressive liberation of soluble Cu(II) from insoluble Cu sulfide salts (Burkhead et al., 2009). Since then, iron in biological molecules has been progressively substituted by Cu which is able to perform similar functions. This explains why many Cu-proteins have a functional counterpart that uses Fe as cofactor and why growth on a substrate with a toxic Cu level is commonly linked to a decreased Fe-content in roots and leaves (Pätsikkä et al., 2002; Burkhead et al., 2009; Festa and Thiele, 2011). Consequently, plant phenotypes associated with Cu toxicity share similarities with those related to Fe-deficiency, such as the presence of leaf chlorosis, decreased leaf chlorophyll content and enhanced oxidative stress (Pätsikkä et al., 2002). In contrast, copper deficient plants develop chlorotic symptoms that appear first at the tip of the young leaves prior forming necrotic lesions. Plants grown under Cu deficiency also show impairment in the photosynthetic transport chain and a reduction in non-photochemical quenching, which is consistent with a lack of plastocyanin (PC) function (Abdel-Ghany and Pilon, 2008).

Under physiological conditions, the transition metal Cu is found in the two common forms, the reduced Cu(I) state and the oxidized Cu(II) state. Depending on this state, Cu can bind different substrates. In living organisms, the main functions of Cu are the transport of electrons in mitochondria and chloroplasts (the most abundant Cu protein is plastocyanin, a photosynthesis-related protein involved in the transfer of electrons from cytochrome f to P700+), the control of the cellular redox state (a major Cu-binding protein is the Cu/Zn superoxide dismutase) but also the remodeling of the cell wall (Cohu and Pilon, 2010).

Fe Iron cofactors also function in oxygen transport, oxygen or iron sensing, or regulation of protein stability. The chloroplasts are particularly rich in iron–sulphur (FeS) proteins such as Photosystem I, ferredoxins and a range of metabolic enzymes. Mitochondria are another hotspot for iron enzymes, such as respiratory complexes containing multiple FeS clusters (complex I and II), a mix of FeS and haem (complex III) or haem and copper (complex IV). The peroxisomes and the endoplasmic reticulum contain haem proteins such as peroxidases and cytochrome P450s, whereas mono- and di-iron enzymes are found in all cell compartments. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5708359/

K Potassium is not only a constituent of the plant structure but it also has a regulatory function in several biochemical processes related to protein synthesis, carbohydrate metabolism, and enzyme activation. Several physiological processes depend on K, such as stomatal regulation and photosynthesis. In recent decades, K was found to provide abiotic stress tolerance. Under salt stress, K helps to maintain ion homeostasis and to regulate the osmotic balance. Under drought stress conditions, K regulates stomatal opening and helps plants adapt to water deficits. Many reports support the notion that K enhances antioxidant defense in plants and
therefore protects them from oxidative stress under various environmental adversities. Hasanuzzaman et al Agronomy 2018.

Mg Mg is particularly important to plants, with some 75% of leaf Mg involved in protein synthesis and 15–20% of total Mg associated with chlorophyll pigments [10], acting mainly as a cofactor of a series of enzymes involved in photosynthetic carbon fixation and metabolism [6], [11], [12]. In particular, Mg plays a central role in plant chlorophyll biosynthesis and carbon fixation as a cofactor of a series of enzymes involved in carbon metabolism.4 Conditions with dry soil and high levels of competing elements, such as potassium and calcium, also result in Mg deficiency.13 Mg is abundant in soil (the 8th most abundant element on earth), has high solubility in water for plant absorption, and has complex functionalities.3

Mn Manganese is used in plants as a major contributor to various biological systems including photosynthesis, respiration, and nitrogen assimilation. Manganese is also involved in pollen germination, pollen tube growth, root cell elongation and resistance to root pathogens. In chloroplasts, Mn is at the center of the oxygen-evolving complex of photosystem II (PSII) and is responsible for turning water into molecular oxygen and electrons to replenish those donated to the electron transport chain by photo-excited PSII. In mitochondria, Mn is the cofactor of superoxide dismutase that detoxifies superoxide (O2−). Like other micronutrients, Mn is absorbed by roots and distributed throughout the plant to sink organelles (plastids and mitochondria) and storage sites (vacuole). https://pdfs.semanticscholar.org/f455/1076c5ea42ece6c227caa47f14f7dad52db0.pdf Mn2+ is the only available metal form for plants. In this phase, Mn2+appears to be adsorbed by the negatively charged cell wall constituents of the root-cell apoplastic spaces (Humphries et al., 2007; Clarkson, 1988). Manganese also plays a role in ATP synthesis (Pfeffer et al., 1986), in RuBP carboxylase reactions (Houtz et al., 1988) and the biosynthesis of fatty acids, acyl lipids and proteins (Ness and Woolhouse, 1980). In addition, Mn plays a primary role in the activation and as cofactor of various enzymes in plants (~35) (Burnell, 1988), such as: Mn-superoxide dismutase, Mn-catalase, pyruvate carboxylase and phospho-enolpyruvate carboxykinase (Ducic and Polle, 2005). Manganese is also essential for the biosynthesis of chlorophyll (through the activation of specific enzymes), aromatic amino acids (tyrosine), secondary products, like lignin and flavonoids (Lidon et al., 2004). It also participates in the biosynthetic pathway of isoprenoids (Lidon et al., 2004) and assimilation of nitrate (Ducic and Polle, 2005). Hence, Mn is involved in metabolic processes such as respiration, photosynthesis, synthesis of aminoacids and hormone activation (indol acetic acid, IAA) throughout the IAA-oxidases (Burnell, 1988). 

Mo Molybdenum is an essential component in two enzymes that convert nitrate into nitrite (a toxic form of nitrogen) and then into ammonia before it is used to synthesize amino acids within the plant. It also needed by symbiotic nitrogen fixing bacteria in legumes to fix atmospheric nitrogen. Apart from Cu, Mo is the least abundant essential micronutrient found in most plant tissues and is often set as the base from which all other nutrients are compared and measured. Molybdenum is utilized by selected enzymes to carry out redox reactions. Enzymes that require molybdenum for activity include nitrate reductase, xanthine dehydrogenase, aldehyde oxidase and sulfite oxidase. (Kaiser et al 2005). Molybdenum is a transition element, which can exist in several oxidation states ranging from zero to VI, where VI is the most common form found in most agricultural soils.There are recent review articles on molybdoenzymes in plants, animals and prokaryotes (Mendel and Haensch, 2002; Williams and Frausto da Silva, 2002; Sauer and Frebort, 2003) that cover the extensive literature on the regulation and formation of Moco and the activity of Moco with molybdenum-dependent apoenzymes. Dissolved molybdenum available to plants is commonly found in the soluble
MoO−4 anion form (Lindsay, 1979). 

Na
P
Rb
Sr
Zn

## Go enrichments

Under physiological conditions, the transition metal Cu is found in the two common forms, the reduced Cu(I) state and the oxidized Cu(II) state. Depending on this state, Cu can bind different substrates. In living organisms, the main functions of Cu are the transport of electrons in mitochondria and chloroplasts (the most abundant Cu protein is plastocyanin, a photosynthesis-related protein involved in the transfer of electrons from cytochrome f to P700+), the control of the cellular redox state (a major Cu-binding protein is the Cu/Zn superoxide dismutase) but also the remodeling of the cell wall (Cohu and Pilon, 2010).

Cu63 QTL regions were significantly enriched for GO ontologies of cell wall biogenesis (p = 3.8x10-8), oxidoreductase activity (7.3x10-3), mitochondrion (7.7x10-6) & ADP binding (4.9x10-6), dna binding, regulation of transcription, and DNA-binding transcription factor activity. 

Fe54 QTL regions were significantly enriched for GO ontologies of lipid biosynthetic process (2.4x10-4), mitocondrion (1.3x10-3), oxidoreductase acitivty, and heme binding (3.5x10-3), amongst 15 ontologies.

K39 QTL regions were significantly enriched for GO ontologies of oxidoreductase activity, carbohydrate binding, and antioxidant activity.

Mg26 QTL regions were significantly enriched for GO ontologies of ADP binding and protein transport.

Mn55 QTL regions were significantly enriched for GO ontologies of multicellular organism development (1.7x-10), ubiquitin-dependent protein catabolic process (1.5x10-7), the mitochondrion, photosystem I, the photosystem I reaction center, and electron transfer activity (x10-3).

Mo98 QTL regions were significantly enriched for GO ontologies of 