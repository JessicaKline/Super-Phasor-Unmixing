# Super-Phasor-Unmixing
General Overview
The program phasor_calculations.m is a supplement to the publication “Chemical complexity of the retina
addressed by novel phasor analysis of unstained multimodal microscopy.” [CHEMPH_111091] In that
publication, we introduced the super-phasor unmixing (SPU) method to unmix seven spectroscopic signals
in the spectral range 300-690 nm. The performance of SPU was found to be significantly superior to fully
constrained linear spectral unmixing especially in regions of high spectral overlap. This program can serve
as a tool for implementing the super-phasor analysis. The reading of [CHEMPH_111091] is highly
encouraged before implementing phasor_calculations.m as that publication provides much greater detail
and context than these instructions. All equation numbers in these instructions reference
[CHEMPH_111091].

phasor_calculations.m has the ability to do two things. The first is to convert spectral coordinates into
phasor coordinates (eqs. 1 and 2) and perform a simple phasor analysis to find percent compositions (eq.
3). This function works with two and three pure components. The second function of phasor_calculations.m
is to perform super phasor analysis (eq 4). This function works with two and three super-phasors. Further
details on the program’s variables are addressed in Variable Specifications. A pseudo-run of the program
for a seven component analysis is provided in Implementation.

Variable Specifications

[phasor_coords, percentages] = phasor_calculations(pure_spectra, exp_spectrum, lambda_min,
lambda_max)

Input Variables

pure_spectra: An n x 4, n x 6, 2 x 2 or 2 x 3 matrix used to hold the flourescence spectra expected
to combine to form the mixed experimental spectrum.

An n x 4 matrix should be used when two pure spectra are involved and an n x 6 matrix
should be used for three pure spectra. These matrices should take the form of [λ ,
I , λ , I , λ , I] across the 4-6 columns with n rows.

The 2 x 2 and 3 x 2 matrices are used for super phasor analysis where each index represents
super-phasor G and S coordinates. These matrices should be of the form [G, S] and
have two columns with each row containg a set of super-phasor coordinates.

exp_spectrum: A m x 2 matrix used to hold the measured mixed flourescence spectrum. This matrix
should take the form [λ , I] with m rows.

lambda_min: the minimum wavelength for calculations

lambda_max: the maximum wavelength for calculations

Output Variables

phasor_coords: a 4 x 2 or 3 x 2 matrix that returns the various calculated phasor coordinates. The
columns in this matrix are of the format [G, S]. The first 2-3 rows are the phasor coordinates for
the pure spectra and are returned according to their input order in the variable pure_spectra. The
last row is always reserved for the set of phasor coordinates belonging to exp_spectrum.

percentages: a 3 x 1 or 2 x 1 matrix that returns the percent composition associated with each pure
component. The order is once again determined by the order in pure_spectra.

Implementation

𝑠𝑝𝑒𝑐1 = [𝜆𝑎 𝐼𝑎 𝜆𝑏 𝐼𝑏 𝜆𝑐
𝐼𝑐
] (nx6)

𝑠𝑝𝑒𝑐2 = [𝜆𝑑 𝐼𝑑 𝜆𝑒
𝐼𝑒] (nx4)

𝑠𝑝𝑒𝑐3 = [𝜆𝑓 𝐼𝑓 𝜆ℎ 𝐼ℎ] (nx4)


𝑒𝑥𝑝𝑆𝑝𝑒𝑐 = [𝜆𝑧
𝐼𝑧
] (nx2)

(𝑐𝑜𝑜𝑟𝑑𝑠1, 𝑝𝑒𝑟𝑐𝑒𝑛𝑡1) = 𝑝ℎ𝑎𝑠𝑜𝑟𝑐𝑎𝑙𝑐𝑢𝑙𝑎𝑡𝑖𝑜𝑛𝑠(𝑠𝑝𝑒𝑐1, 𝑒𝑥𝑝𝑆𝑝𝑒𝑐, 300, 500)

𝑐𝑜𝑜𝑟𝑑𝑠1 = [
𝐺𝑎 𝑆𝑎
𝐺𝑏 𝑆𝑏
𝐺𝑐
𝐺𝑧1
𝑆𝑐
𝑆𝑧1
]

𝑝𝑒𝑟𝑐𝑒𝑛𝑡1 = [
𝑎
𝑏
𝑐
]

(𝑐𝑜𝑜𝑟𝑑𝑠2, 𝑝𝑒𝑟𝑐𝑒𝑛𝑡2) = 𝑝ℎ𝑎𝑠𝑜𝑟𝑐𝑎𝑙𝑐𝑢𝑙𝑎𝑡𝑖𝑜𝑛𝑠(𝑠𝑝𝑒𝑐2, 𝑒𝑥𝑝𝑆𝑝𝑒𝑐, 500, 650)

𝑐𝑜𝑜𝑟𝑑𝑠2 = [
𝐺𝑑 𝑆𝑑
𝐺𝑒 𝑆𝑒
𝐺𝑧2 𝑆𝑧2
]

𝑝𝑒𝑟𝑐𝑒𝑛𝑡2 = [
𝑑
𝑒
]

(𝑐𝑜𝑜𝑟𝑑𝑠3, 𝑝𝑒𝑟𝑐𝑒𝑛𝑡3) = 𝑝ℎ𝑎𝑠𝑜𝑟𝑐𝑎𝑙𝑐𝑢𝑙𝑎𝑡𝑖𝑜𝑛𝑠(𝑠𝑝𝑒𝑐3, 𝑒𝑥𝑝𝑆𝑝𝑒𝑐, 650, 800)

𝑐𝑜𝑜𝑟𝑑𝑠3 = [
𝐺𝑓 𝑆𝑓
𝐺ℎ 𝑆ℎ
𝐺𝑧3 𝑆𝑧3
]

𝑝𝑒𝑟𝑐𝑒𝑛𝑡3 = [
𝑓
ℎ
]

𝑠𝑢𝑝𝑒𝑟 = [
𝐺𝑎𝑏𝑐 𝑆𝑎𝑏𝑐
𝐺𝑑𝑒 𝑆𝑑𝑒
𝐺𝑓ℎ 𝑆𝑓ℎ
]

𝐺𝑎𝑏𝑐 = 𝑎 · 𝐺𝑎 + 𝑏 · 𝐺𝑏 + 𝑐 · 𝐺𝑐
*the other G and S values follow the same pattern
*the formation of super is not included in the program and must be done by the user

(𝑐𝑜𝑜𝑟𝑑𝑠4, 𝑝𝑒𝑟𝑐𝑒𝑛𝑡4) = 𝑝ℎ𝑎𝑠𝑜𝑟𝑐𝑎𝑙𝑐𝑢𝑙𝑎𝑡𝑖𝑜𝑛𝑠(𝑠𝑢𝑝𝑒𝑟, 𝑒𝑥𝑝𝑆𝑝𝑒𝑐, 300, 800)

𝑐𝑜𝑜𝑟𝑑𝑠4 = [
𝐺𝑎𝑏𝑐 𝑆𝑎𝑏𝑐
𝐺𝑑𝑒 𝑆𝑑𝑒
𝐺𝑓ℎ
𝐺𝑧4
𝑆𝑓ℎ
𝑆𝑧4
]

𝑝𝑒𝑟𝑐𝑒𝑛𝑡4 = [
𝑎𝑏𝑐
𝑑𝑒
𝑓ℎ
]

𝑡𝑜𝑡𝑎𝑙𝑐𝑜𝑛𝑡𝑟𝑖𝑏𝑢𝑡𝑖𝑜𝑛𝑠𝑎 = 𝑎𝑏𝑐 · 𝑎 * the other total contribtion values follow the same pattern
