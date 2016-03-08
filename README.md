This is IDL software that I developed my first year of grad school in order to separate the spectral components of a close unresolved M dwarf + White Dwarf (WD+dM) companion.

I have relied on this code numerous times throughout my graduate school career and yet, it hasn't been updated since I first wrote it. In my free time I have started to translate this code to Python as well as make upgrades. 

Upgrades:
1) Document and package.
2) Make more user friendly.
3) Upgrade the white dwarf models.
4) Implement markov-chain monte carlo to estimate uncertainties in parameter fits.

At present, it uses minimum chi-squared estimation to find the best-fit M dwarf template and best-fit white dwarf model. At present, the method doesn't match the two components simultaneously, but rather fits the white dwarf first, subtracts the best white dwarf fit from the raw spectrum, fits the M dwarf, subtracts the best M dwarf fit from the raw spectrum, and then fits the white dwarf again. This process iterates until it converges to consistent solutions for each the white dwarf and M dwarf or until it hits 10 iterations. Typically the method converges in less than five iterations.