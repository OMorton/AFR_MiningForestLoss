## Mining induced deforestation across sub-Saharan Africa

Code associated with Morton O.\*, Bousfield C G.\*, Dégny Valé P., Lamb I., Maus V., Bryant R G., & Edwards D P. (\*Joint first authors)  as **"Mining triggers extensive additional deforestation across sub-Saharan Africa"** at *Nature* (*Accepted in Principle*).

### File summary

**Functions.R**
Assorted helper functions to read in data, fit the heterogeneity robust DiD models, and convenience functions for plotting and summarising.

**01.national.DiD.R**
Fit country level DiD models to estimate the total additional deforestation (direct and indirect) associated with mine establishment. Across a suite of buffer sizes.

**02.SSA.wide.DiD.R**
Fit SSA wide DiD models to estimate the total additional deforestation (direct and indirect) associated with mine establishment. Across a suite of buffer sizes.

**03.direct.indirect.DiD.R**
Estimate additional direct and indirect deforestation associated with mine establishment, within a 0-5 km buffer.

**04.commodity.impact.R**
Using a subset of the data (see Methods) fit commodity specific DiD models to estimate the total additional deforestation (direct and indirect) associated with mine establishment. Across a suite of buffer sizes.

**05.main.text.figures.R**
Plotting and producing the main text figures and key variations of those included in the SM.