# README

### Version 1 BDSKY -  : 

	- change times 3.3 6.6 (xml is prelim_gastruloid_8targetBCs_skyline.xml)
	- origin 10

### Version 2 BDSKY: 

    - (NOT in repo: xml is identical to prelim_gastruloid_8targetBCs_skyline.xml, except for origin and changetimes)
    - origin 14
	- change times 7.0 9.0 
	- meant to mimick clonal formation, aggregate formation, elongation

### Version 3 BDSKY: 3-mGASv2-skyline-ou.xml

CORRECT timeline as shared by Sam (see Below)

 In total, 22 days of experiment, of which 10 are the gastruloid growth itself? (D12 -> D22 if clonal expansion part of the process?)

	- Transfection of DNA Typewriter machinery and stable integration: D0-D11
	- From single cell to clonal expansion: D12-D15 (1 cell)
	- Aggregate growth upstream of chiron treatment: D16-D18 
	- Chiron treatment: D19 (300 cells)
	- Symmetry breaking and elongation: D20-D22 (30K/ gastruloid)

change times: 0 4 7 8

Assuming that during the day of chiron treatment, there may be a
separate birth and death rate. 


### Version 4 BDSKY: 4-mGASv2-skyline-ou.xml


change times: 0 4 7.5

Chiron treatment begins at day 7. Here we assume that for half a day
nothing changes to before and then at half a day the rates start to
change. Effectively, this means we have 2 parameters less to estimate
and are a bit less prone to overfitting anything happening during a 1
day time bin (as done in version 3). 
