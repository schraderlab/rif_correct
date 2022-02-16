Columns of the tab delimited output file
A.	Gene name
B.	Slope of the linear regression from the natural log transformed data
C.	Y intercept of the linear regression from the natural log transformed data
D.	R value of the linear regression from the natural log transformed data
E.	Half-Life Calc (minutes)
F.	Number of data pts. After user inputted cut off was applied
G.	Number of leading data points affected by lag calc, if no poly elongation rate then “0”, if any rate then at least “-1”
H.	Natural Log adjusted value of last data point used in the regression
I.	In operon (y or n): if you did not use the optional transcriptome file, will always be “n”
J.	Operon number
K.	Complex or simple: if you did not use the optional transcriptome file, will always be “simple”
L.	Total number of data points inputted for that particular gene
M.	R value Flag: if the r value of the linear regression is above user inputted r value for a gene, this flag is raised
N.	Half-Life Flag: if the half-life is more than twice as long as the last time point used in the regression, or the half-life is negative (indicating the level of mRNA is not decreasing, but rather increasing) this flag is raised

Columns 0. Through AH. Are headered by the user inputted “time points” and for each gene will report the natural log transformed value of your data, including data below cut off

