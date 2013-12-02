SMARTS validation
-----------------

The smarts_validation.py script is used to validate SMARTS matching. The examples
wre taken from http://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html

One SMARTS (i.e. [H]) is excluded since the data/*.hel molecule files do not
contian explicit hydrogens.

These are excluded since OpenBabel does not support multi-component SMARTS queries:

[NX3;H2,H1;!$(NC=O)].[NX3;H2,H1;!$(NC=O)]
[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8].[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]
[#16X2H0][!#16].[#16X2H0][!#16]
[F,Cl,Br,I].[F,Cl,Br,I].[F,Cl,Br,I]
[$([*R2]([*R])([*R])([*R]))].[$([*R2]([*R])([*R])([*R]))]
[cR1]1[cR1][cR1][cR1][cR1][cR1]1.[cR1]1[cR1][cR1][cR1][cR1][cR1]1
[#16X2H0][!#16].[#16X2H0][!#16]
([Cl!$(Cl~c)].[c!$(c~Cl)])
([Cl]).([c])
([Cl].[c])
[NX3;H2,H1;!$(NC=O)].[NX3;H2,H1;!$(NC=O)]
([!-0!-1!-2!-3!-4].[!+0!+1!+2!+3!+4])


SMARTS parser error in helium: [N&X4&+,N&X3&+0]
