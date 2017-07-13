         
nudel o_180020.trl
nudel o_180020.lnz
nudel o_180020.tmp
nudel o_180020.lev
nudel o_180120.trl
nudel o_180120.lnz
nudel o_180120.tmp
nudel o_180120.lev
nudel o_180220.trl
nudel o_180220.lnz
nudel o_180220.tmp
nudel o_180220.lev
nudel o_180320.trl
nudel o_180320.lnz
nudel o_180320.tmp
nudel o_180320.lev
nudel o_180420.trl
nudel o_180420.lnz
nudel o_180420.tmp
nudel o_180420.lev
         
nuaddint o_180bb.addint o_180.int o_180.ant
nucp  o_180.ppar p.par
nucp  o_180.npar n.par
nudel o_180sl.inf
NuShellX < o_180.modelx > o_180mod.out
NuShellX < o_180.levelx > o_180lev.out
nucp nul > opd.dat
source shellx.bat > o_180.cpu
nuren shellx.bat shellx1.bat
nudel opd.dat
nucp  bb0200.lph o_180020.eng o_180.modelx o_180.levelx bb0200.lpe > nucp.txt
nuren o_180020.lp bb0200.lp
nuren o_180020.ls bb0200.ls
nuren o_180020.nhw bb0200.nhw
nucp  bb0202.lph o_180120.eng o_180.modelx o_180.levelx bb0202.lpe > nucp.txt
nuren o_180120.lp bb0202.lp
nuren o_180120.ls bb0202.ls
nuren o_180120.nhw bb0202.nhw
nucp  bb0204.lph o_180220.eng o_180.modelx o_180.levelx bb0204.lpe > nucp.txt
nuren o_180220.lp bb0204.lp
nuren o_180220.ls bb0204.ls
nuren o_180220.nhw bb0204.nhw
nucp  bb0206.lph o_180320.eng o_180.modelx o_180.levelx bb0206.lpe > nucp.txt
nuren o_180320.lp bb0206.lp
nuren o_180320.ls bb0206.ls
nuren o_180320.nhw bb0206.nhw
nucp  bb0208.lph o_180420.eng o_180.modelx o_180.levelx bb0208.lpe > nucp.txt
nuren o_180420.lp bb0208.lp
nuren o_180420.ls bb0208.ls
nuren o_180420.nhw bb0208.nhw
nulev o_18b     
levp o_18b     
