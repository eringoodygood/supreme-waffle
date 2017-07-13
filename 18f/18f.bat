         
nudel f_180011.trl
nudel f_180011.lnz
nudel f_180011.tmp
nudel f_180011.lev
nudel f_180111.trl
nudel f_180111.lnz
nudel f_180111.tmp
nudel f_180111.lev
nudel f_180211.trl
nudel f_180211.lnz
nudel f_180211.tmp
nudel f_180211.lev
nudel f_180311.trl
nudel f_180311.lnz
nudel f_180311.tmp
nudel f_180311.lev
nudel f_180411.trl
nudel f_180411.lnz
nudel f_180411.tmp
nudel f_180411.lev
nudel f_180511.trl
nudel f_180511.lnz
nudel f_180511.tmp
nudel f_180511.lev
         
nuaddint f_180bb.addint f_180.int f_180.ant
nucp  f_180.ppar p.par
nucp  f_180.npar n.par

nudel f_180sl.inf
NuShellX < f_180.modelx > f_180mod.out
NuShellX < f_180.levelx > f_180lev.out
source shellx.bat > f_180.cpu
nuren shellx.bat shellx1.bat
nucp  bb1200.lph f_180011.eng f_180.modelx f_180.levelx bb1200.lpe > nucp.txt
nuren f_180011.lp bb1200.lp
nuren f_180011.ls bb1200.ls
nuren f_180011.nhw bb1200.nhw
nucp  bb1202.lph f_180111.eng f_180.modelx f_180.levelx bb1202.lpe > nucp.txt
nuren f_180111.lp bb1202.lp
nuren f_180111.ls bb1202.ls
nuren f_180111.nhw bb1202.nhw
nucp  bb1204.lph f_180211.eng f_180.modelx f_180.levelx bb1204.lpe > nucp.txt
nuren f_180211.lp bb1204.lp
nuren f_180211.ls bb1204.ls
nuren f_180211.nhw bb1204.nhw
nucp  bb1206.lph f_180311.eng f_180.modelx f_180.levelx bb1206.lpe > nucp.txt
nuren f_180311.lp bb1206.lp
nuren f_180311.ls bb1206.ls
nuren f_180311.nhw bb1206.nhw
nucp  bb1208.lph f_180411.eng f_180.modelx f_180.levelx bb1208.lpe > nucp.txt
nuren f_180411.lp bb1208.lp
nuren f_180411.ls bb1208.ls
nuren f_180411.nhw bb1208.nhw
nucp  bb120a.lph f_180511.eng f_180.modelx f_180.levelx bb120a.lpe > nucp.txt
nuren f_180511.lp bb120a.lp
nuren f_180511.ls bb120a.ls
nuren f_180511.nhw bb120a.nhw
nulev f_18b     
levp f_18b     
