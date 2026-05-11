set psf ../HV1.POPC.Wat.box.ion.ModelN264RFA_4A.hmr.psf 
set dcd ../eq6.gather.dcd
set outname PMEPOT.dx
set xsc ../eq6.xsc
set eFact 0.257952
set Res 1.0
#set frames 200:10:1200
set frames "all"
proc main { } {
    package require pmepot
    global psf dcd xsc eFact
    global outname frames Res
    
    mol load psf $psf
    mol addfile $dcd waitfor all
    pmepot -mol top -frames $frames -xscfile $xsc -ewaldfactor $eFact -grid $Res -dxfile $outname
}
main


