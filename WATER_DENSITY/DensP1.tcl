## INPUT ##
### Run example :
# vmd -e WaterDensP1.tcl
# 


#Input PSF
set psf HV1.POPC.Wat.box.ion.ModelA.hmr.psf ;
#set reference_pdb ../HV1.POPC.Wat.box.ion.ModelIR.160.ASN.pdb ; # pdb with centered channel aligned in z for alignment
set reference_pdb HV1.POPC.Wat.box.ion.ModelA.hmr.pdb ; # pdb with centered channel aligned in z for alignment
#Number of input dcd
set firstDCD 5
set lastDCD  9
#inputs dcd (asumes that dcd start with "eq")
for {set i $firstDCD} {$i <= $lastDCD} {incr i} {
    #set dcd($i) 0.125.gather.dcd
    set dcd($i) ../eq$i.gather.dcd
}

#Path to bigdcd script
set bigdcd /home/jgarate/work/HV1/MODELS_Ci-HV1_paperCG_with_S0/MDs/ANA/bigdcd.tcl
#Dipole Descriptors

set first 0; #First snapshot

set binNumZ  60;  #Numbers of bins along axial Dim
set binNumR  15;   #Numbers of bins along radial Dim

set minZ -30; #Min Z value for pore
set maxZ  30; #Max Z value for pore


set rad 10; #Radius of pore
set WatNum  17314; # Total Number of waters
set BoxVol 548576; #in A^3 box volume minus membrane volume
set MaxLoad 300;# Could be any number >= Minimum Loads of pore

set descriptors {SEL FIT}

set names(SEL) "name OH2"; #
set names(FIT) "chain H and resid 160 255 258 261 264"; #Change Depending pore selection

## END INPUT ##
#################################################################################################################################### 
### Do not Change!!!
set steps 0
array set indexes {}
set rad2 [expr $rad*$rad]
set PI  3.14159265359;
# Bin sizes 
# Z: axial bins
# R: Radial bins
set BinSizeZ  [expr 1.000*($maxZ-$minZ)/$binNumZ]
set BinSizeR  [expr 1.000*($rad)/$binNumR]
# Set arrays of P1 values along z axis
# And loads along z axis
for { set i 0}  {$i < $binNumZ} {incr i} {
    set P1Axis($i) 0
    set LoadsCounter($i) 0
    for { set j 0}  {$j < $binNumR} {incr j} {
	set AxialRad($i,$j) 0
    }
}
## PROCEDURES ##
#Transforms name selections into indexes selections
proc SetIndex { descriptors  &arrName } {
    upvar 1 ${&arrName} names
    foreach  descriptor $descriptors {
	#bad if descriport varies with order (e.g. torsion)
	set temp [ [atomselect top "$names($descriptor)"] get index]
	set indexes($descriptor) $temp
	#$temp delete
    }
    return [array get indexes]
}
#Structural Fit against a reference structure
#reference_pdb
#selection (selection indexes) 
proc RMSD {selection } {
    #set ref [atomselect top "index $selection" frame 0]
    #align against the same reference
    #frame 0 contains a centred chain A, same reference for all monomers
    #Only TMD domain
    #Change to employ same referece for all momomers
    #set ref [atomselect top "chain A and resid 210 213" frame 0]
    set ref [atomselect top "index $selection" frame 0]
    set sel [atomselect top "index $selection"]
    set all [atomselect top all]
    $all move [measure fit $sel $ref]
    set rmsd [measure rmsd $sel $ref]
    $ref delete
    $sel delete
    $all delete
    return $rmsd
} 
#Main Procedure that collects al data of MD 
#trajectory
proc SelecCollect { selection &arrName1 &arrName2 &arrName3 } {
    global minZ maxZ rad rad2; #Binning parameters
    global BinSizeZ BinSizeR;  #Binning parameters
    global binNumZ binNumR;    #Binning parameters
    upvar 1 ${&arrName1}  LoadsCounter;# Total Counter of observations binned in z
    upvar 1 ${&arrName2}  AxialRad;# Array of list in 2D for binnig Axial Radial histograms
    upvar 1 ${&arrName3}  P1Axis;
    set WatLoadCounter 0
    #loop all water molecules
    set indexes [ [atomselect top "index $selection and z> $minZ and z < $maxZ and (x^2 + y^2)< $rad2"] get index]
    set WatLoadCounter [llength $indexes]
    foreach index $indexes {
	# Do selections and collect vectors
	set indexH1 [expr $index +1]
	set indexH2 [expr $index +2]
	set sel [atomselect top "index $index $indexH1 $indexH2"]
	set coord [$sel get {x y z}]
	set x [lindex [lindex $coord 0] 0]
	set y [lindex [lindex $coord 0] 1]
	set z [lindex [lindex $coord 0] 2]
	set rDist2 [expr $x*$x + $y*$y]
	# Binning
	set stepZ [expr int(( ($z-$minZ)/($BinSizeZ) ))]
	set rDist [expr sqrt($rDist2)]
	set stepR [expr int(($rDist)/($BinSizeR))]
	# Accumulate Axial Radial histogram
	if { $stepZ < $binNumZ && $stepZ >= 0  && $stepR < $binNumR } {
	    #incr WatLoadCounter
	    incr AxialRad($stepZ,$stepR)
	    #Accumulate P1 values, binned along z axis
	    set vectorDip [vecnorm [measure dipole $sel -masscenter]]
	    set P1Axis($stepZ) [expr $P1Axis($stepZ) +[lindex $vectorDip 2]]
	    incr LoadsCounter($stepZ); # and accumulate for averages
	}
	$sel delete
    }
    return $WatLoadCounter;
}

#Ubin Average P1 along axial dim
proc UnbinDipZ { &arrName1 &arrName2 } {
    upvar 1 ${&arrName1} P1Axis
    upvar 1 ${&arrName2} LoadsCounter
    global binNumZ BinSizeZ minZ
    set out [open "P1_axial.dat" w]
    puts $out "#AvgP1 along z axis"
    puts $out "#z          <P1>"
    for { set i 0}  {$i < $binNumZ} {incr i} { 
	set z [format {%8.2f} [expr $i*$BinSizeZ + $minZ + $BinSizeZ*0.5]]
	if { $LoadsCounter($i) > 0} {
	    set AverageP1 [format {%8.2f} [expr 1.00*$P1Axis($i)/$LoadsCounter($i)] ]
	} else {
	    set AverageP1 [format {%8.2f} 0 ]
	}
	puts -nonewline $out $z
	puts -nonewline $out $AverageP1
	puts $out ""
    }
    close $out
}
#Ubin Average Load along axial dim
proc UnbinLoadZ {&arrName } {
    upvar 1 ${&arrName} LoadsCounter
    global steps minZ BinSizeZ binNumZ
    ######
    set out [open "Loads_axial.dat" w]
    puts $out "#AvgLoad along z axis"
    puts $out "#z           <Load>"
    for { set i 0}  {$i < $binNumZ} {incr i} { 
	set z [format {%8.2f} [expr $i*$BinSizeZ + $minZ + $BinSizeZ*0.5]]
	set AverageLoad [format {%8.2f} [expr 1.00*$LoadsCounter($i)/$steps] ]
	puts -nonewline $out $z
	puts -nonewline $out $AverageLoad
	puts $out ""
    }
    close $out
}
#Unbin Axial Radial 2D arrays
proc AxRadUnbin {&arrName} {
    upvar 1 ${&arrName} AxialRad
    global steps WatNum Volume
    global minZ maxZ rad PI BoxVol
    global BinSizeZ BinSizeR
    global binNumZ binNumR
    # Set normalization factor for Cylindrical shells 
    set Total_density [expr 1.00*$WatNum/$BoxVol];
    set const_factor  [expr $PI*$BinSizeR*$BinSizeR];
    set norm_factor [expr $Total_density*$const_factor*$steps];
    # Write output
    set out [open "AxialRadDens.dat" w]
    puts $out "#Axial Radial density"
    puts $out "#Z          RAD      P/P0"
    for {set i  0} {$i < $binNumZ} {incr i} {
	set Z [format {%8.2f} [expr $i*$BinSizeZ + $minZ + $BinSizeZ*0.5]]
	for {set j  0} { $j < $binNumR} {incr j} {
	    set RAD [ format {%8.2f} [expr $j*$BinSizeR + $BinSizeR*0.5]];
	    set shell [expr 2*$j + 1 ]
	    set radial_norm_fact [expr $norm_factor*$shell*$BinSizeZ]
	    set 2Density [format {%8.2f} [expr $AxialRad($i,$j)/$radial_norm_fact ]]
	    puts -nonewline $out $Z
	    puts -nonewline $out $RAD
	    puts -nonewline $out $2Density
	    puts $out ""
	}
    }
    close $out
}

#Write Initial file: The avoids deletion when loading multiple dcd in a for loop
proc WriteInit {descriptors outname} { 
    set out [open $outname w]
    puts -nonewline $out "# Frame "
    foreach  descriptor $descriptors {
	puts -nonewline $out [ format {%10s} $descriptor]
    }
    puts $out ""
    close $out
}
#Writes results into a single File
proc WriteResult {descriptors outname frame &arrName} {
    upvar 1 ${&arrName} results
    set out [open $outname a+]
    puts -nonewline $out [format {%10s} $frame ]
    foreach  descriptor $descriptors {
	puts -nonewline $out [format {%10.3f} $results($descriptor)]
    }
    puts $out ""
    close $out
} 

# Run analyses of defined descriptors
# arrays names must be declared global if
# a function within a functions uses them
# as arguments
proc RunAna {&arrName1 } {
    upvar 1 ${&arrName1} indexes
    global P1Axis LoadsCounter AxialRad
    set results(FIT) [RMSD $indexes(FIT)]
    set results(SEL) [SelecCollect $indexes(SEL) LoadsCounter AxialRad P1Axis]
    return [array get results]
}
# Procedure to be run with bigdcd
proc RunBigDCD {frame} {
    global names indexes first steps descriptors
    if {$steps >= $first } {
	array set results [RunAna  indexes ]
	WriteResult $descriptors "LoadsRMSD.dat" $steps results
    }
    incr steps
}
## END PROCEDURES ##

## MAIN ##
proc main {&arrName } {
    upvar 1 ${&arrName} dcd 
    global psf bigdcd reference_pdb
    global firstDCD lastDCD descriptors 
    global P1Axis LoadsCounter AxialRad
    global indexes names
    WriteInit $descriptors "LoadsRMSD.dat"
    mol load psf $psf
    animate read pdb $reference_pdb
    #Do index selections once
    array set indexes [SetIndex $descriptors names]
    source $bigdcd
    for {set i $firstDCD} {$i <= $lastDCD} {incr i} {
    	bigdcd RunBigDCD $dcd($i)
    	bigdcd_wait
    }
    #Unbin P1 along pore axis
    UnbinDipZ P1Axis LoadsCounter
    #Unbin loads along pore axis
    UnbinLoadZ LoadsCounter
    #Unbin P/P0 Radial Axial
    AxRadUnbin AxialRad
}

#### RUN ##
main dcd
exit
