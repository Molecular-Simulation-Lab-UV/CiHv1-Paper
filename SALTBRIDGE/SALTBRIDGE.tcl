## INPUT ##
### Run example :
# vmd -e SaltBridge.tcl
# 


#Input PSF
set psf /home/dbasaez/NAS-data/RESPALDO_BRUTO_CIHV1_VOLTAGE/NVT150mV/WT/MD1/Hv1_MD1_WT-FA_400ns.psf;
set reference_pdb /home/dbasaez/NAS-data/RESPALDO_BRUTO_CIHV1_VOLTAGE/NVT150mV/WT/MD1/Hv1_MD1_WT-FA_400ns.pdb; # pdb with centered channel aligned in z for alignment
#Number of input dcd
set firstDCD 8
set lastDCD  9
#inputs dcd (asumes that dcd start with "eq")
for {set i $firstDCD} {$i <= $lastDCD} {incr i} {
    set dcd($i) /home/dbasaez/NAS-data/RESPALDO_BRUTO_CIHV1_VOLTAGE/FA_WT_MD1_4us/eq$i.center.gather.dcd
}

#Path to bigdcd script
set bigdcd /home/jgarate/work/HV1/MODELS_Ci-HV1_paperCG_with_S0/MDs/ANA/bigdcd.tcl
#Salt Bridges Descriptors

set first 0; #First snapshot

#set E for N264E, R for N264R and anything for WT
set system WT

#atoms for residues  are:
#For Glutamics (E) name CD
#For Aspartics (D) name CG
#For Arginines (R) name CZ
#For Histidines (H) name CE1
#For Lysines   (K) name NZ

#For WT

set descriptors {FIT D171-R255 D171-R258 D160-R261 D233-R261 D171-K173 D222-K273 D222-K271 D222-K205 D233-H188 E167-R258 E167-R261 E219-R137 E201-R137 E219-K205 E201-K205 E185-H188 E167-H188}
set names(D171-R255) "(protein and resid 171 and name CG) or (protein and resid 255 and name CZ)"; #
set names(D171-R258) "(protein and resid 171 and name CG) or (protein and resid 258 and name CZ)"; #
set names(D160-R261) "(protein and resid 160 and name CG) or (protein and resid 261 and name CZ)"; #
set names(D233-R261) "(protein and resid 233 and name CG) or (protein and resid 261 and name CZ)"; #
set names(D171-K173) "(protein and resid 171 and name CG) or (protein and resid 173 and name NZ)"; #
set names(D222-K273) "(protein and resid 222 and name CG) or (protein and resid 273 and name NZ)"; #
set names(D222-K271) "(protein and resid 222 and name CG) or (protein and resid 271 and name NZ)"; #
set names(D222-K205) "(protein and resid 222 and name CG) or (protein and resid 205 and name NZ)"; #
set names(D233-H188) "(protein and resid 233 and name CG) or (protein and resid 188 and name CE1)"; #
set names(E167-R258) "(protein and resid 167 and name CD) or (protein and resid 258 and name CZ)"; #
set names(E167-R261) "(protein and resid 167 and name CD) or (protein and resid 261 and name CZ)"; #
set names(E219-R137) "(protein and resid 219 and name CD) or (protein and resid 137 and name CZ)"; #
set names(E201-R137) "(protein and resid 201 and name CD) or (protein and resid 137 and name CZ)"; #
set names(E219-K205) "(protein and resid 219 and name CD) or (protein and resid 205 and name NZ)"; #
set names(E201-K205) "(protein and resid 201 and name CD) or (protein and resid 205 and name NZ)"; #
set names(E185-H188) "(protein and resid 185 and name CD) or (protein and resid 188 and name CE1)"; #
set names(E167-H188) "(protein and resid 167 and name CD) or (protein and resid 188 and name CE1)"; #

if {$system eq "R"} {
    #For N264R
    set descriptors {FIT D160-R258 D160-R261 D160-R255 D222-R264 E201-R264 D233-R264 D160-R264}
    set names(D160-R258) "(protein and resid 160 and name CG) or (protein and resid 258 and name CZ)"; #
    set names(D160-R261) "(protein and resid 160 and name CG) or (protein and resid 261 and name CZ)"; #
    set names(D160-R255) "(protein and resid 160 and name CG) or (protein and resid 255 and name CZ)"; #
    # Posible salt bridges of R264 with neighbour negatives residues
    set names(D222-R264) "(protein and resid 222 and name CG) or (protein and resid 264 and name CZ)"; #
    set names(E201-R264) "(protein and resid 201 and name CD) or (protein and resid 264 and name CZ)"; #
    set names(D233-R264) "(protein and resid 233 and name CG) or (protein and resid 264 and name CZ)"; #
    set names(D160-R264) "(protein and resid 160 and name CG) or (protein and resid 264 and name CZ)"; #
}
if {$system eq "E"} {
    #For N264E
    set descriptors {FIT D160-R258 D160-R261 D160-R255 E264-R261 E264-R258 E264-R255}
    set names(D160-R258) "(protein and resid 160 and name CG) or (protein and resid 258 and name CZ)"; #
    set names(D160-R261) "(protein and resid 160 and name CG) or (protein and resid 261 and name CZ)"; #
    set names(D160-R255) "(protein and resid 160 and name CG) or (protein and resid 255 and name CZ)"; #
    # Posible salt bridges of E264  with neighbour positive residues
    set names(E264-R261) "(protein and resid 264 and name CD) or (protein and resid 261 and name CZ)"; #
    set names(E264-R258) "(protein and resid 264 and name CD) or (protein and resid 258 and name CZ)"; #
    set names(E264-R255) "(protein and resid 264 and name CD) or (protein and resid 255 and name CZ)"; #
}




set names(FIT) "protein and resid 160 193 232 258 and name CA"; #Change Depending pore selection

## END INPUT ##
#################################################################################################################################### 
### Do not Change!!!
#init array for selections
array set indexes {}
set steps 0

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

#Computes euclidian distance from a list of two atom indexes
proc distance { distance } {
    set index [lindex $distance 0]
    set coord1 [measure center [atomselect top "index $index"]]
    set index [lindex $distance 1]
    set coord2 [measure center [atomselect top "index $index"]]
    set distance [veclength [vecsub $coord1 $coord2]]
    return $distance
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
####
proc RunAna {descriptors &arrName } {
    upvar 1 ${&arrName} indexes
    foreach  descriptor $descriptors {
	if { [llength $indexes($descriptor)] > 2} {
	    set results($descriptor)  [RMSD $indexes(FIT)]
	}
	if { [llength $indexes($descriptor)] == 2} {
	    set results($descriptor) [distance $indexes($descriptor)]
	} 
    }
    return [array get results]
}


# Procedure to be run with bigdcd
proc RunBigDCD {frame} {
    global names indexes first steps descriptors
    if {$steps >= $first } {
	array set results [RunAna   $descriptors indexes]
	WriteResult $descriptors "Distances.dat" $steps results
    }
    incr steps
}

## END PROCEDURES ##

## MAIN ##
proc main {&arrName } {
    upvar 1 ${&arrName} dcd 
    global psf bigdcd reference_pdb names
    global firstDCD lastDCD descriptors indexes
    WriteInit $descriptors "Distances.dat"
    mol load psf $psf
    animate read pdb $reference_pdb
    #Do selections once
    array set indexes [SetIndex $descriptors names]
    source $bigdcd
    for {set i $firstDCD} {$i <= $lastDCD} {incr i} {
    	bigdcd RunBigDCD $dcd($i)
    	bigdcd_wait
    }
}

#### RUN ##
main dcd
exit
