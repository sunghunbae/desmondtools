# 
#  usage: 
#  vmdcli vmd_desmond.tcl -args md/r01/desmond_md_job_...-out.cms
#

package require pbctools

set outcms [lindex $argv 0]
set incms [string map {\-out.cms \-in.cms} $outcms]
set outtrj [string map {\-out.cms _trj} $outcms]

# tcl script
# file tail <filename>
# file rootname <filename>

set outpdb [string map {\-out.cms \-out.pdb} $outcms]
set outpdb [file tail $outpdb ]

set outdrydcd [string map {\-out.cms \-out-nowat.dcd} $outcms]
set outdrydcd [file tail $outdrydcd ]

set outdrypdb [string map {\-out.cms \-out-nowat.pdb} $outcms]
set outdrypdb [file tail $outdrypdb ]

set inpdb [string map {\-in.cms \-in.pdb} $incms]
set inpdb [file tail $inpdb]

set indrypdb [string map {\-in.cms \-in-nowat.pdb} $incms]
set indrypdb [file tail $indrypdb]

puts ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
puts "convert desmond  in trajectory into pdb"

mol new $incms
set sys [atomselect top all]
set off [vecsub {0 0 0} [measure center $sys]]
$sys moveby $off
$sys writepdb $inpdb
set dry [atomselect top {not (resname T3P NA CL)}]
$dry writepdb $indrypdb
mol delete all

puts ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
puts "convert desmond out trajectory into dcd"
puts ""
puts "cms= $outcms"
puts "trj= $outtrj"

mol new  $outcms
mol addfile $outtrj type {dtr} first 0 last -1 step 1 waitfor all top

set lastframe [ expr [ molinfo top get numframes ] - 1 ]
puts "lastframe= $lastframe"


#pbc wrap -centersel "not (resname T3P NA CL)" -center bb -compound chain -all
#pbc wrap -centersel "chain A" -center unitcell -compound chain -all
pbc wrap -compound segid -all
#pbc unwrap -sel "not (resname T3P NA CL)" -all


set sys [atomselect top all]
set off [vecsub {0 0 0} [measure center $sys]]
$sys moveby $off
puts "translated= $off"

$sys writepdb $outpdb
puts "written to $outpdb"

set ref [atomselect top "chain A" frame 0]
set sel [atomselect top "chain A"]
set all [atomselect top all]
for { set i 1 } { $i <= $lastframe } { incr i } {
    $sel frame $i
    $all frame $i
    $all move [ measure fit $sel $ref ]
}

#set dry [atomselect top {not (resname T3P NA CL)}]
set dry [atomselect top {not (resname SPC T3P NA CL)}]

$dry writepdb $outdrypdb
puts "nowater written to $outdrypdb"


animate write dcd $outdrydcd beg 0 end $lastframe skip 1 waitfor all sel $dry top
puts "nowater written to $outdrydcd"

exit
