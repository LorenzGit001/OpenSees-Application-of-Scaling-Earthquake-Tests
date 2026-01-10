# source MDOF_elcentro_Froude_artificialmass.tcl
# El Centro 
# Three story frame model with Froude Similitude
# same Material E, rho steel 
# L=	6,40E-02m    one story 
# m1=m2=m3=	1,12E-01 kg
# h=	4,00E-04	m
# b=	2,54E-02	m
# E=	1,96E+11	N/m^2
# I=	1,35E-13	m^4 one colum 
#           2688 steps * 0.02s = 53.76s 
#     simplefied lumped colum model   
# first version i put the massen on top of node so 64mm
#            x node 4 
#            |   
#            |
#            x node 3 
#            | 
#.           |
#            x node 2
#            |
#            |
#     _______x__node 1_
#   
#.  10.4170.  Eigenfrequency calculated in matlab 
#   29.1878
#   42.1776
# 
wipe;
model basic -ndm 2 -ndf 3;
file mkdir data_MODF_Froude_mass;

# --- Scale factors Froude similtude ---
set Ls 0.7   ;# length scale (e.g. 1/100)
set Ms [expr $Ls**3] ;# mass scale (rho constant)
set As [expr $Ls**2] ;# area scale
set Is [expr $Ls**4] ;# inertia scale
# dummy/artificial mass increase (1.4286 re-calculated through the desired eingenfrequencies) 
  # extra time scaling will distort the kinematic Froude conditions!
  # acceleration scale should be = 1 or close to one 
set gamma_m 1.4286;
# for sqrt(Ls/($gamma_m) considering the increase 
# or sqrt(Ls) no extra time scaling!
set Ts [expr sqrt($Ls)] ; 


# --- Node coordinates ( )
node 1 0.0 0.0
node 2 0.0 [expr 0.064*$Ls]
node 3 0.0 [expr 0.128*$Ls]
node 4 0.0 [expr 0.192*$Ls]


# fix node one in all dof 
fix 1 1 1 1
# fix 2 0 1 1
# fix 3 0 1 1
# fix 4 0 1 1

# fix and equal DOF seem to result in the same 
# no vertical displacemnts and rotations at nodes) end nodes of the beams are constrained to each other in the 2nd DOF (vertical displacement) and the 3rd DOF (rotation)
equalDOF 1 2 2 3 
equalDOF 2 3 2 3 
equalDOF 3 4 2 3


# --- Nodal mass (1,12E-01 kg * Ls^3) ---
mass 2 [expr 1.12e-01*$Ms*$gamma_m]  0.0 0.0
mass 3 [expr 1.12e-01*$Ms*$gamma_m]  0.0 0.0
mass 4 [expr 1.12e-01*$Ms*$gamma_m]  0.0 0.0


# --- Geometric transformation 
geomTransf Linear 1

# --- Element properties (prototype values) ---
set A0 [expr {1.02E-05*2}]  ;# singe 1,02E-05 , double (frame)2.0320e-05 full-scale cross-sectional area [m^2]
set E0 1.96e+11  ;# modulus (unchanged) N/m^2
set Iz0 1.35e-13 ;#  2.71e-13 Iz [m^4] double the area 




# Apply scaling
set A  [expr $A0*$As]
set Iz [expr $Iz0*$Is]
set dt_ground_motion [expr 0.02*$Ts]; 
set dt_analysis [expr 0.02*$Ts];

# combined inertia for two beams 
set Icomb [expr {2.0 * $Iz}] ;  # = 2.71e-13

# -- Element with material prop --
# compute ks and add translational springs between nodes in DOF1
set Lseg [expr {0.064*$Ls}] ;# segment length used in MATLAB
set ks [expr {12.0 * $E0 * $Icomb / ($Lseg**3)}] ; # using Icomb for double area

#uniaxialMaterial Elastic 1 $ks    #can work with ks 
#element zeroLength 1 1 2 -mat 1 -dir 1
#element zeroLength 2 2 3 -mat 1 -dir 1
# element zeroLength 3 3 4 -mat 1 -dir 1


# -- elements -- 
 element elasticBeamColumn 1 1 2 $A $E0 $Icomb 1
 element elasticBeamColumn 2 2 3 $A $E0 $Icomb 1
 element elasticBeamColumn 3 3 4 $A $E0 $Icomb 1

# --- Recorders ---
recorder Node -file data_MODF_Froude_mass/mass134FstoryNode2_disp_L07.txt -time -node 2 -dof 1 2 3 disp
recorder Node -file data_MODF_Froude_mass/mass134FstoryNode3_disp_L07.txt -time -node 3 -dof 1 2 3 disp
recorder Node -file data_MODF_Froude_mass/mass134FstoryNode4_disp_L07.txt -time -node 4 -dof 1 2 3 disp
# recorder Node -file Data4/Gomes_10th_dispNode4_TOG.txt -time -node 4 -dof 1 2 3 disp
# recorder Node -file Data2/Node2Accv7_5_centro_100th_eig.txt -time -node 2 -dof 1 accel
# recorder Node -file Data/Node1Accvv7_5_centro_100th_eig.txt -time -node 1 -dof 1 accel
# recorder Node -file Data/RBasevv7_5_centro_100th_eig.txt -time -node 1 -dof 1 2 3 reaction
# recorder Drift -file Data/Driftvv7_5_centro_100th_eig.txt -time -iNode 1 -jNode 2 -dof 1 -perpDirn 2
# recorder Element -file Data/FColvv7_5_centro_100th_eig.txt -time -ele 1 force
# recorder Element -file data_MODF_Froude/Fele1_forceL07.txt -time -ele 1 force
recorder Node -file data_MODF_Froude_mass/mass134Fvelocity_L07.txt -time -nodeRange 2 4 -dof 1 2 3 vel
recorder Node -file data_MODF_Froude_mass/mass134Fbase_shear07.txt -time -node 1 -dof 1 reaction

# --- GRAVITY load ---
# load is not written for this script
#   kg * 9.81 on node 2
set G 9.81
timeSeries Linear 1
pattern Plain 1 1 {
    load 2 0.0 [expr -1.12e-01*$Ms*$G] 0.0
    load 3 0.0 [expr -1.12e-01*$Ms*$G] 0.0
    load 4 0.0 [expr -1.12e-01*$Ms*$G] 0.0
}
constraints Plain
numberer Plain
system BandGeneral
algorithm Linear
integrator LoadControl 0.1
analysis Static
analyze 10
loadConst -time 0.0

# --- Eigenvalue and Eigenfrequency Analysis ---
# Description: This section calculates the natural frequencies of the structure in Hertz (Hz).

puts "Starting Eigenvalue Analysis..."

# Set the number of eigenvalues (modes)
set numModes 3

# The 'eigen' command returns a Tcl list of eigenvalues (lambda = omega^2)
set eigenvalues [eigen -fullGenLapack $numModes]

puts "-------------------------------------"
puts "Eigenvalue and Eigenfrequency Results"
puts "-------------------------------------"

# Open a file to store the results
set resultsFile [open "Eigenfrequencies_Hz.txt" "w"]
puts $resultsFile "Mode\tEigenvalue (rad²/s²)\tFrequency (Hz)"

# Loop through each of the calculated eigenvalues
for {set i 0} {$i < $numModes} {incr i} {

    # Get the eigenvalue for the current mode from the list
    set lambda [lindex $eigenvalues $i]

    # Eigenvalues for stable structures should be positive.
    if {$lambda > 0} {
        # Calculate the circular frequency omega (in rad/s) by taking the square root
        set omega [expr {sqrt($lambda)}]

        # Calculate the natural frequency f (in Hz) by dividing by 2*pi
        # Using [expr {acos(-1.0)}] is a precise way to get the value of pi in Tcl
        set frequencyHz [expr {$omega / (2.0 * acos(-1.0))}]

        # Print the results to the console for immediate feedback
        puts "Mode [expr {$i + 1}]: Eigenvalue = ${lambda}, Frequency = ${frequencyHz} Hz"

        # Write the formatted results to the output file
        puts $resultsFile "[expr {$i + 1}]\t$lambda\t$frequencyHz"

    } else {
        # Handle cases with non-positive eigenvalues, which can indicate model instability
        puts "Mode [expr {$i + 1}]: Non-positive eigenvalue found ($lambda). Cannot calculate frequency."
        puts $resultsFile "[expr {$i + 1}]\t$lambda\tNaN"
    }
}

# Close the output file
close $resultsFile

puts "-------------------------------------"
puts "Eigenfrequency analysis complete. Results saved to Eigenfrequencies_Hz.txt"
# ----------------------------------------------------


# DYNAMIC ground-motion analysis
set G 9.81
# El Centro record: 2688 points, dt = 0.02s original
timeSeries Path 2 -dt [expr 0.02*$Ts] -filePath ElCentro_acc.tcl -factor $G;
pattern UniformExcitation 2 1 -accel 2;

# damping via first eigenmode
# -- Gomes et. al calculated with no damping --
set freq [expr [eigen -fullGenLapack 1]**0.5]
set dampRatio 0.05
rayleigh 0. 0. 0. [expr 2*$dampRatio/$freq] 


wipeAnalysis;
constraints Plain;
numberer Plain;
system BandGeneral;
algorithm Linear
integrator Newmark 0.5 0.25 ;
analysis Transient;
# total time = 53.74s -> 
# alternative analyze every 2nd step 1344Pts and dt*2 (=0.04s)
# is not necessery 
analyze 2688 $dt_analysis;
# analyze 2688 0.02; 

puts "Length scale Ls = $Ls"
puts "Time scale Ts = Ls = $Ts"
puts "Prototype earthquake dt = 0.02 s"
puts "Model ground motion dt = 0.02×Ts = [expr 0.02*$Ts] s"
puts "analysis dt = $dt_analysis s"
puts "Mismatch (>1<) = [expr $dt_analysis/(0.02*$Ts)]"
puts "Time compression factor (prototype/model) = 1/Ts = [expr 1.0/$Ts]"
puts "El Centro has been shaking!"

wipe

