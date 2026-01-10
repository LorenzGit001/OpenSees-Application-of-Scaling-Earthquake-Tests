# source SDOF_elcentro_Froude.tcl
# VERSION v7.5 EL Centro
# Streetlight Model with Froude Similitude
# same Material E, rho
# damping \zeta should be 2% 

wipe;
model basic -ndm 2 -ndf 3;
file mkdir data_SDOF_Froude;

# --- Scale factors ---
set Ls 0.7   ;# length scale (e.g. 1/100)
set Ms [expr $Ls**3] ;# mass scale (rho constant)
set As [expr $Ls**2] ;# area scale
set Is [expr $Ls**4] ;# 2nd inertia scale
set Ts [expr sqrt($Ls)]; 


# --- Node coordinates (6m -> 6*$Ls m)
node 1 0.0 0.0
node 2 0.0 [expr 6.0*$Ls]

fix 1 1 1 1

# --- Nodal mass (20 kg * Ls^3) ---
mass 2 [expr 20.0*$Ms] 0.0 0.0

# --- Geometric transformation 
geomTransf Linear 1

# --- Element properties (prototype values) ---
set A0 0.017671    ;# full-scale cross-sectional area [m^2]
set E0 2.1e+11     ;# modulus (unchanged)
set Iz0 2.485e-5   ;# fullscale Iz [m^4]

# Apply scaling
set A  [expr $A0*$As]
set Iz [expr $Iz0*$Is]
set dt_ground_motion [expr 0.02*$Ts]; 
set dt_analysis [expr 0.02*$Ts];

# element 
element elasticBeamColumn 1 1 2 $A $E0 $Iz 1

# --- Recorders ---
recorder Node -file data_SDOF_Froude/Node2_disp_L07.txt -time -node 2 -dof 1 2 3 disp
# recorder Node -file data_SDOF_Froude/Node2_acc_L1.txt -time -node 2 -dof 1 accel
# recorder Node -file Data/Node1Accvv7_5_centro_100th_eig.txt -time -node 1 -dof 1 accel
# recorder Node -file Data/RBasevv7_5_centro_100th_eig.txt -time -node 1 -dof 1 2 3 reaction
# recorder Drift -file Data/Driftvv7_5_centro_100th_eig.txt -time -iNode 1 -jNode 2 -dof 1 -perpDirn 2
# recorder Element -file Data/FColvv7_5_centro_100th_eig.txt -time -ele 1 force

# --- GRAVITY load ---
set G 9.81
timeSeries Linear 1
pattern Plain 1 1 {
    load 2 0.0 [expr -20.0*$Ms*$G] 0.0
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
set numModes 1

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
set Ts [expr sqrt($Ls)]; 
set G 9.81
# El Centro record: 2688 points, dt = 0.02s original
timeSeries Path 2 -dt [expr 0.02*$Ts] -filePath ElCentro_acc.tcl -factor $G;
pattern UniformExcitation 2 1 -accel 2;

# damping via first eigenmode

set freq [expr [eigen -fullGenLapack 1]**0.5]
set dampRatio 0.02
rayleigh 0. 0. 0. [expr 2*$dampRatio/$freq]


wipeAnalysis;
constraints Plain;
numberer Plain;
system BandGeneral;
algorithm Linear
integrator Newmark 0.5 0.25 ;
analysis Transient;
# total time = 53.74s -> 2688 steps * 0.02s = 53.76s
# alternative analyze every 2nd step 1344Pts and dt*2 (=0.04s)
# is not necessery 
analyze 2688 $dt_analysis;

puts "Length scale Ls = $Ls"
puts "Time scale Ts = sqrt(Ls) = $Ts"
puts "Prototype earthquake dt = 0.02 s"
puts "Model ground motion dt = 0.02xTs = [expr 0.02*$Ts] s"
puts "analysis dt = $dt_analysis s"
puts "Mismatch (>1<) = [expr $dt_analysis/(0.02*$Ts)]"
puts "Time compression factor (prototype/model) = 1/Ts = [expr 1.0/$Ts]"
puts "El Centro has been shaking!"
wipe
