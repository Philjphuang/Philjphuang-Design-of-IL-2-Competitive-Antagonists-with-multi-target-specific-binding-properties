
sed -i 's/CYX/CYS/g'  wt.pdb

sed -i 's/HIE/HIS/g'  wt.pdb
sed -i 's/HID/HIS/g'  wt.pdb
sed -i 's/HIP/HIS/g'  wt.pdb

fixbb.default.linuxgccrelease -in:file:s wt.pdb -out:path:all ./ -resfile resfile.in -overwrite -ex1:level 2 -ex2:level 2 -ex3:level 0 -ex4:level 0 -nstruct 1


