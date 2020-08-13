#!/bin/csh
set sourceroot = "$HOME/repositories/UCI-ACME-ECP/components/cam/src/physics/crm"
foreach file ( ` find $sourceroot -name "*.[fF]*" ` )
#replacements:
sed -i 's/\!\$acc parallel loop/\!\$omp parallel do/g' $file
sed -i 's/\!\$acc atomic update/\!\$omp atomic update/g' $file
sed -i 's/\!\$acc wait(asyncid)/\!\$omp taskwait/g' $file
#deletions:
sed  -i 's/async(asyncid)//g' $file
sed  -i 's/gang vector//g' $file
sed  -i 's/vector_length(128)//g' $file
sed  -i 's/num_gangs(numgangs)//g' $file
echo "Done $file"
end
