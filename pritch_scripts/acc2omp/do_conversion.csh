#!/bin/csh
set sourceroot = "$HOME/repositories/UCI-ACME-ECP/components/cam/src/physics/crm"
foreach file ( ` find $sourceroot -name "*.[fF]*" ` )
sed -i 's/\!\$acc parallel/\!\$omp parallel/g' $file
sed  -i 's/async(asyncid)//g' $file
echo "Done $file"
end
