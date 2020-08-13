#!/bin/csh
set sourceroot = "$HOME/repositories/UCI-ACME-ECP/components"
foreach file ( ` find $sourceroot -name "*.[fF]*" ` )
sed -i '/\!\$omp/d' $file
sed -i '/\!\$OMP/d' $file
echo "Done $file"
end
