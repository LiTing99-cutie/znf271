
mkdir human mouse rhesus
pushd human
Rscript run.R ../../../evo_compare/h/h.cnt.bed
popd

pushd mouse
Rscript run.R ../../../evo_compare/m/m.cnt.bed
popd

pushd rhesus
Rscript run.R ../../../evo_compare/r/r.cnt.bed
popd