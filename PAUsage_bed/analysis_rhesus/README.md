# cp scripts
cp -r ../analysis_mouse/scripts ./

# config
./config.txt

# run
nohup bash -c "time bash run.sh" &>run.log &
source ./config.txt
bash scripts/Step_4_case.sh $md $frag $case