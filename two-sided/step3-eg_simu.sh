dirs=(
    "wild"
)

for group in "${dirs[@]}"
do
    data_dir="${group}/cool_outputs"
    python eg_simu.py "$data_dir" cr139.pkl "$group" 
    python eg_zoom.py "$data_dir" cr139.pkl "$group"
done


