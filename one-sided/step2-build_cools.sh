dirs=(
    "wild"
)

for group in "${dirs[@]}"
do
    output_dir="${group}/cool_outputs"
    mkdir -p "$output_dir"
    python build_cools.py "$group" "$output_dir"
done
