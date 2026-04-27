ls ./raw_data/*/*/*/*/*_genomic.fna | xargs -I @ grep "^>" @ | python3 ./extract_names.py > names.txt ;
#creates list of samples with chromosomes
#for f in ./raw_data/*/*/*/*/*_genomic.fna; do seqtk subseq $f names.txt > $f.filtered ; done
# extracts chromosomes from samples
for f in ./raw_data/*/*/*/*/*_genomic.fna; do
    seqtk subseq "$f" names.txt | \
    awk '/^>/ {
            if (match($0, /chromosome[[:space:]]+([0-9XYM]+)/, a))
                print ">chr" a[1];
            next
         }
         {print}' \
    > "$f.filtered"
done;

for file in ./raw_data/*/*/*/*/*_genomic.fna.filtered; do
  for i in {1..22}; do
    out="$file.chr$i"
    seqtk subseq "$file" <(echo "chr$i") > "$out"
    if ! grep -q "^>" "$out"; then
        rm -f "$out"
    fi;
  done;
done;
for file in hg002.?at.fasta; do
  base="${file%.fasta}"
  for i in {1..22}; do
    out="${base}_chr$i.fasta"
    seqtk subseq "$file" <(echo "chr$i") > "$out"
    if ! grep -q "^>" "$out"; then
        rm "$out"
    fi;
  done;
done;
