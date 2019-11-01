for file in  data/*.mtx; do 
	./lu $file
done

echo "done with files that should work"

for file in data/*.takes_forever; do 
	./lu $file
done
