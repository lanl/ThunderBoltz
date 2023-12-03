for i in *; do pdftoppm -png ${i%%.*}.pdf > ${i%%.*}.png; done
