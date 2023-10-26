cd docs
# Clear sphinx build tree
make clean
# Compile ThunderBoltz documentation
sphinx-build -M latex source build -D "root_doc=short_index"

cd ..
# Fix formating in C++ latex docs to prevent overrun tables
# and exclude API params
mkdir docs/build/latex/tables
python docs/reformat_latex.py runoff remove-members copy-table add_afil
cp docs/build/latex/thunderboltz.toc docs/tb_manual/tb_manual.toc

cd docs/tb_manual
pdflatex tb_manual
bibtex tb_manual
pdflatex tb_manual
pdflatex tb_manual
cd ../../
cp docs/tb_manual/tb_manual.pdf tb_manual.pdf

# Build API html docs
cd docs
make html && open build/html/index.html
# Build API latex docs
make latexpdf
cd ..
# Fix formating in API latex docs to prevent overrun tables.
python docs/reformat_latex.py runoff add_afil
cp docs/build/latex/thunderboltz.pdf pytb_manual.pdf
# Open latex doc
open pytb_manual.pdf
open tb_manual.pdf
