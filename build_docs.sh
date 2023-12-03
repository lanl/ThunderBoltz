cd docs
# Clear sphinx build tree
make clean
# Compile ThunderBoltz documentation
sphinx-build -M latex source build -D "root_doc=short_index"

cd ..
# Fix formating in C++ latex docs to prevent overrun tables
# and exclude API params
mkdir docs/build/latex/tables
python docs/reformat_latex.py runoff remove-members copy-table add_afil \
    --title "ThunderBoltz C++ Manual" --laur "23-31893"
cp docs/build/latex/thunderboltz.toc docs/cpp_manual/cpp_manual.toc

cd docs/cpp_manual
pdflatex cpp_manual
bibtex cpp_manual
pdflatex cpp_manual
pdflatex cpp_manual
cd ../../
cp docs/cpp_manual/cpp_manual.pdf cpp_manual.pdf

# Build API html docs
cd docs
sphinx-build -M html source build
open build/html/index.html
# Build API latex docs
make latexpdf
cd ..
# Fix formating in API latex docs to prevent overrun tables.
python docs/reformat_latex.py runoff add_afil \
    --title "ThunderBoltz API Manual" --laur "23-32356"
cp docs/build/latex/thunderboltz.pdf api_manual.pdf
# Open latex doc
open api_manual.pdf
open cpp_manual.pdf
