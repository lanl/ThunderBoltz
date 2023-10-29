# If previously installed, uninstall
python -c "import thunderboltz" &> /dev/null
if [[ $? == 0 ]]; then
    pip uninstall -y thunderboltz
fi
# Update pip
pip install --user --upgrade pip
# Ensure upgraded build
pip install --user --upgrade build
# Ensure we are running from thunderboltz project root
FILE=pyproject.toml
if test -f "$FILE"; then
    if [[ $(sed -n 6p $FILE) == 'name = "thunderboltz"' ]]; then
        :
    else
        echo "Aborting - installation occured outside of setup directory"
    fi
else
    echo "Aborting - installation occured outside of setup directory"
fi
# Build and remove build tmp backend
pip install --user .
if test -d "build"; then
    rm -r build
fi
# Compile from source into bin
g++ -std=c++17 src/thunderboltz/cpp/DSMC0D.cpp -o bin/thunderboltz.bin -Wall -Werror -Wsign-compare
